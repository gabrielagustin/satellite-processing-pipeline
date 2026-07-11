"""Measuring displacement between two images, robustly enough to trust.

Phase correlation on windows, with two deliberate choices that make the difference
between a measurement and a plausible-looking number.

**Match on gradients, not on radiance.** Two spectral bands look at the same ground
and disagree about it: vegetation is dark in blue and bright in the near infrared,
water is the reverse, and a coastline's contrast can invert outright between bands.
Correlating raw radiance across bands therefore biases towards whatever the two
happen to share, and over a scene with strong spectral structure it can lock onto the
wrong peak entirely. Edge *structure* is largely spectrally invariant — a shoreline is
an edge in every band, whatever its sign — so the gradient magnitude is the feature
that survives the comparison.

**Insist on texture, and say so when there is none.** Phase correlation always returns
a peak. Over featureless open water it returns a peak too, and that peak is noise. A
window is used only if it carries real structure, and a band with too few usable
windows yields *no estimate* rather than a confident average of nonsense.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class Match:
    """One window's measured displacement.

    Attributes
    ----------
    line, column:
        Centre of the window in the reference image, in pixels.
    dline, dcolumn:
        Displacement of the moving image relative to the reference, in pixels.
        Sub-pixel.
    peak:
        Height of the correlation peak, in ``[0, 1]``. A sharp, isolated peak is
        near 1; noise is near 0.
    """

    line: float
    column: float
    dline: float
    dcolumn: float
    peak: float


def gradient_magnitude(image: np.ndarray) -> np.ndarray:
    """Edge strength — the feature that survives a change of spectral band."""
    filled = np.nan_to_num(image.astype(np.float64))
    gy, gx = np.gradient(filled)
    return np.hypot(gy, gx)


def has_texture(image: np.ndarray, *, min_relative_contrast: float = 0.08) -> bool:
    """Whether a window carries enough structure to be worth correlating.

    Open water has none. Correlating it anyway produces a displacement — phase
    correlation always produces one — and that displacement is noise dressed as a
    measurement.
    """
    finite = image[np.isfinite(image)]
    if finite.size < image.size // 2:
        return False
    mean = float(finite.mean())
    if mean <= 0.0:
        return False
    return float(finite.std()) / mean >= min_relative_contrast


def phase_correlate(reference: np.ndarray, moving: np.ndarray) -> tuple[float, float, float]:
    """Sub-pixel displacement of ``moving`` relative to ``reference``.

    Implemented on NumPy's FFT rather than pulling in an image-processing
    dependency: the algorithm is twenty lines, and a dependency added for twenty
    lines is a dependency to keep working forever.

    Returns
    -------
    tuple
        ``(dline, dcolumn, peak)``. The displacement is in pixels, refined to
        sub-pixel by fitting a parabola to the correlation peak and its neighbours.

    Notes
    -----
    **Sign convention**, which is pinned by a test because getting it backwards does
    not produce a weaker correction — it produces one that *doubles* the error it was
    meant to remove:

        ``moving[line, column]`` shows the feature that
        ``reference[line + dline, column + dcolumn]`` shows.

    Equivalently: if the moving image is the reference translated by ``+d``, this
    returns ``-d``.
    """
    a = gradient_magnitude(reference)
    b = gradient_magnitude(moving)
    a -= a.mean()
    b -= b.mean()

    # A window taper: without it the FFT sees the tile edges as enormous step edges,
    # and their correlation swamps the scene's.
    taper = np.outer(np.hanning(a.shape[0]), np.hanning(a.shape[1]))
    fa = np.fft.fft2(a * taper)
    fb = np.fft.fft2(b * taper)

    cross = fa * np.conj(fb)
    cross /= np.abs(cross) + 1e-12  # whiten: phase only, amplitude discarded
    surface = np.fft.ifft2(cross).real

    peak_index = np.unravel_index(np.argmax(surface), surface.shape)
    shift = [_parabolic_peak(surface, peak_index, axis) for axis in (0, 1)]
    return shift[0], shift[1], float(surface.max())


def _parabolic_peak(surface: np.ndarray, peak: tuple[int, ...], axis: int) -> float:
    """Refine an integer correlation peak to sub-pixel, along one axis."""
    size = surface.shape[axis]
    index = peak[axis]

    def value(offset: int) -> float:
        position = list(peak)
        position[axis] = (index + offset) % size
        return float(surface[tuple(position)])

    before, centre, after = value(-1), value(0), value(+1)
    denominator = before - 2.0 * centre + after
    sub = 0.5 * (before - after) / denominator if abs(denominator) > 1e-12 else 0.0

    # The FFT's origin is at index 0; the second half of the axis is negative shift.
    signed = index - size if index > size // 2 else index
    return float(signed + sub)


def match_windows(
    reference: np.ndarray,
    moving: np.ndarray,
    *,
    window: int = 384,
    stride: int | None = None,
    min_peak: float = 0.02,
    min_relative_contrast: float = 0.08,
) -> list[Match]:
    """Measure displacement on every textured window of a pair of images.

    Parameters
    ----------
    reference, moving:
        Same-shape images. ``NaN`` marks invalid samples.
    window:
        Side of the correlation window, pixels.
    stride:
        Spacing between window centres. Defaults to the window size (no overlap).
    min_peak:
        Reject windows whose correlation peak is weaker than this — a flat
        correlation surface means the match did not lock onto anything.
    min_relative_contrast:
        Reject windows without this much structure (see :func:`has_texture`).

    Returns
    -------
    list of Match
        Possibly empty. An empty list is a legitimate answer over a featureless
        scene, and callers must treat it as "no estimate", never as "zero shift".
    """
    if reference.shape != moving.shape:
        raise ValueError(
            f"Images must be the same shape to be matched, got {reference.shape} "
            f"and {moving.shape}"
        )
    stride = stride or window

    matches: list[Match] = []
    for row in range(0, reference.shape[0] - window + 1, stride):
        for col in range(0, reference.shape[1] - window + 1, stride):
            ref_tile = reference[row : row + window, col : col + window]
            mov_tile = moving[row : row + window, col : col + window]

            if not has_texture(ref_tile, min_relative_contrast=min_relative_contrast):
                continue
            if np.isnan(mov_tile).mean() > 0.02:
                continue

            dline, dcolumn, peak = phase_correlate(ref_tile, mov_tile)
            if peak < min_peak:
                continue

            matches.append(
                Match(
                    line=row + window / 2,
                    column=col + window / 2,
                    dline=dline,
                    dcolumn=dcolumn,
                    peak=peak,
                )
            )
    return matches


def robust_median(matches: list[Match]) -> tuple[np.ndarray, np.ndarray, int]:
    """Consensus displacement, with outliers rejected by median absolute deviation.

    The median, not the mean: a handful of windows will lock onto a cloud edge, a
    sun-glint patch or a moving ship, and a mean would carry them into the answer.

    Returns
    -------
    tuple
        ``(displacement, scatter, n_used)`` — the consensus ``(dline, dcolumn)``,
        the robust scatter of the inliers, and how many windows survived.
    """
    if not matches:
        return np.array([np.nan, np.nan]), np.array([np.nan, np.nan]), 0

    shifts = np.array([[m.dline, m.dcolumn] for m in matches])
    centre = np.median(shifts, axis=0)

    deviation = np.abs(shifts - centre)
    mad = np.median(deviation, axis=0)
    # 3 scaled MADs ~ 3 sigma for a normal distribution; the 1.4826 makes the MAD
    # an unbiased estimator of the standard deviation.
    tolerance = np.maximum(3.0 * 1.4826 * mad, 0.5)  # never reject inside half a pixel

    inliers = np.all(deviation <= tolerance, axis=1)
    kept = shifts[inliers]
    if kept.size == 0:
        kept = shifts

    return np.median(kept, axis=0), kept.std(axis=0), int(len(kept))
