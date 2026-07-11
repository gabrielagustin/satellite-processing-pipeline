"""Estimating the interior orientation the instrument never shipped.

The Calibration Parameter File's geometric section is a placeholder: an identity
boresight, and line-of-sight coefficients identical for all eight bands. That cannot
be physically true — the bands sit on different detector rows and demonstrably look in
different directions — so the sensor model reconstructs the interior orientation from
the pinhole geometry, and the calibration's line-of-sight term enters as an additive
angular correction that is **zero as delivered**.

Zero is not the truth. Measured against the imagery, the model places each band wrong
by a *constant* angular offset — about 8 detector pixels for the extreme band, which is
an ordinary amount of optical distortion and detector placement for an instrument whose
geometric calibration was never filled in. That offset is exactly what the empty term
is for.

So the level does not work around the missing calibration; it **estimates** it. The
band-to-band displacement measured in the imagery is inverted, through the sensor model,
into the per-band angular offsets that would produce it — and those offsets are written
out in the same schema the input file left unpopulated.

Why this is legitimate, and when it would not be
------------------------------------------------
Feeding a residual back into a model can correct a missing calibration or **conceal a
bug**, and nothing about the arithmetic distinguishes the two. What distinguishes them
is the residual's *shape*: an interior-orientation error is a fixed angular offset per
band, constant across the strip and the swath, while a modelling error varies with
position. On the reference acquisition the residual's correlation with along-track and
cross-track position is −0.09, −0.05, +0.01 and −0.17 — that is, none.

The check is not decoration. :func:`estimate` measures the same correlations and
**refuses to fit** a band whose residual is position-dependent, because for such a band
the fit would be papering over something else.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np

from spp.geometry.camera import Camera
from spp.geometry.sensor_model import SensorModel
from spp.geometry.terrain import TerrainModel
from spp.refine.matcher import Match, robust_median

logger = logging.getLogger(__name__)

PERTURBATION_RAD = 1e-5
"""Step used to probe the model's sensitivity to a line-of-sight offset.

Small enough that the response is linear (a 10 µrad offset moves the ground point by
about 4 m, a fraction of a pixel), large enough to stay well clear of floating-point
noise in the geodetic conversions.
"""

MAX_POSITION_CORRELATION = 0.4
"""Above this, the residual varies with position and is not an interior-orientation error."""


@dataclass(frozen=True)
class BandCorrection:
    """The line-of-sight offset estimated for one band.

    Attributes
    ----------
    band:
        Band name.
    along_rad, across_rad:
        Additive angular corrections, radians. These populate the calibration term
        the input file ships empty.
    residual_before_px, residual_after_px:
        The **systematic** model-versus-imagery displacement — the magnitude of the
        median residual vector — before and after the correction. This is what a
        per-band angular offset can remove, and driving it to zero is what the
        correction does. The pair is the evidence: a correction that does not shrink it
        has not corrected anything.
    scatter_px:
        The window-to-window **scatter** of the residual, which a constant per-band
        offset cannot touch. It is the noise floor of the co-registration — matching
        noise, attitude jitter between the bands' 0.78 s separation, elevation-model
        error — and reporting it separately from the systematic part is the difference
        between an honest figure and a flattering one. A correction that drove the
        systematic term to zero while leaving 3 px of scatter has not delivered 0 px of
        co-registration.
    n_windows:
        Number of image windows the estimate rests on.
    position_correlation:
        Largest absolute correlation between the residual and the sample's position
        in the raster. Near zero means a fixed angular offset — an interior-orientation
        error. Large means something position-dependent, which this correction must not
        be used to hide.
    fitted:
        Whether a correction was accepted. ``False`` when there was too little
        evidence, or when the residual was position-dependent.
    reason:
        Why, when ``fitted`` is ``False``.
    """

    band: str
    along_rad: float
    across_rad: float
    residual_before_px: float
    residual_after_px: float
    scatter_px: float
    n_windows: int
    position_correlation: float
    fitted: bool
    reason: str = ""

    @property
    def detector_px(self) -> float:
        """The offset expressed in detector pixels, which is how an optician reads it."""
        return float(np.hypot(self.along_rad, self.across_rad))


def _ground_offset(
    model: SensorModel,
    band: str,
    reference_band: str,
    matches: list[Match],
    terrain: TerrainModel | None,
    *,
    lat0: float,
) -> np.ndarray:
    """Where the model puts each matched sample, minus where the imagery says it is.

    The imagery's answer is the reference band's own geolocation, evaluated at the
    matched position: this is a **relative** measurement, correcting the bands onto the
    reference rather than onto the ground. Absolute accuracy is a separate problem with
    a separate reference (an external orthoimage), and conflating the two would let an
    absolute bias leak into the interior orientation.
    """
    lines = np.array([m.line for m in matches], dtype=np.float64)
    columns = np.array([m.column for m in matches], dtype=np.float64)
    dline = np.array([m.dline for m in matches], dtype=np.float64)
    dcolumn = np.array([m.dcolumn for m in matches], dtype=np.float64)

    modelled = model.locate(band, lines, columns, terrain=terrain)

    # The matcher's sign convention (pinned by a test in `matcher`): the feature seen
    # at (line, column) in the moving band is the one at (line + dline, column +
    # dcolumn) in the reference band. Reading this backwards does not weaken the
    # correction, it *doubles* the error — which is what happened the first time, and
    # what `residual_after_px` caught.
    reference = model.locate(
        reference_band, lines + dline, columns + dcolumn, terrain=terrain
    )
    return _to_metres(modelled[:, :2] - reference[:, :2], lat0)


def _to_metres(delta_degrees: np.ndarray, lat0: float) -> np.ndarray:
    """Longitude/latitude differences to metres on a local plane."""
    metres = np.empty_like(delta_degrees)
    metres[:, 0] = delta_degrees[:, 0] * 111_320.0 * np.cos(np.radians(lat0))
    metres[:, 1] = delta_degrees[:, 1] * 110_574.0
    return metres


def estimate(
    model: SensorModel,
    band: str,
    reference_band: str,
    matches: list[Match],
    *,
    terrain: TerrainModel | None = None,
    gsd_m: float = 3.77,
    min_windows: int = 8,
) -> BandCorrection:
    """Invert a band's measured displacement into a line-of-sight offset.

    The relationship between an angular offset and the ground displacement it produces
    is derived **numerically**, by perturbing the model and watching where the ground
    point goes. That is deliberate: the alternative is to reason through the chain of
    sign conventions from detector axis to camera frame to body to orbital frame to
    Earth-fixed, and a single sign error there produces a correction that confidently
    doubles the very error it was meant to remove.

    Parameters
    ----------
    model:
        The sensor model, with the band's current (usually zero) line-of-sight term.
    band:
        Band to correct.
    reference_band:
        The band everything is corrected *onto*. This is a relative measurement: it
        aligns the bands with each other, not with the ground.
    matches:
        Displacements measured against the reference band, in the raster's pixels.
    terrain:
        Elevation model. Terrain parallax between two bands is real (they are a
        short-baseline stereo pair), so leaving it out would push relief into the
        estimated calibration.
    gsd_m:
        Ground sample distance, for reporting residuals in pixels.
    min_windows:
        Refuse to fit on fewer windows than this.

    Returns
    -------
    BandCorrection
        With ``fitted=False`` and a reason, when the evidence does not support a fit.
    """
    if len(matches) < min_windows:
        return _refused(band, len(matches), f"only {len(matches)} usable windows")

    lat0 = float(
        np.nanmean(
            model.locate(
                band,
                np.array([m.line for m in matches]),
                np.array([m.column for m in matches]),
                terrain=terrain,
            )[:, 1]
        )
    )

    residual = _ground_offset(model, band, reference_band, matches, terrain, lat0=lat0)

    # The systematic part -- the magnitude of the median vector -- is what a constant
    # angular offset can remove. The scatter about it is what it cannot.
    before = float(np.linalg.norm(np.median(residual, axis=0)))
    scatter = float(np.median(np.abs(residual - np.median(residual, axis=0))) * 1.4826)

    # Is this a fixed angular offset, or something that varies across the raster?
    lines = np.array([m.line for m in matches])
    columns = np.array([m.column for m in matches])
    correlation = max(
        abs(np.corrcoef(axis, residual[:, component])[0, 1])
        for axis in (lines, columns)
        for component in (0, 1)
    )
    if correlation > MAX_POSITION_CORRELATION:
        return _refused(
            band,
            len(matches),
            f"residual varies with position (correlation {correlation:.2f}); this is "
            "not an interior-orientation offset, and fitting it would hide whatever "
            "it is",
            before_px=before / gsd_m,
            correlation=correlation,
        )

    # The model's sensitivity to each angular offset, measured rather than derived.
    jacobian = np.empty((2, 2))
    for index, axis in enumerate(("along", "across")):
        perturbed = SensorModel(
            _with_offset(model.camera, band, axis, PERTURBATION_RAD),
            model.ephemeris,
            model.attitude,
            model.timing,
            convention=model.convention,
            rate_aware=model.rate_aware,
        )
        moved = _ground_offset(perturbed, band, reference_band, matches, terrain, lat0=lat0)
        jacobian[:, index] = np.median(moved - residual, axis=0) / PERTURBATION_RAD

    if abs(np.linalg.det(jacobian)) < 1e-9:
        return _refused(
            band, len(matches), "the model is insensitive to this band's line of sight",
            before_px=before / gsd_m, correlation=correlation,
        )

    # Solve for the offset that drives the median residual to zero.
    along, across = np.linalg.solve(jacobian, -np.median(residual, axis=0))

    corrected = SensorModel(
        _with_offset(model.camera, band, "along", along, across=across),
        model.ephemeris,
        model.attitude,
        model.timing,
        convention=model.convention,
        rate_aware=model.rate_aware,
    )
    corrected_residual = _ground_offset(
        corrected, band, reference_band, matches, terrain, lat0=lat0
    )
    after = float(np.linalg.norm(np.median(corrected_residual, axis=0)))

    logger.info(
        "Band %s: line-of-sight offset (%.1f, %.1f) urad; systematic residual "
        "%.2f -> %.2f px, scatter %.2f px, over %d windows",
        band,
        along * 1e6,
        across * 1e6,
        before / gsd_m,
        after / gsd_m,
        scatter / gsd_m,
        len(matches),
    )
    return BandCorrection(
        band=band,
        along_rad=float(along),
        across_rad=float(across),
        residual_before_px=before / gsd_m,
        residual_after_px=after / gsd_m,
        scatter_px=scatter / gsd_m,
        n_windows=len(matches),
        position_correlation=float(correlation),
        fitted=True,
    )


def _with_offset(
    camera: Camera, band: str, axis: str, value: float, *, across: float | None = None
) -> Camera:
    """A copy of the camera with this band's line-of-sight term shifted.

    The baseline is the calibration **in force** (:meth:`Camera.effective_los`), not the
    coefficients as stored. The delivered file stores an unpopulated sentinel whose
    leading term is ``1.0``; read literally that is a one-radian — 57 degree —
    line-of-sight offset. Inheriting it as the baseline for a correction produces a
    perfectly self-consistent answer that is catastrophically wrong, which is what the
    ``residual_after_px`` check exists to catch.
    """
    along_coefficients, across_coefficients = camera.effective_los(band)
    along_base = along_coefficients[0] if along_coefficients else 0.0
    across_base = across_coefficients[0] if across_coefficients else 0.0

    if across is not None:
        along_new, across_new = along_base + value, across_base + across
    elif axis == "along":
        along_new, across_new = along_base + value, across_base
    else:
        along_new, across_new = along_base, across_base + value

    return camera.with_los({band: ((along_new,), (across_new,))})


def _refused(
    band: str,
    n_windows: int,
    reason: str,
    *,
    before_px: float = float("nan"),
    correlation: float = float("nan"),
) -> BandCorrection:
    """No estimate, and why. Never a silent zero — that would read as 'no error'."""
    logger.warning("Band %s: no line-of-sight correction estimated — %s", band, reason)
    return BandCorrection(
        band=band,
        along_rad=0.0,
        across_rad=0.0,
        residual_before_px=before_px,
        residual_after_px=before_px,
        scatter_px=float("nan"),
        n_windows=n_windows,
        position_correlation=correlation,
        fitted=False,
        reason=reason,
    )
