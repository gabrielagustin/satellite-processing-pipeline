"""The part of the geolocation error that grows along the strip.

The boresight correction (:mod:`spp.refine.absolute`) is a rigid rotation. It removes a
**constant** pointing bias, and by construction it can remove nothing else.

On the reference acquisition, what it leaves behind is not constant. Measured against the
reference orthoimage, the residual **grows from one end of the strip to the other**:

    lines      0 –  7,000     7.4 m
    lines  7,000 – 14,000    16.5 m
    lines 14,000 – 21,000    30.5 m

— a correlation of +0.68 with along-track position, and diagonal on the ground because
both components grow together. Some 30 m accumulated over 79 km of travel is an
along-track scale error of about 0.03%: a clock-rate error, or an attitude drift, or some
combination. From a single strip they are not separable, and it does not matter for the
product, because the *effect* is what is corrected.

Where this belongs, and where it does not
-----------------------------------------
In the **product**, never in the calibration. A drift is a property of *this pass on this
day* — the clock ran slightly fast during these sixteen seconds, or the platform rotated
slowly during them. Writing it into the instrument model would make every other
acquisition wrong while making this one look perfect, which is the same trap
:mod:`spp.refine.scene` exists to avoid, at a different scale.

Applied to every band identically
---------------------------------
Unlike the per-band corrections, this one is **common**: a clock or attitude drift moves
the whole product, not one band relative to another. Applying it per band would be fitting
the same physical effect several times over, and would silently absorb band-to-band
differences into what claims to be a global drift.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np

from spp.geometry.geoloc_grid import GeolocationGrid

logger = logging.getLogger(__name__)

DEFAULT_DEGREE = 1
"""Degree of the along-track model.

Linear by default, because a clock-rate error and a constant attitude drift are both
linear in time and that is what the evidence shows. A higher degree would start following
the measurement scatter, and the cross-validation in ``docs/l1c/limitations.md`` shows
there is nothing above linear that generalises.
"""

MIN_WINDOWS = 20
"""Below this, there is not enough evidence to separate a drift from noise."""


@dataclass(frozen=True)
class DriftCorrection:
    """A ground displacement that varies along the strip, common to every band.

    Attributes
    ----------
    east_coefficients, north_coefficients:
        Polynomial coefficients (highest power first) giving the displacement to **add**
        to the modelled position, in metres, as a function of the normalised raster line.
    residual_before_m, residual_after_m:
        Root-mean-square distance to the reference, before and after. The pair is the
        evidence; a correction that does not shrink it has corrected nothing.
    n_windows:
        Windows the fit rests on.
    fitted:
        Whether a correction was produced.
    reason:
        Why not, when ``fitted`` is ``False``.
    """

    east_coefficients: np.ndarray
    north_coefficients: np.ndarray
    residual_before_m: float
    residual_after_m: float
    n_windows: int
    fitted: bool
    reason: str = ""

    def apply(self, grid: GeolocationGrid, *, lat0: float) -> GeolocationGrid:
        """Bend a geolocation grid by this drift."""
        if not self.fitted:
            return grid

        lines = np.arange(grid.lon.shape[0], dtype=np.float64) * grid.step
        x = lines / max(grid.n_lines - 1, 1)
        east = np.polyval(self.east_coefficients, x)
        north = np.polyval(self.north_coefficients, x)

        metres_per_deg_lon = 111_320.0 * np.cos(np.radians(lat0))
        return GeolocationGrid(
            lon=grid.lon + (east / metres_per_deg_lon)[:, None],
            lat=grid.lat + (north / 110_574.0)[:, None],
            step=grid.step,
            n_lines=grid.n_lines,
            n_columns=grid.n_columns,
        )


def estimate(
    lines: np.ndarray,
    east_m: np.ndarray,
    north_m: np.ndarray,
    *,
    n_lines: int,
    degree: int = DEFAULT_DEGREE,
) -> DriftCorrection:
    """Fit the along-track drift from measured displacements.

    Parameters
    ----------
    lines:
        Raster line of each measurement.
    east_m, north_m:
        Measured ground displacement of the product **from** the reference, in metres.
        The correction is the negative of this.
    n_lines:
        Raster height, to normalise the along-track coordinate.
    degree:
        Degree of the along-track model.
    """
    finite = np.isfinite(east_m) & np.isfinite(north_m)
    lines, east_m, north_m = lines[finite], east_m[finite], north_m[finite]

    if lines.size < MIN_WINDOWS:
        return _refused(
            int(lines.size),
            f"{lines.size} windows is too few to separate a drift from scatter "
            f"(need {MIN_WINDOWS})",
        )

    # Reject outliers first. The field is measured window by window, and a window over a
    # cloud edge, a ship, or a stretch of featureless water locks onto noise and reports it
    # with the same confidence as one that nailed a coastline. Fitting a drift through
    # those is fitting the noise: on the reference acquisition the raw field has a
    # root-mean-square of 327 m, against a true drift of ~30 m.
    displacement = np.column_stack([east_m, north_m])
    centre = np.median(displacement, axis=0)
    deviation = np.abs(displacement - centre)
    mad = np.median(deviation, axis=0) * 1.4826
    inliers = np.all(deviation <= np.maximum(3.0 * mad, 5.0), axis=1)

    if inliers.sum() < MIN_WINDOWS:
        return _refused(
            int(inliers.sum()),
            f"only {int(inliers.sum())} windows survived outlier rejection "
            f"(of {lines.size})",
        )
    rejected = int(lines.size - inliers.sum())
    lines, east_m, north_m = lines[inliers], east_m[inliers], north_m[inliers]

    x = lines / max(n_lines - 1, 1)
    before = float(np.sqrt(np.mean(east_m**2 + north_m**2)))

    # Fit through BINNED MEDIANS, not through individual windows.
    #
    # The drift on the reference acquisition is ~25 m end to end, and the scatter of a
    # single window's measurement is 35-55 m. The signal is smaller than the noise on any
    # one window, so a least-squares fit through the raw points is dominated by the
    # scatter and comes out shallow: it recovered 18 m of a 25 m trend, leaving most of it
    # in the product.
    #
    # The trend is only visible in aggregate -- which is how it was found in the first
    # place, by looking at the median of each segment of the strip. So that is what is
    # fitted. Binning trades away resolution the data does not support anyway, and the
    # median inside each bin is robust to the outliers that survive the gate above.
    bins = _bin_medians(x, east_m, north_m, n_bins=max(degree + 2, 12))
    if bins is None:
        return _refused(int(lines.size), "too few populated along-track bins")
    bin_x, bin_east, bin_north = bins

    # The correction is the negative of the measured displacement: move the model to
    # where the reference says it should be.
    east = np.polyfit(bin_x, -bin_east, degree)
    north = np.polyfit(bin_x, -bin_north, degree)

    after = float(
        np.sqrt(
            np.mean(
                (east_m + np.polyval(east, x)) ** 2 + (north_m + np.polyval(north, x)) ** 2
            )
        )
    )
    trend_before = float(np.sqrt(np.mean(bin_east**2 + bin_north**2)))
    trend_after = float(
        np.sqrt(
            np.mean(
                (bin_east + np.polyval(east, bin_x)) ** 2
                + (bin_north + np.polyval(north, bin_x)) ** 2
            )
        )
    )
    logger.info(
        "      the along-track TREND (binned medians): %.0f m -> %.0f m", trend_before, trend_after
    )

    # Judge the fit on the TREND, not on the per-window root-mean-square. The scatter is
    # irreducible here -- it is scene, not geometry (see docs/l1c/limitations.md) -- and a
    # metric dominated by it would reject a correction that is doing exactly its job.
    if trend_after >= trend_before:
        return _refused(
            int(lines.size),
            f"the drift model does not reduce the along-track trend "
            f"({trend_before:.0f} m -> {trend_after:.0f} m)",
        )

    slope = float(np.hypot(east[-2] if degree >= 1 else 0.0, north[-2] if degree >= 1 else 0.0))
    logger.info(
        "Along-track drift: %.0f m accumulated end to end; residual %.0f m -> %.0f m "
        "over %d windows (%d rejected as outliers). This describes THIS pass (a clock "
        "rate, an attitude drift, or both -- not separable from one strip) and is applied "
        "to the product, not the calibration.",
        slope,
        before,
        after,
        lines.size,
        rejected,
    )
    return DriftCorrection(
        east_coefficients=east,
        north_coefficients=north,
        residual_before_m=before,
        residual_after_m=after,
        n_windows=int(lines.size),
        fitted=True,
    )


def _bin_medians(
    x: np.ndarray, east: np.ndarray, north: np.ndarray, *, n_bins: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """Median displacement in each along-track bin, dropping the empty ones."""
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    index = np.clip(np.digitize(x, edges) - 1, 0, n_bins - 1)

    centres, easts, norths = [], [], []
    for b in range(n_bins):
        inside = index == b
        if inside.sum() < 3:
            continue
        centres.append(0.5 * (edges[b] + edges[b + 1]))
        easts.append(float(np.median(east[inside])))
        norths.append(float(np.median(north[inside])))

    if len(centres) < 4:
        return None
    return np.array(centres), np.array(easts), np.array(norths)


def _refused(n_windows: int, reason: str) -> DriftCorrection:
    logger.warning("Along-track drift: not corrected — %s", reason)
    return DriftCorrection(
        east_coefficients=np.zeros(1),
        north_coefficients=np.zeros(1),
        residual_before_m=float("nan"),
        residual_after_m=float("nan"),
        n_windows=n_windows,
        fitted=False,
        reason=reason,
    )
