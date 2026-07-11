"""The part of the band misregistration that belongs to the scene, not the instrument.

Two corrections, and keeping them apart is the whole point of this module.

The **interior orientation** (:mod:`spp.refine.relative`) is a fixed angular offset per
band. It is a property of the optics and the detector, it is the same for every
acquisition this instrument ever makes, and it belongs in the calibration file. It is
estimated once and reused.

What remains is **not** that. On the reference acquisition the residual for the bands
furthest from the reference varies along the strip, and the size of that variation grows
with each band's *time separation* from the reference (correlation +0.82). A fixed
optical offset cannot do that. Attitude can: the bands are up to 0.47 s apart, the
attitude jitters by ~0.005° about a smooth fit, and that jitter does not cancel across
half a second.

That residual is a property of **this pass, on this day**. Writing it into the
calibration would produce a per-band "optical" constant that is really a snapshot of one
scene's pointing noise — wrong for every other acquisition, while looking, on this one,
like a triumph. So :mod:`spp.refine.relative` **refuses** to fit it.

This module corrects it anyway — in the *product*, never in the calibration. The
distinction is not pedantry. It is the difference between a co-registered image and a
corrupted instrument model, and the only thing keeping them apart is where the number is
written down.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np

from spp.geometry.geoloc_grid import GeolocationGrid
from spp.geometry.sensor_model import SensorModel
from spp.geometry.terrain import TerrainModel
from spp.refine.matcher import Match

logger = logging.getLogger(__name__)

DEFAULT_DEGREE = 3
"""Degree of the along-track polynomial.

Low, deliberately. The correction exists to follow a slow attitude drift, not to chase
every matching error: a high-degree fit would interpolate the noise, and the product
would be "co-registered" against nothing but its own measurement scatter.
"""

MIN_WINDOWS_PER_DEGREE = 8
"""Windows required per polynomial degree, before a fit is allowed at all."""


@dataclass(frozen=True)
class SceneCorrection:
    """A per-band ground displacement that varies along the strip.

    **This is not a calibration.** It describes one acquisition's pointing behaviour and
    must never be fed back into the instrument model. It is applied to this product's
    geolocation and reported as a scene property.

    Attributes
    ----------
    band:
        Band being corrected.
    east_coefficients, north_coefficients:
        Polynomial coefficients (highest power first, as :func:`numpy.polyval` expects)
        giving the ground displacement to **add** to the model's position, in metres, as
        a function of the normalised raster line.
    degree:
        Degree of the polynomial actually fitted.
    residual_before_px, residual_after_px:
        Root-mean-square band-to-band displacement, before and after.
    n_windows:
        Windows the fit rests on.
    fitted:
        Whether a correction was produced.
    reason:
        Why not, when ``fitted`` is ``False``.
    """

    band: str
    east_coefficients: np.ndarray
    north_coefficients: np.ndarray
    degree: int
    residual_before_px: float
    residual_after_px: float
    n_windows: int
    fitted: bool
    reason: str = ""

    def displacement(self, lines: np.ndarray, n_lines: int) -> tuple[np.ndarray, np.ndarray]:
        """Ground displacement to add, in metres, at the given raster lines."""
        normalised = np.asarray(lines, dtype=np.float64) / max(n_lines - 1, 1)
        return (
            np.polyval(self.east_coefficients, normalised),
            np.polyval(self.north_coefficients, normalised),
        )

    def apply(self, grid: GeolocationGrid, *, lat0: float) -> GeolocationGrid:
        """Bend a geolocation grid by this correction.

        The grid says where the model thinks each lattice node landed. This moves those
        nodes to where the imagery says they landed — which, once the band is resampled
        through the moved grid, is what makes it line up with the reference.
        """
        if not self.fitted:
            return grid

        lines = np.arange(grid.lon.shape[0], dtype=np.float64) * grid.step
        east, north = self.displacement(lines, grid.n_lines)

        metres_per_deg_lon = 111_320.0 * np.cos(np.radians(lat0))
        lon = grid.lon + (east / metres_per_deg_lon)[:, None]
        lat = grid.lat + (north / 110_574.0)[:, None]

        return GeolocationGrid(
            lon=lon, lat=lat, step=grid.step, n_lines=grid.n_lines, n_columns=grid.n_columns
        )


def estimate(
    model: SensorModel,
    band: str,
    reference_band: str,
    matches: list[Match],
    *,
    terrain: TerrainModel | None = None,
    gsd_m: float = 3.77,
    degree: int = DEFAULT_DEGREE,
) -> SceneCorrection:
    """Fit the along-track-varying part of a band's misregistration.

    Run **after** the constant interior-orientation offset has been applied, so that what
    this sees is only what a constant cannot explain.

    Parameters
    ----------
    model:
        Sensor model, already carrying whatever line-of-sight calibration is in force.
    band, reference_band:
        The band to correct, and the one to correct it onto.
    matches:
        Displacements measured against the reference band, in raster pixels.
    terrain:
        Elevation model — terrain parallax between two bands is real and must not end up
        inside this fit.
    gsd_m:
        Ground sample distance, for reporting in pixels.
    degree:
        Polynomial degree along-track.
    """
    from spp.refine.relative import _ground_offset  # shared measurement, one definition

    required = MIN_WINDOWS_PER_DEGREE * (degree + 1)
    if len(matches) < required:
        return _refused(
            band, len(matches), f"{len(matches)} windows is too few for degree {degree} "
            f"(needs {required}); a fit would interpolate the noise"
        )

    lines = np.array([m.line for m in matches], dtype=np.float64)
    columns = np.array([m.column for m in matches], dtype=np.float64)

    lat0 = float(np.nanmean(model.locate(band, lines, columns, terrain=terrain)[:, 1]))
    residual = _ground_offset(model, band, reference_band, matches, terrain, lat0=lat0)

    valid = np.isfinite(residual).all(axis=1)
    if valid.sum() < required:
        return _refused(band, int(valid.sum()), "too few windows located successfully")

    n_lines = int(np.nanmax(lines)) + 1
    normalised = lines[valid] / max(n_lines - 1, 1)

    # The correction is the negative of the residual: move the model to the imagery.
    east = np.polyfit(normalised, -residual[valid, 0], degree)
    north = np.polyfit(normalised, -residual[valid, 1], degree)

    before = float(np.sqrt(np.mean(np.sum(residual[valid] ** 2, axis=1))))
    corrected = residual[valid] + np.stack(
        [np.polyval(east, normalised), np.polyval(north, normalised)], axis=-1
    )
    after = float(np.sqrt(np.mean(np.sum(corrected**2, axis=1))))

    logger.info(
        "Band %s: scene correction (degree %d) over %d windows; band-to-band RMS "
        "%.2f -> %.2f px. This describes THIS acquisition's pointing, not the "
        "instrument — it is not written to the calibration.",
        band,
        degree,
        int(valid.sum()),
        before / gsd_m,
        after / gsd_m,
    )
    return SceneCorrection(
        band=band,
        east_coefficients=east,
        north_coefficients=north,
        degree=degree,
        residual_before_px=before / gsd_m,
        residual_after_px=after / gsd_m,
        n_windows=int(valid.sum()),
        fitted=True,
    )


def _refused(band: str, n_windows: int, reason: str) -> SceneCorrection:
    logger.warning("Band %s: no scene correction — %s", band, reason)
    return SceneCorrection(
        band=band,
        east_coefficients=np.zeros(1),
        north_coefficients=np.zeros(1),
        degree=0,
        residual_before_px=float("nan"),
        residual_after_px=float("nan"),
        n_windows=n_windows,
        fitted=False,
        reason=reason,
    )
