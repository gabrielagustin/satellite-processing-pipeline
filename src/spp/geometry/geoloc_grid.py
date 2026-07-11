"""The geolocation grid: the sensor model, sampled onto a lattice.

The forward model maps sensor coordinates to the ground. Resampling needs the
*inverse* — for each output pixel, which detector sample saw it? — and for a
pushbroom sensor that inverse has no closed form: each output pixel would require
solving for the line time that observed it, some 10⁹ solves across eight bands.

The standard escape, and the one used here, is to evaluate the forward model on a
**subsampled lattice** and hand that to GDAL as a *geolocation array*. GDAL builds
the inverse transformer from it, interpolating between the nodes. This is the same
machinery that grids swath products from instruments like MODIS and Sentinel-3.

Why the spacing is 8 and not 64
-------------------------------
The tempting choice is a coarse lattice: the orbit and attitude are smooth, so a
node every 64 pixels would capture them perfectly and cost a sixty-fourth of the
work.

That reasoning is wrong, and quietly so. The geolocation field has a smooth
component (orbit, attitude) **plus a terrain component that varies at the elevation
model's own resolution**. Sampling every 64 pixels (~240 m) low-passes the terrain
correction — it would silently undo part of the orthorectification, leaving a
product that looks orthorectified and is not.

So the node spacing tracks the **DEM**, not the smoothness of the orbit: 8 pixels
is ~30 m on the ground, matching the elevation model. The cost is trivial (a few
tens of megabytes per band) and the alternative is a subtly broken product.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np

from spp.geometry.sensor_model import SensorModel
from spp.geometry.terrain import TerrainModel

logger = logging.getLogger(__name__)

DEFAULT_STEP = 8
"""Lattice spacing, in detector samples. Chosen to match the DEM — see above."""


@dataclass(frozen=True)
class GeolocationGrid:
    """Ground coordinates of a subsampled lattice over one band's raster.

    Attributes
    ----------
    lon, lat:
        Longitude and latitude of each lattice node, degrees, shape
        ``(rows, cols)``. **``float64``, and not negotiable**: ``float32``
        resolves longitude to only ~0.7 m at mid latitudes, which would inject
        noise of the same order as the terrain correction being applied.
    step:
        Lattice spacing in detector samples, along both axes.
    n_lines, n_columns:
        Dimensions of the raster the lattice covers.
    """

    lon: np.ndarray
    lat: np.ndarray
    step: int
    n_lines: int
    n_columns: int

    def __post_init__(self) -> None:
        if self.lon.dtype != np.float64 or self.lat.dtype != np.float64:
            raise ValueError(
                "Geolocation arrays must be float64; float32 resolves longitude to "
                "only ~0.7 m, comparable to the corrections being applied"
            )

    @property
    def as_array(self) -> np.ndarray:
        """Stacked ``(2, rows, cols)``, the layout GDAL's warper expects."""
        return np.stack([self.lon, self.lat])

    @property
    def bounds(self) -> tuple[float, float, float, float]:
        """``(min_lon, min_lat, max_lon, max_lat)`` of the valid nodes."""
        return (
            float(np.nanmin(self.lon)),
            float(np.nanmin(self.lat)),
            float(np.nanmax(self.lon)),
            float(np.nanmax(self.lat)),
        )

    @property
    def valid_fraction(self) -> float:
        """Fraction of nodes whose ray reached the ground."""
        return float(np.isfinite(self.lon).mean())


def build(
    model: SensorModel,
    band: str,
    *,
    n_lines: int,
    n_columns: int,
    terrain: TerrainModel | None = None,
    step: int = DEFAULT_STEP,
) -> GeolocationGrid:
    """Evaluate the sensor model on a lattice over one band's raster.

    Parameters
    ----------
    model:
        The forward sensor model.
    band:
        Band name. Each band gets its **own** grid: its own line times, its own
        view angles. That is what co-registers them — they are not offsets applied
        to a shared grid, they are independently located.
    n_lines, n_columns:
        Raster dimensions of the band.
    terrain:
        Elevation model to intersect. Without one this is georeferencing on the
        ellipsoid, and the lattice spacing (§ module docstring) is over-fine.
    step:
        Lattice spacing in detector samples.

    Returns
    -------
    GeolocationGrid
    """
    # The lattice must reach the far edge of the raster, not stop short of it:
    # GDAL extrapolates beyond the last node, and extrapolated geolocation at the
    # image border is exactly where a swath product goes wrong.
    lines = np.arange(0, n_lines + step, step, dtype=np.float64)
    columns = np.arange(0, n_columns + step, step, dtype=np.float64)
    mesh_lines, mesh_columns = np.meshgrid(lines, columns, indexing="ij")

    ground = model.locate(
        band,
        mesh_lines.ravel(),
        mesh_columns.ravel(),
        terrain=terrain,
    )
    lon = ground[:, 0].reshape(mesh_lines.shape).astype(np.float64)
    lat = ground[:, 1].reshape(mesh_lines.shape).astype(np.float64)

    grid = GeolocationGrid(
        lon=lon, lat=lat, step=step, n_lines=n_lines, n_columns=n_columns
    )
    if grid.valid_fraction < 1.0:
        logger.warning(
            "Band %s: %.2f%% of geolocation nodes have no ground intersection",
            band,
            100.0 * (1.0 - grid.valid_fraction),
        )
    return grid


def interpolation_error_px(
    model: SensorModel,
    grid: GeolocationGrid,
    band: str,
    *,
    terrain: TerrainModel | None = None,
    n_samples: int = 500,
    gsd_m: float = 3.77,
    seed: int = 0,
) -> float:
    """Largest error the lattice's interpolation introduces, in pixels.

    The geolocation grid is an approximation: GDAL interpolates between its nodes,
    and the truth is the model evaluated at full resolution. This measures the gap
    at random points *between* nodes — where the error is, by construction, largest.

    It is reported rather than assumed. A lattice too coarse for the terrain would
    pass every other check in this pipeline and still produce a product whose
    orthorectification has been quietly smoothed away (§ module docstring), and this
    is the measurement that would catch it.

    Returns
    -------
    float
        Maximum discrepancy between the interpolated lattice and the model, in
        ground pixels.
    """
    rng = np.random.default_rng(seed)
    lines = rng.uniform(0, grid.n_lines - 1, n_samples)
    columns = rng.uniform(0, grid.n_columns - 1, n_samples)

    truth = model.locate(band, lines, columns, terrain=terrain)

    # Bilinear interpolation of the lattice, exactly as the warper will do it.
    fl, fc = lines / grid.step, columns / grid.step
    i0 = np.clip(np.floor(fl).astype(int), 0, grid.lon.shape[0] - 2)
    j0 = np.clip(np.floor(fc).astype(int), 0, grid.lon.shape[1] - 2)
    dy, dx = fl - i0, fc - j0

    def sample(field: np.ndarray) -> np.ndarray:
        top = field[i0, j0] * (1 - dx) + field[i0, j0 + 1] * dx
        bottom = field[i0 + 1, j0] * (1 - dx) + field[i0 + 1, j0 + 1] * dx
        return top * (1 - dy) + bottom * dy

    lat0 = float(np.nanmean(truth[:, 1]))
    metres_per_deg_lat = 110_574.0
    metres_per_deg_lon = 111_320.0 * np.cos(np.radians(lat0))

    dlon = (sample(grid.lon) - truth[:, 0]) * metres_per_deg_lon
    dlat = (sample(grid.lat) - truth[:, 1]) * metres_per_deg_lat
    error_m = np.hypot(dlon, dlat)

    return float(np.nanmax(error_m) / gsd_m)
