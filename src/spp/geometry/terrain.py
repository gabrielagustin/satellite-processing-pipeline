"""Terrain: where the viewing ray actually meets the ground.

Georeferencing intersects the ray with an ellipsoid. **Ortho**rectification
intersects it with the *terrain*, which is what removes relief displacement — and,
for this instrument, what holds the spectral bands registered to each other over
relief (see ``docs/l1c/spec.md`` §3.4: the bands are a short-baseline stereo pair,
so terrain displaces them relative to one another far more than it displaces the
product absolutely).

The height datum is the trap
----------------------------
A digital elevation model (DEM) reports **orthometric** height — height above the
geoid, i.e. above mean sea level. The platform ephemeris is referenced to the
**ellipsoid**. In many regions those differ by tens of metres, and mixing them is a
silent, systematic error that no downstream check will attribute to its cause.

Worse, the conversion fails *quietly*. Asked to convert orthometric to ellipsoidal
height without the geoid grid installed, PROJ does not raise — it performs a
"ballpark" transformation that returns the height unchanged, i.e. an undulation of
exactly zero. Over the reference acquisition the true undulation is about −31 m, so
the silent failure costs 31 m of height everywhere. This module refuses to accept
that: it demands the real grid and fails loudly if it cannot have it.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np

from spp.geometry.frames import (
    ecef_to_geodetic,
    ray_ellipsoid_intersection,
)

logger = logging.getLogger(__name__)

WGS84_3D = "EPSG:4979"
"""WGS84 with ellipsoidal height."""

EGM2008_HEIGHT = "EPSG:9518"
"""WGS84 horizontal with EGM2008 orthometric height — the Copernicus DEM's datum."""


def check_undulation(undulation: np.ndarray, source_crs: str = EGM2008_HEIGHT) -> None:
    """Refuse a geoid model that is secretly not a geoid model.

    Asked to convert orthometric height to ellipsoidal height without the geoid
    grid installed, PROJ does **not** raise. It performs a "ballpark vertical
    transformation", which returns the height unchanged — an undulation of exactly
    zero, everywhere. Nothing downstream would notice: the heights are finite, the
    rays converge, the product looks fine, and every elevation is wrong by the local
    undulation (about −31 m over the reference acquisition).

    A silent, uniform, tens-of-metres height error is precisely the kind of defect
    that survives to production, so it is made loud here.

    Raises
    ------
    RuntimeError
        If the undulation is non-finite, or identically zero over the whole area.
    """
    if not np.isfinite(undulation).all():
        raise RuntimeError(
            f"Geoid model returned non-finite undulations; the vertical datum grid "
            f"for {source_crs} is unusable"
        )
    if np.allclose(undulation, 0.0, atol=1e-6):
        raise RuntimeError(
            "Geoid undulation is identically zero over the whole footprint. PROJ has "
            "fallen back to a 'ballpark' vertical transformation, which returns the "
            "height unchanged, because the geoid grid is unavailable. Elevations "
            "would be wrong by the local undulation (tens of metres) with no other "
            "symptom. Enable PROJ's network grid access "
            "(`pyproj.network.set_network_enabled(True)`) or install the grid with "
            "`projsync`."
        )


class TerrainModel(ABC):
    """A surface for viewing rays to intersect."""

    @abstractmethod
    def height(self, lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
        """Ellipsoidal height of the surface, metres. Shape follows the inputs."""

    def intersect(
        self,
        origins: np.ndarray,
        directions: np.ndarray,
        *,
        tolerance_m: float = 0.4,
        max_iterations: int = 8,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Intersect rays with the terrain, iteratively.

        The ray's ground point depends on the terrain height there, and the height
        depends on the ground point — so the two are solved together:

        1. Intersect the ellipsoid at the current height estimate.
        2. Sample the terrain at that point.
        3. Re-intersect at the new height, and repeat.

        Convergence is fast for a narrow-field instrument: the ellipsoid guess is
        already within a metre or two of the answer, so two iterations usually
        suffice. It is not assumed, though — it is measured and reported.

        Parameters
        ----------
        origins, directions:
            ECEF rays, shape ``(n, 3)``.
        tolerance_m:
            Stop when the ground point moves less than this between iterations.
            The default is a tenth of a pixel at the reference acquisition's
            ground sample distance.
        max_iterations:
            Give up after this many. Non-convergence is reported, not hidden.

        Returns
        -------
        tuple of numpy.ndarray
            ``(ground_ecef, converged)`` — the ECEF intersection points, shape
            ``(n, 3)``, and a boolean mask of which rays converged. Rays that
            missed the Earth are ``NaN`` and marked not converged.
        """
        point = ray_ellipsoid_intersection(origins, directions, height=0.0)
        converged = np.zeros(len(np.atleast_2d(point)), dtype=bool)

        for _ in range(max_iterations):
            geodetic = ecef_to_geodetic(point)
            valid = np.isfinite(geodetic).all(axis=-1)
            if not valid.any():
                break

            height = np.zeros(len(geodetic))
            height[valid] = self.height(geodetic[valid, 0], geodetic[valid, 1])

            new_point = ray_ellipsoid_intersection(origins, directions, height=height)
            moved = np.linalg.norm(new_point - point, axis=-1)
            point = new_point

            converged = np.where(valid, moved < tolerance_m, False)
            if converged[valid].all():
                break
        else:
            n_bad = int((~converged).sum())
            if n_bad:
                logger.warning(
                    "Terrain intersection: %d ray(s) did not converge within %d "
                    "iterations; their ground points are approximate",
                    n_bad,
                    max_iterations,
                )

        return point, converged


class EllipsoidTerrain(TerrainModel):
    """A surface of constant height above the ellipsoid.

    This is *georeferencing*, not orthorectification. It is exact over water at
    height zero, correct nowhere else, and useful as the reference against which
    the terrain correction is measured.
    """

    def __init__(self, height_m: float = 0.0) -> None:
        self.height_m = float(height_m)

    def height(self, lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
        return np.full(np.shape(lon), self.height_m, dtype=np.float64)


@dataclass(frozen=True)
class GeoidModel:
    """Geoid undulation: the offset between orthometric and ellipsoidal height.

    Evaluated on a coarse cached grid and interpolated. The geoid is smooth — it
    varies by centimetres per kilometre — so a coarse grid costs nothing in
    accuracy and saves a per-pixel call into PROJ.

    Attributes
    ----------
    lons, lats:
        Grid coordinates, degrees.
    undulation:
        Undulation on that grid, metres. ``ellipsoidal = orthometric + undulation``.
    """

    lons: np.ndarray
    lats: np.ndarray
    undulation: np.ndarray

    @classmethod
    def build(
        cls,
        bounds: tuple[float, float, float, float],
        *,
        spacing_deg: float = 0.05,
        source_crs: str = EGM2008_HEIGHT,
    ) -> GeoidModel:
        """Sample the geoid over a bounding box.

        Parameters
        ----------
        bounds:
            ``(min_lon, min_lat, max_lon, max_lat)``, degrees.
        spacing_deg:
            Grid spacing. At 0.05° (~5 km) the interpolation error is
            millimetric, the geoid gradient being of order 0.1 m per km.
        source_crs:
            The DEM's vertical datum.

        Raises
        ------
        RuntimeError
            If the geoid grid is unavailable, rather than silently returning a
            zero undulation — see the module docstring.
        """
        import pyproj
        from pyproj import Transformer

        pyproj.network.set_network_enabled(True)

        min_lon, min_lat, max_lon, max_lat = bounds
        lons = np.arange(min_lon - spacing_deg, max_lon + 2 * spacing_deg, spacing_deg)
        lats = np.arange(min_lat - spacing_deg, max_lat + 2 * spacing_deg, spacing_deg)
        mesh_lon, mesh_lat = np.meshgrid(lons, lats)

        transformer = Transformer.from_crs(source_crs, WGS84_3D, always_xy=True)
        _, _, ellipsoidal = transformer.transform(
            mesh_lon.ravel(), mesh_lat.ravel(), np.zeros(mesh_lon.size)
        )
        undulation = np.asarray(ellipsoidal, dtype=np.float64).reshape(mesh_lon.shape)
        check_undulation(undulation, source_crs)

        logger.info(
            "Geoid: undulation over the footprint spans %.2f to %.2f m",
            undulation.min(),
            undulation.max(),
        )
        return cls(lons=lons, lats=lats, undulation=undulation)

    def __call__(self, lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
        """Undulation at the given coordinates, metres, bilinearly interpolated."""
        return _bilinear(self.undulation, self.lons, self.lats, lon, lat)


class DEMTerrain(TerrainModel):
    """A raster digital elevation model, on the ellipsoid.

    The DEM's orthometric heights are converted to ellipsoidal heights with the
    geoid, so the surface lives in the same datum as the platform ephemeris.

    Parameters
    ----------
    elevation:
        Orthometric heights, shape ``(rows, cols)``, metres. Row 0 is the
        northernmost.
    lons, lats:
        Cell-centre coordinates of the grid, degrees. ``lats`` descending.
    geoid:
        Undulation model covering the same area.
    voids:
        Boolean mask of cells with no valid elevation, filled with sea level.
    """

    def __init__(
        self,
        elevation: np.ndarray,
        lons: np.ndarray,
        lats: np.ndarray,
        geoid: GeoidModel,
        *,
        voids: np.ndarray | None = None,
    ) -> None:
        self.elevation = np.asarray(elevation, dtype=np.float64)
        self.lons = np.asarray(lons, dtype=np.float64)
        self.lats = np.asarray(lats, dtype=np.float64)
        self.geoid = geoid
        self.voids = voids if voids is not None else np.zeros(self.elevation.shape, bool)

    @property
    def void_fraction(self) -> float:
        """Fraction of the DEM window with no valid elevation."""
        return float(self.voids.mean())

    def height(self, lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
        """Ellipsoidal height of the terrain: orthometric plus undulation."""
        orthometric = _bilinear(self.elevation, self.lons, self.lats, lon, lat)
        return orthometric + self.geoid(lon, lat)

    @classmethod
    def from_raster(
        cls,
        path: str,
        bounds: tuple[float, float, float, float],
        *,
        margin_deg: float = 0.05,
        geoid: GeoidModel | None = None,
    ) -> DEMTerrain:
        """Load the DEM window covering ``bounds`` (plus a margin) into memory.

        The margin matters: terrain correction *moves* the ground point, so the
        ray can land outside the footprint computed on the ellipsoid. Loading
        exactly the ellipsoid footprint would sample off the edge of the array
        precisely where the correction is largest.

        Parameters
        ----------
        path:
            Any raster GDAL can open, including a remote Cloud-Optimized GeoTIFF
            (``/vsicurl/...``) or a virtual mosaic (``.vrt``).
        bounds:
            ``(min_lon, min_lat, max_lon, max_lat)`` of the area needed.
        margin_deg:
            Extra margin loaded around ``bounds``.
        geoid:
            Undulation model. Built over the same area if not supplied.
        """
        import rasterio
        from rasterio.windows import from_bounds

        min_lon, min_lat, max_lon, max_lat = bounds
        box = (
            min_lon - margin_deg,
            min_lat - margin_deg,
            max_lon + margin_deg,
            max_lat + margin_deg,
        )

        with rasterio.open(path) as src:
            window = from_bounds(*box, transform=src.transform).round_offsets().round_lengths()
            elevation = src.read(1, window=window, boundless=True, fill_value=np.nan).astype(
                np.float64
            )
            transform = src.window_transform(window)
            nodata = src.nodata

        voids = ~np.isfinite(elevation)
        if nodata is not None:
            voids |= elevation == nodata

        # Sea level over voids. The Copernicus DEM is void over open water, and
        # the physically correct surface there is the geoid — height zero
        # *orthometric*, which is NOT height zero ellipsoidal. Filling with the
        # ellipsoid instead would put the sea surface tens of metres off.
        elevation[voids] = 0.0

        rows, cols = elevation.shape
        lons = transform.c + transform.a * (np.arange(cols) + 0.5)
        lats = transform.f + transform.e * (np.arange(rows) + 0.5)

        if geoid is None:
            geoid = GeoidModel.build(bounds)

        if voids.any():
            logger.info(
                "DEM: %.1f%% of the window is void (open water, most likely); "
                "filled with sea level",
                100.0 * voids.mean(),
            )
        return cls(elevation, lons, lats, geoid, voids=voids)


def _bilinear(
    grid: np.ndarray,
    grid_lons: np.ndarray,
    grid_lats: np.ndarray,
    lon: np.ndarray,
    lat: np.ndarray,
) -> np.ndarray:
    """Bilinearly sample a lon/lat grid. Latitudes may ascend or descend.

    Coordinates outside the grid are clamped to its edge rather than extrapolated:
    a ray that lands just off the loaded window should get the nearest real height,
    not a fabricated one.
    """
    lon = np.atleast_1d(np.asarray(lon, dtype=np.float64))
    lat = np.atleast_1d(np.asarray(lat, dtype=np.float64))

    def axis_index(coords: np.ndarray, values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        n = coords.size
        step = coords[1] - coords[0]
        pos = np.clip((values - coords[0]) / step, 0.0, n - 1.0)
        i0 = np.clip(np.floor(pos).astype(int), 0, n - 2)
        return i0, pos - i0

    ix, fx = axis_index(grid_lons, lon)
    iy, fy = axis_index(grid_lats, lat)

    top = grid[iy, ix] * (1 - fx) + grid[iy, ix + 1] * fx
    bottom = grid[iy + 1, ix] * (1 - fx) + grid[iy + 1, ix + 1] * fx
    return top * (1 - fy) + bottom * fy
