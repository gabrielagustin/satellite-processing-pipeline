"""The target grid: the one map grid every band is resampled onto.

Choosing it is where band co-registration actually *happens*. Each band is located
independently by the sensor model — its own line times, its own view angles, its own
terrain intersection — and then all eight are resampled onto **this single grid**.
They come out aligned not because anything shifted them into alignment, but because
each was put where it belongs.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass

import numpy as np
from rasterio.crs import CRS
from rasterio.transform import Affine
from rasterio.warp import transform as warp_transform

logger = logging.getLogger(__name__)


def utm_epsg(lon: float, lat: float) -> int:
    """EPSG code of the Universal Transverse Mercator zone containing a point.

    Never hard-code a zone: the framework is mission-agnostic, and the zone is a
    property of where the platform happened to be looking.
    """
    if abs(lat) > 80.0:
        # UTM is not defined at the poles; fall back to polar stereographic.
        return 3413 if lat > 0 else 3031
    zone = int(math.floor((lon + 180.0) / 6.0) % 60) + 1
    return (32600 if lat >= 0 else 32700) + zone


@dataclass(frozen=True)
class TargetGrid:
    """A projected raster grid: coordinate system, resolution, extent.

    Attributes
    ----------
    crs:
        Coordinate reference system of the output.
    transform:
        Affine transform mapping pixel indices to projected coordinates.
    width, height:
        Raster dimensions, pixels.
    gsd_m:
        Ground sample distance, metres.
    native_gsd_m:
        The instrument's *actual* ground sample distance for this acquisition,
        derived from the ephemeris. Reported so that a consumer can see how far
        the output grid departs from the native sampling — and so that nobody has
        to trust a data-sheet figure that may describe a different orbit entirely.
    zone_straddle:
        Whether the footprint crosses a UTM zone boundary. The centroid's zone is
        used regardless; the flag exists so the distortion is visible rather than
        silent.
    """

    crs: CRS
    transform: Affine
    width: int
    height: int
    gsd_m: float
    native_gsd_m: float
    zone_straddle: bool = False

    @property
    def bounds(self) -> tuple[float, float, float, float]:
        """``(left, bottom, right, top)`` in projected coordinates."""
        left, top = self.transform * (0, 0)
        right, bottom = self.transform * (self.width, self.height)
        return left, bottom, right, top

    @property
    def shape(self) -> tuple[int, int]:
        return self.height, self.width

    @classmethod
    def covering(
        cls,
        footprints: dict[str, tuple[float, float, float, float]],
        *,
        native_gsd_m: float,
        gsd_m: float | None = None,
        crs: CRS | str | None = None,
    ) -> TargetGrid:
        """Build the grid covering the **union** of the per-band footprints.

        Parameters
        ----------
        footprints:
            Band name to ``(min_lon, min_lat, max_lon, max_lat)``.
        native_gsd_m:
            The acquisition's true ground sample distance, from the ephemeris.
        gsd_m:
            Output resolution. Defaults to ``native_gsd_m`` rounded **up** to the
            next half metre — rounding up slightly oversamples rather than
            decimating, which avoids aliasing the source.
        crs:
            Output coordinate system. Defaults to the UTM zone of the footprint
            centroid.

        Notes
        -----
        The **union** is used, not the intersection. The bands' footprints very
        nearly coincide — the product already staggers them — so the union is barely
        larger than any one band, and it preserves every sample. A consumer wanting
        only fully co-registered coverage can crop to the intersection, which the
        quality report records.
        """
        if not footprints:
            raise ValueError("Cannot build a target grid with no footprints")

        lons = [v for f in footprints.values() for v in (f[0], f[2])]
        lats = [v for f in footprints.values() for v in (f[1], f[3])]
        min_lon, max_lon = min(lons), max(lons)
        min_lat, max_lat = min(lats), max(lats)
        centre_lon = 0.5 * (min_lon + max_lon)
        centre_lat = 0.5 * (min_lat + max_lat)

        straddle = utm_epsg(min_lon, centre_lat) != utm_epsg(max_lon, centre_lat)
        if crs is None:
            crs = CRS.from_epsg(utm_epsg(centre_lon, centre_lat))
            if straddle:
                logger.warning(
                    "Footprint straddles a UTM zone boundary; using the centroid's "
                    "zone (%s). Pass an explicit CRS to override.",
                    crs.to_string(),
                )
        crs = CRS.from_user_input(crs)

        if gsd_m is None:
            gsd_m = math.ceil(native_gsd_m * 2.0) / 2.0

        # Project the footprint corners, then take their bounding box. The corners
        # alone are not enough for a long strip -- the projected edges bow -- so the
        # boundary is densified before projecting.
        lon_edge, lat_edge = _densified_boundary(min_lon, min_lat, max_lon, max_lat)
        xs, ys = warp_transform(CRS.from_epsg(4326), crs, lon_edge, lat_edge)

        left = math.floor(min(xs) / gsd_m) * gsd_m
        right = math.ceil(max(xs) / gsd_m) * gsd_m
        bottom = math.floor(min(ys) / gsd_m) * gsd_m
        top = math.ceil(max(ys) / gsd_m) * gsd_m

        width = int(round((right - left) / gsd_m))
        height = int(round((top - bottom) / gsd_m))
        transform = Affine.translation(left, top) * Affine.scale(gsd_m, -gsd_m)

        logger.info(
            "Target grid: %s, %.2f m (native %.2f m), %d x %d px",
            crs.to_string(),
            gsd_m,
            native_gsd_m,
            width,
            height,
        )
        return cls(
            crs=crs,
            transform=transform,
            width=width,
            height=height,
            gsd_m=float(gsd_m),
            native_gsd_m=float(native_gsd_m),
            zone_straddle=straddle,
        )


def _densified_boundary(
    min_lon: float, min_lat: float, max_lon: float, max_lat: float, n: int = 32
) -> tuple[list[float], list[float]]:
    """Points along the edges of a lon/lat box, not merely its corners.

    A projected strip's edges are curved. Taking the bounding box of four projected
    corners can clip the bulge in between — a small error, but one that crops real
    data off the edge of the product, which is not a defensible way to lose it.
    """
    lon = np.linspace(min_lon, max_lon, n)
    lat = np.linspace(min_lat, max_lat, n)
    edge_lon = np.concatenate([lon, np.full(n, max_lon), lon[::-1], np.full(n, min_lon)])
    edge_lat = np.concatenate([np.full(n, min_lat), lat, np.full(n, max_lat), lat[::-1]])
    return edge_lon.tolist(), edge_lat.tolist()


def native_gsd(
    altitude_m: float, pixel_size_mm: float, focal_length_mm: float
) -> float:
    """Ground sample distance at nadir, from the geometry that actually applies.

    Derived from the ephemeris, never from a data sheet: the reference acquisition
    was taken from ~398 km against a nominal ~501 km, so the nominal figure is wrong
    by 25%. A constant here would be a bug that no test on this scene could catch,
    because the scene would simply be resampled onto a grid of the wrong size.
    """
    return altitude_m * pixel_size_mm / focal_length_mm
