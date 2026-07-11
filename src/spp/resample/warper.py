"""Resampling a band onto the target grid, through its geolocation array.

The geolocation array says where each lattice node landed on the ground. GDAL
inverts that mapping and pulls each output pixel from the source. This module is
the thin, careful layer around it.

Careful, because resampling a *physical* quantity is not the same as resampling a
picture.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from pathlib import Path

import numpy as np
import rasterio
from rasterio.crs import CRS
from rasterio.enums import Resampling
from rasterio.warp import reproject

from spp.geometry.geoloc_grid import GeolocationGrid
from spp.resample.grid import TargetGrid

logger = logging.getLogger(__name__)

WGS84 = CRS.from_epsg(4326)

COG_PROFILE = {
    "driver": "GTiff",
    "tiled": True,
    "blockxsize": 512,
    "blockysize": 512,
    "compress": "deflate",
    "predictor": 3,  # floating-point predictor
    "BIGTIFF": "IF_SAFER",
}
"""Cloud-Optimized GeoTIFF layout, matching the L1B products."""


class Warper(ABC):
    """Resamples one band onto a target grid."""

    @abstractmethod
    def warp(
        self,
        source_path: Path,
        geoloc: GeolocationGrid,
        target: TargetGrid,
        output_path: Path,
    ) -> Path:
        """Resample ``source_path`` onto ``target``, writing ``output_path``."""


class GeolocWarper(Warper):
    """Warp through a geolocation array.

    Parameters
    ----------
    resampling:
        Kernel. **Bilinear** by default, and cubic on request. Lanczos is
        deliberately *not* offered: its negative lobes ring around saturated pixels
        and NoData edges, and ringing in a radiance product is not a cosmetic
        artefact — it manufactures values the sensor never measured.
    nodata:
        Output NoData. ``NaN``, matching the L1B convention.
    """

    def __init__(
        self,
        *,
        resampling: Resampling = Resampling.bilinear,
        nodata: float = float("nan"),
    ) -> None:
        if resampling not in (Resampling.bilinear, Resampling.cubic, Resampling.nearest):
            raise ValueError(
                f"{resampling.name} is not offered for radiance: only bilinear, "
                "cubic and nearest preserve the physics. Kernels with negative "
                "lobes ring around saturation and NoData."
            )
        self.resampling = resampling
        self.nodata = nodata

    def warp(
        self,
        source_path: Path,
        geoloc: GeolocationGrid,
        target: TargetGrid,
        output_path: Path,
    ) -> Path:
        with rasterio.open(source_path) as src:
            if (src.height, src.width) != (geoloc.n_lines, geoloc.n_columns):
                raise ValueError(
                    f"Geolocation grid describes a {geoloc.n_lines} x "
                    f"{geoloc.n_columns} raster but {source_path.name} is "
                    f"{src.height} x {src.width}. A grid built for a different band "
                    "would silently place this one on the wrong ground."
                )
            source = src.read(1).astype(np.float32)
            src_nodata = src.nodata

        # NoData must not be averaged into its neighbours: an interpolation kernel
        # straddling the edge of valid data would blend a sentinel into a physical
        # radiance. Marking it NaN first makes GDAL exclude it instead.
        if src_nodata is not None and np.isfinite(src_nodata):
            source[source == src_nodata] = np.nan

        destination = np.full(target.shape, self.nodata, dtype=np.float32)

        reproject(
            source,
            destination,
            src_geoloc_array=geoloc.as_array,
            src_crs=WGS84,
            src_nodata=np.nan,
            dst_crs=target.crs,
            dst_transform=target.transform,
            dst_nodata=self.nodata,
            resampling=self.resampling,
        )

        profile = {
            **COG_PROFILE,
            "height": target.height,
            "width": target.width,
            "count": 1,
            "dtype": "float32",
            "crs": target.crs,
            "transform": target.transform,
            "nodata": self.nodata,
        }
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with rasterio.open(output_path, "w", **profile) as dst:
            dst.write(destination, 1)
            dst.build_overviews([2, 4, 8, 16, 32], Resampling.average)

        valid = float(np.isfinite(destination).mean())
        logger.info(
            "Warped %s -> %s (%.1f%% of the grid carries data)",
            source_path.name,
            output_path.name,
            100.0 * valid,
        )
        return output_path
