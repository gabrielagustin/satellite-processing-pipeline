"""Resampling a band onto the target grid, through its geolocation array.

The geolocation array says where each lattice node landed on the ground. GDAL
inverts that mapping and pulls each output pixel from the source. This module is
the thin, careful layer around it.

Careful, because resampling a *physical* quantity is not the same as resampling a
picture.
"""

from __future__ import annotations

import logging
import os
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
    "compress": "zstd",
    "zstd_level": 1,
    "predictor": 3,  # floating-point predictor
    "BIGTIFF": "IF_SAFER",
}
"""Cloud-Optimized GeoTIFF layout.

**Zstandard, not Deflate.** Measured on a real warped band (9,773 x 29,403, 39% data):
writing plus overviews takes **5.7 s with Zstandard level 1 against 20.4 s with Deflate**,
and the files are the same size to within half a percent (432 MB against 434 MB). There is
no trade here to think about — Deflate was simply costing 15 seconds a band for nothing.

Level 1, not 9: level 9 costs another 2.6 s and saves 2% of the file. The floating-point
predictor is what is actually doing the compressing, and it does it before the codec sees
the data.
"""

WARP_THREADS = max(1, (os.cpu_count() or 2) - 1)
"""Threads for the resampling itself.

GDAL's warper is single-threaded unless told otherwise, and it was. On the reference
acquisition the resampling of one band took **24.6 s on one core and 6.5 s on ten** — a
free 3.8x that had been left on the table because the default is 1.
"""


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

    def resample(
        self,
        source_path: Path,
        geoloc: GeolocationGrid,
        target: TargetGrid,
    ) -> np.ndarray:
        """Resample onto the target grid and **return the array**, writing nothing.

        Exposed separately from :meth:`warp` because the product is a multi-band stack:
        writing each band to its own file and then reading all of them back to build the
        stack costs a compress, a decompress and two passes over 1.7 GB, for data that was
        already in memory. On the reference acquisition that round trip was 42 s of a
        187 s run.
        """
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
            num_threads=WARP_THREADS,
        )

        valid = float(np.isfinite(destination).mean())
        logger.info(
            "      resampled %s (%.1f%% of the grid carries data)",
            source_path.name,
            100.0 * valid,
        )
        return destination

    def warp(
        self,
        source_path: Path,
        geoloc: GeolocationGrid,
        target: TargetGrid,
        output_path: Path,
        *,
        overviews: bool = True,
    ) -> Path:
        """Resample and write a single-band raster.

        The array lives and dies inside this call. Returning it instead — to assemble the
        stack from memory and skip a disk round trip — keeps the previous band's 1.15 GB
        destination alive while the next one is allocated, and the process gets killed. The
        round trip is cheaper than the memory.
        """
        destination = self.resample(source_path, geoloc, target)

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
            if overviews:
                dst.build_overviews([2, 4, 8, 16, 32], Resampling.average)
        return output_path
