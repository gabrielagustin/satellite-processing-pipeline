"""GeoTIFF writer for calibrated raster bands.

Writes each band as a tiled, compressed GeoTIFF, one block of lines at a time,
so that large rasters never need to be held in memory in full. Outputs carry no
CRS or geotransform: an L1B product is still in sensor (detector) coordinates.
Internal overviews are built on close to mirror the input COGs and keep
quicklooks fast.
"""

from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
from typing import Iterator

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.windows import Window


class GeoTIFFWriter:
    """Write calibrated bands as tiled, compressed GeoTIFFs.

    Parameters
    ----------
    output_dir:
        Directory the band rasters are written to (created if missing).
    compress:
        GeoTIFF compression (e.g. ``"deflate"``, ``"lzw"``).
    blocksize:
        Internal tile size (square). Matching the input COG block size (512)
        keeps windowed read/write aligned.
    build_overviews:
        Whether to build internal overviews when a band is finalised.
    overview_levels:
        Decimation factors for the overviews.
    suffix:
        File extension for the written rasters.
    """

    def __init__(
        self,
        output_dir: Path | str,
        *,
        compress: str = "deflate",
        blocksize: int = 512,
        build_overviews: bool = True,
        overview_levels: tuple[int, ...] = (2, 4, 8, 16, 32),
        suffix: str = ".tif",
    ) -> None:
        self.output_dir = Path(output_dir)
        self.compress = compress
        self.blocksize = blocksize
        self.build_overviews = build_overviews
        self.overview_levels = overview_levels
        self.suffix = suffix

    @contextmanager
    def open_band(
        self,
        name: str,
        *,
        width: int,
        height: int,
        dtype: np.dtype | str,
        nodata: float | None,
    ) -> Iterator["_BandWriter"]:
        """Open ``name`` for windowed writing, yielding a band-writer handle."""
        self.output_dir.mkdir(parents=True, exist_ok=True)
        path = self.output_dir / f"{name}{self.suffix}"
        profile = {
            "driver": "GTiff",
            "width": width,
            "height": height,
            "count": 1,
            "dtype": np.dtype(dtype).name,
            "nodata": nodata,
            "tiled": True,
            "blockxsize": self.blocksize,
            "blockysize": self.blocksize,
            "compress": self.compress,
            "BIGTIFF": "IF_SAFER",
        }
        dataset = rasterio.open(path, "w", **profile)
        try:
            yield _BandWriter(dataset, path, width)
            if self.build_overviews:
                dataset.build_overviews(
                    list(self.overview_levels), Resampling.average
                )
                dataset.update_tags(ns="rio_overview", resampling="average")
        finally:
            dataset.close()


class _BandWriter:
    """Handle for writing one band raster block by block."""

    def __init__(self, dataset, path: Path, width: int) -> None:
        self._ds = dataset
        self.path = path
        self.width = width

    def write_block(self, array: np.ndarray, line_start: int) -> None:
        """Write a ``(lines, width)`` block whose first row is ``line_start``."""
        if array.ndim != 2:
            raise ValueError(f"Expected a 2-D block, got shape {array.shape}")
        if array.shape[1] != self.width:
            raise ValueError(
                f"Block width {array.shape[1]} != raster width {self.width}"
            )
        window = Window(0, line_start, self.width, array.shape[0])
        self._ds.write(array, 1, window=window)
