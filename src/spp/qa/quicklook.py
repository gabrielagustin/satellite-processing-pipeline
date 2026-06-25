"""Generate an RGB quicklook PNG from calibrated band rasters.

A quicklook is a small, contrast-stretched RGB preview for visual inspection —
it has no quantitative use. It is read at a decimated resolution (via the
rasters' internal overviews, so it stays fast and low-memory) and each channel
is independently percentile-stretched to bring out detail.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.errors import NotGeoreferencedWarning


def write_quicklook(
    band_paths: dict[str, Path],
    out_path: Path | str,
    *,
    rgb: tuple[str, str, str] = ("R", "G", "B"),
    max_size: int = 2048,
    stretch: tuple[float, float] = (2.0, 98.0),
) -> Path:
    """Write an RGB quicklook PNG built from three bands.

    Parameters
    ----------
    band_paths:
        Mapping of band name to raster path (must include the ``rgb`` bands).
    out_path:
        Destination PNG path.
    rgb:
        Band names to map to the red, green and blue channels.
    max_size:
        Maximum size (pixels) of the longest axis of the preview.
    stretch:
        Lower/upper percentiles for the per-channel contrast stretch.

    Returns
    -------
    pathlib.Path
        The written PNG path.
    """
    missing = [b for b in rgb if b not in band_paths]
    if missing:
        raise KeyError(f"Missing bands for quicklook: {missing}")

    channels = [
        _stretch(_read_decimated(band_paths[name], max_size), *stretch)
        for name in rgb
    ]
    image = np.dstack(channels)
    image = np.nan_to_num(image, nan=0.0)  # NoData -> black
    rgb_uint8 = (np.clip(image, 0.0, 1.0) * 255).astype(np.uint8)

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    _save_png(rgb_uint8, out_path)
    return out_path


def _read_decimated(path: Path, max_size: int) -> np.ndarray:
    """Read a single band at a decimated resolution, masking NoData to NaN."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", NotGeoreferencedWarning)
        with rasterio.open(path) as ds:
            scale = max(max(ds.height, ds.width) / max_size, 1.0)
            out_h = max(1, round(ds.height / scale))
            out_w = max(1, round(ds.width / scale))
            data = ds.read(
                1,
                out_shape=(out_h, out_w),
                resampling=Resampling.average,
            ).astype(np.float32)
            nodata = ds.nodata
    if nodata is not None and not np.isnan(nodata):
        data = np.where(data == nodata, np.nan, data)
    return data


def _stretch(channel: np.ndarray, p_low: float, p_high: float) -> np.ndarray:
    """Percentile contrast stretch of a single channel to [0, 1]."""
    finite = np.isfinite(channel)
    if not finite.any():
        return np.zeros_like(channel)
    lo, hi = np.percentile(channel[finite], [p_low, p_high])
    if hi <= lo:
        hi = lo + 1.0
    return np.clip((channel - lo) / (hi - lo), 0.0, 1.0)


def _save_png(rgb_uint8: np.ndarray, out_path: Path) -> None:
    """Save an HxWx3 uint8 array as a PNG."""
    # Imported lazily so the QA layer stays import-light; matplotlib is a core
    # dependency, so the import is guaranteed to succeed.
    from matplotlib import image as mpimage

    mpimage.imsave(out_path, rgb_uint8)
