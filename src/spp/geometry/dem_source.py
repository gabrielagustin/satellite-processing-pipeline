"""Locating an elevation model for a footprint, without making the user find one.

The geometric level needs terrain. Asking a user to go and assemble the right tiles is
a good way to have the level run without terrain, so this module builds a virtual mosaic
over the footprint from a public, unauthenticated source — and falls back gracefully,
and loudly, when it cannot.
"""

from __future__ import annotations

import logging
import math
from pathlib import Path

logger = logging.getLogger(__name__)

COPERNICUS_GLO30 = (
    "/vsicurl/https://copernicus-dem-30m.s3.amazonaws.com/"
    "Copernicus_DSM_COG_10_{ns}{lat:02d}_00_{ew}{lon:03d}_00_DEM/"
    "Copernicus_DSM_COG_10_{ns}{lat:02d}_00_{ew}{lon:03d}_00_DEM.tif"
)
"""Copernicus GLO-30: global, 30 m, EGM2008-referenced, no authentication."""


def tile_urls(bounds: tuple[float, float, float, float]) -> list[str]:
    """Every one-degree elevation tile the footprint touches."""
    min_lon, min_lat, max_lon, max_lat = bounds
    urls = []
    for lat in range(math.floor(min_lat), math.floor(max_lat) + 1):
        for lon in range(math.floor(min_lon), math.floor(max_lon) + 1):
            urls.append(
                COPERNICUS_GLO30.format(
                    ns="N" if lat >= 0 else "S",
                    lat=abs(lat),
                    ew="E" if lon >= 0 else "W",
                    lon=abs(lon),
                )
            )
    return urls


def build_mosaic(bounds: tuple[float, float, float, float], destination: Path) -> Path | None:
    """Assemble a virtual mosaic of the elevation tiles covering ``bounds``.

    Tiles that do not exist are **skipped, not fatal**: the source omits tiles that are
    entirely open ocean, and a scene over water legitimately has none. A footprint with
    no tiles at all yields ``None``, and the caller falls back to the ellipsoid — which
    is exactly right over water and wrong everywhere else, so it is said out loud.

    Returns
    -------
    pathlib.Path or None
        Path to a GDAL virtual raster, or ``None`` if no tile could be reached.
    """
    import rasterio
    from rasterio.errors import RasterioIOError

    reachable = []
    for url in tile_urls(bounds):
        try:
            with rasterio.open(url):
                reachable.append(url)
        except RasterioIOError:
            logger.debug("Elevation tile absent (open ocean, most likely): %s", url)

    if not reachable:
        logger.warning(
            "No elevation tiles could be reached for this footprint. Falling back to "
            "the ellipsoid: correct over water, and displacing every metre of relief "
            "elsewhere. Pass --dem to supply one, or check network access."
        )
        return None

    destination.parent.mkdir(parents=True, exist_ok=True)
    _write_vrt(reachable, destination)
    logger.info("Elevation: %d tile(s) mosaicked into %s", len(reachable), destination.name)
    return destination


def _write_vrt(sources: list[str], destination: Path) -> None:
    """Build the virtual mosaic with GDAL, through rasterio's own GDAL."""
    from rasterio.shutil import copy as rio_copy
    from rasterio.vrt import WarpedVRT
    import rasterio

    # `gdalbuildvrt` is not importable, and shelling out would add a binary dependency,
    # so the mosaic is assembled as a plain VRT document. Every Copernicus tile shares
    # the same CRS and resolution, which is what makes this safe.
    with rasterio.open(sources[0]) as first:
        res = first.res
        crs = first.crs
        dtype = first.dtypes[0]
        nodata = first.nodata

    bounds = []
    for source in sources:
        with rasterio.open(source) as src:
            bounds.append(src.bounds)

    left = min(b.left for b in bounds)
    right = max(b.right for b in bounds)
    bottom = min(b.bottom for b in bounds)
    top = max(b.top for b in bounds)

    width = int(round((right - left) / res[0]))
    height = int(round((top - bottom) / res[1]))

    lines = [
        f'<VRTDataset rasterXSize="{width}" rasterYSize="{height}">',
        f"  <SRS>{crs.to_wkt()}</SRS>",
        f"  <GeoTransform>{left}, {res[0]}, 0.0, {top}, 0.0, {-res[1]}</GeoTransform>",
        f'  <VRTRasterBand dataType="{_gdal_type(dtype)}" band="1">',
    ]
    if nodata is not None:
        lines.append(f"    <NoDataValue>{nodata}</NoDataValue>")

    for source, b in zip(sources, bounds):
        x_off = int(round((b.left - left) / res[0]))
        y_off = int(round((top - b.top) / res[1]))
        x_size = int(round((b.right - b.left) / res[0]))
        y_size = int(round((b.top - b.bottom) / res[1]))
        lines += [
            "    <SimpleSource>",
            f'      <SourceFilename relativeToVRT="0">{source}</SourceFilename>',
            "      <SourceBand>1</SourceBand>",
            f'      <SrcRect xOff="0" yOff="0" xSize="{x_size}" ySize="{y_size}" />',
            f'      <DstRect xOff="{x_off}" yOff="{y_off}" xSize="{x_size}" ySize="{y_size}" />',
            "    </SimpleSource>",
        ]

    lines += ["  </VRTRasterBand>", "</VRTDataset>"]
    destination.write_text("\n".join(lines))


def _gdal_type(numpy_dtype: str) -> str:
    """NumPy dtype name to the GDAL type name a VRT expects."""
    return {
        "float32": "Float32",
        "float64": "Float64",
        "int16": "Int16",
        "uint16": "UInt16",
        "int32": "Int32",
    }.get(numpy_dtype, "Float32")
