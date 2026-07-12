"""Finding an external reference orthoimage, so absolute accuracy can be measured.

Every other check in this level is *internal*. The bands are co-registered onto each
other. The footprint is compared against a catalogue geometry the provider derived from
the same telemetry we consume. The terrain's coastline is a real external reference, but
a coarse one — it is good to the elevation model's own resolution and to how sharply a
shoreline happens to be defined.

None of that can tell you where the product actually is to within a pixel. Only a
**reference orthoimage** can: an independently georeferenced picture of the same ground,
with real texture, matched against ours.

Choosing one
------------
The reference is searched by catalogue over the footprint, filtered for cloud, and ranked
by **how close it is in time to the acquisition**. Time matters more than it looks:
coastlines move with the tide, fields are ploughed, rivers shift. A reference from the
same day is worth far more than a sharper one from a different season, and the search
says which it found rather than silently taking the first hit.

Bands are matched by **wavelength**, not by name. Our red sits at 665 nm and so does the
reference's; our near-infrared at 842 nm and so does theirs. Matching a red band against
a near-infrared one would correlate two different pictures of the same place.
"""

from __future__ import annotations

import json
import logging
import ssl
import urllib.request
from dataclasses import dataclass
from datetime import date
from pathlib import Path

logger = logging.getLogger(__name__)

STAC_API = "https://earth-search.aws.element84.com/v1/search"
"""A public, unauthenticated catalogue of Sentinel-2 surface-reflectance orthoimages."""

COLLECTION = "sentinel-2-l2a"

BAND_BY_WAVELENGTH = {
    490: "blue",
    560: "green",
    665: "red",
    842: "nir",
}
"""Reference asset for a central wavelength, in nanometres.

Matched by wavelength because that is what makes two pictures comparable. A 30 nm
mismatch is fine; matching red against near-infrared is not a mismatch, it is a different
scene.
"""


@dataclass(frozen=True)
class ReferenceImage:
    """An external orthoimage to measure absolute accuracy against.

    Attributes
    ----------
    path:
        A virtual mosaic of the reference tiles, openable by GDAL.
    acquired:
        Date of the reference.
    days_from_target:
        How far it is, in days, from the acquisition being corrected. **Report this.** A
        reference from another season measures a landscape that has changed.
    cloud_cover:
        Mean cloud cover of the mosaicked tiles, percent.
    wavelength_nm:
        Central wavelength of the band fetched.
    n_tiles:
        How many tiles were mosaicked.
    """

    path: Path
    acquired: date
    days_from_target: int
    cloud_cover: float
    wavelength_nm: int
    n_tiles: int


def find(
    bounds: tuple[float, float, float, float],
    target_date: date,
    destination: Path,
    *,
    wavelength_nm: int = 842,
    max_cloud_cover: float = 15.0,
    search_days: int = 60,
) -> ReferenceImage | None:
    """Search for the reference orthoimage closest in time to the acquisition.

    Parameters
    ----------
    bounds:
        ``(min_lon, min_lat, max_lon, max_lat)`` of the footprint.
    target_date:
        Date of the acquisition being corrected.
    destination:
        Where to write the virtual mosaic.
    wavelength_nm:
        Central wavelength to match. Near-infrared by default: water is near-black in it,
        so the land/water boundary is sharp, and vegetation is bright, so there is texture.
    max_cloud_cover:
        Reject candidates cloudier than this.
    search_days:
        How far either side of ``target_date`` to look.

    Returns
    -------
    ReferenceImage or None
        ``None`` when nothing usable was found — which is a legitimate outcome, and the
        caller must treat it as "absolute accuracy not measured", never as "no error".
    """
    asset = BAND_BY_WAVELENGTH.get(wavelength_nm)
    if asset is None:
        closest = min(BAND_BY_WAVELENGTH, key=lambda w: abs(w - wavelength_nm))
        logger.info(
            "No reference band at %d nm; using the closest, %d nm.", wavelength_nm, closest
        )
        wavelength_nm, asset = closest, BAND_BY_WAVELENGTH[closest]

    features = _search(bounds, target_date, max_cloud_cover, search_days)
    if not features:
        logger.warning(
            "No reference orthoimage found over the footprint within %d days at under "
            "%.0f%% cloud. Absolute accuracy cannot be measured against an independent "
            "image.",
            search_days,
            max_cloud_cover,
        )
        return None

    # Group by date, and prefer the date closest in time that covers the footprint.
    by_date: dict[str, list] = {}
    for feature in features:
        by_date.setdefault(feature["properties"]["datetime"][:10], []).append(feature)

    best_day = min(by_date, key=lambda d: abs((date.fromisoformat(d) - target_date).days))
    chosen = by_date[best_day]
    days = abs((date.fromisoformat(best_day) - target_date).days)

    hrefs = [f["assets"][asset]["href"] for f in chosen if asset in f["assets"]]
    if not hrefs:
        logger.warning("Reference scenes carry no %s asset; cannot measure absolute accuracy.", asset)
        return None

    cloud = sum(f["properties"].get("eo:cloud_cover", 0.0) for f in chosen) / len(chosen)
    _write_vrt([f"/vsicurl/{h}" for h in hrefs], destination)

    if days > 30:
        logger.warning(
            "The closest reference orthoimage is %d days from the acquisition. The "
            "landscape will have changed; treat the absolute figure as indicative.",
            days,
        )
    logger.info(
        "Reference orthoimage: %s (%d day(s) from the acquisition, %.1f%% cloud, "
        "%d tile(s), %d nm)",
        best_day,
        days,
        cloud,
        len(hrefs),
        wavelength_nm,
    )
    return ReferenceImage(
        path=destination,
        acquired=date.fromisoformat(best_day),
        days_from_target=days,
        cloud_cover=cloud,
        wavelength_nm=wavelength_nm,
        n_tiles=len(hrefs),
    )


def _search(
    bounds: tuple[float, float, float, float],
    target_date: date,
    max_cloud_cover: float,
    search_days: int,
) -> list[dict]:
    """Query the catalogue. Network failure is a miss, not a crash."""
    from datetime import timedelta

    start = target_date - timedelta(days=search_days)
    end = target_date + timedelta(days=search_days)
    body = {
        "collections": [COLLECTION],
        "bbox": list(bounds),
        "datetime": f"{start.isoformat()}T00:00:00Z/{end.isoformat()}T23:59:59Z",
        "query": {"eo:cloud_cover": {"lt": max_cloud_cover}},
        "limit": 100,
    }
    request = urllib.request.Request(
        STAC_API,
        data=json.dumps(body).encode(),
        headers={"Content-Type": "application/json"},
    )
    try:
        import certifi

        context = ssl.create_default_context(cafile=certifi.where())
    except ImportError:  # pragma: no cover - certifi ships with rasterio's stack
        context = ssl.create_default_context()

    try:
        with urllib.request.urlopen(request, context=context, timeout=60) as response:
            return json.load(response).get("features", [])
    except Exception as error:  # noqa: BLE001 - any network failure is simply a miss
        logger.warning("Reference catalogue unreachable (%s); skipping.", error)
        return []


def _write_vrt(sources: list[str], destination: Path) -> None:
    """Mosaic the reference tiles into a virtual raster."""
    import rasterio

    with rasterio.open(sources[0]) as first:
        res, crs, dtype, nodata = first.res, first.crs, first.dtypes[0], first.nodata

    boxes = []
    for source in sources:
        with rasterio.open(source) as src:
            if src.crs != crs:
                logger.warning(
                    "Reference tiles span more than one projection; using only those in "
                    "%s.",
                    crs.to_string(),
                )
                continue
            boxes.append((source, src.bounds))

    left = min(b.left for _, b in boxes)
    right = max(b.right for _, b in boxes)
    bottom = min(b.bottom for _, b in boxes)
    top = max(b.top for _, b in boxes)
    width = int(round((right - left) / res[0]))
    height = int(round((top - bottom) / res[1]))

    lines = [
        f'<VRTDataset rasterXSize="{width}" rasterYSize="{height}">',
        f"  <SRS>{crs.to_wkt()}</SRS>",
        f"  <GeoTransform>{left}, {res[0]}, 0.0, {top}, 0.0, {-res[1]}</GeoTransform>",
        f'  <VRTRasterBand dataType="{"UInt16" if dtype == "uint16" else "Float32"}" band="1">',
    ]
    if nodata is not None:
        lines.append(f"    <NoDataValue>{nodata}</NoDataValue>")
    for source, box in boxes:
        x_off = int(round((box.left - left) / res[0]))
        y_off = int(round((top - box.top) / res[1]))
        x_size = int(round((box.right - box.left) / res[0]))
        y_size = int(round((box.top - box.bottom) / res[1]))
        lines += [
            "    <SimpleSource>",
            f'      <SourceFilename relativeToVRT="0">{source}</SourceFilename>',
            "      <SourceBand>1</SourceBand>",
            f'      <SrcRect xOff="0" yOff="0" xSize="{x_size}" ySize="{y_size}" />',
            f'      <DstRect xOff="{x_off}" yOff="{y_off}" xSize="{x_size}" ySize="{y_size}" />',
            "    </SimpleSource>",
        ]
    lines += ["  </VRTRasterBand>", "</VRTDataset>"]

    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("\n".join(lines))
