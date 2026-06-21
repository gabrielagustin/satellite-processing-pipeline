"""Emit a STAC Item describing a generated product.

The STAC Item is the product's catalogue/metadata record: it lists the band
rasters as assets (with their spectral and raster properties), carries the
processing lineage, and links back to the source acquisition. An L1B product is
still in **sensor coordinates**, so the Item has a ``null`` geometry and no
``bbox`` — there is no ground footprint yet; geolocation appears at L1C.

The Item is built as a plain dict (``build_stac_item``) so it can be inspected or
tested without touching the filesystem; ``write_stac_item`` serialises it.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

from spp.core.product import Product

STAC_VERSION = "1.0.0"
EO_EXTENSION = "https://stac-extensions.github.io/eo/v1.1.0/schema.json"
RASTER_EXTENSION = "https://stac-extensions.github.io/raster/v1.1.0/schema.json"
PROCESSING_EXTENSION = "https://stac-extensions.github.io/processing/v1.1.0/schema.json"
GEOTIFF_MEDIA_TYPE = "image/tiff; application=geotiff"


def write_stac_item(
    product: Product,
    out_path: Path | str,
    *,
    quicklook_href: str | None = None,
    qa_href: str | None = None,
    created: str | None = None,
) -> Path:
    """Build and write the STAC Item for ``product`` to ``out_path``.

    Asset and link ``href``s are written **relative** to the Item, so the output
    directory stays self-contained and movable.
    """
    out_path = Path(out_path)
    item = build_stac_item(
        product,
        item_href=out_path.name,
        quicklook_href=quicklook_href,
        qa_href=qa_href,
        created=created,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(item, indent=2))
    return out_path


def build_stac_item(
    product: Product,
    *,
    item_href: str | None = None,
    quicklook_href: str | None = None,
    qa_href: str | None = None,
    created: str | None = None,
) -> dict:
    """Return the STAC Item (as a dict) describing ``product``."""
    scene_id = product.scene_id
    item_id = f"{scene_id}_{product.level}"
    band_meta = product.metadata.get("bands", {})
    units = product.units
    software = product.provenance.get("software")

    assets: dict[str, dict] = {}
    for name, path in product.bands.items():
        eo_band = {"name": name}
        cwl_nm = band_meta.get(name, {}).get("cwl_nm")
        if cwl_nm is not None and cwl_nm == cwl_nm:  # not None and not NaN
            # The EO extension expresses centre wavelength in micrometres.
            eo_band["center_wavelength"] = round(cwl_nm / 1000.0, 4)
        assets[name] = {
            "href": Path(path).name,
            "type": GEOTIFF_MEDIA_TYPE,
            "title": f"{name} {product.level} radiance",
            "roles": ["data"],
            "eo:bands": [eo_band],
            "raster:bands": [
                {"data_type": "float32", "unit": units, "nodata": "nan"}
            ],
        }

    if quicklook_href is not None:
        assets["thumbnail"] = {
            "href": quicklook_href,
            "type": "image/png",
            "title": "RGB quicklook",
            "roles": ["thumbnail", "overview"],
        }
    if qa_href is not None:
        assets["qa"] = {
            "href": qa_href,
            "type": "application/json",
            "title": "Quality assessment report",
            "roles": ["metadata"],
        }

    properties = {
        "datetime": product.metadata.get("acquired_at"),
        "created": created or _utc_now(),
        "processing:level": product.level,
        "product_type": "radiance",
    }
    if software is not None:
        properties["processing:software"] = {software["name"]: software["version"]}

    links = [
        {
            "rel": "self",
            "href": item_href or f"{item_id}.json",
            "type": "application/json",
        },
        {
            "rel": "derived_from",
            "href": f"{scene_id}.json",
            "type": "application/json",
            "title": "Source acquisition",
        },
    ]

    return {
        "type": "Feature",
        "stac_version": STAC_VERSION,
        "stac_extensions": [EO_EXTENSION, RASTER_EXTENSION, PROCESSING_EXTENSION],
        "id": item_id,
        "geometry": None,
        "properties": properties,
        "assets": assets,
        "links": links,
    }


def _utc_now() -> str:
    """Current UTC time as an RFC 3339 / STAC-style timestamp."""
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
