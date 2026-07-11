"""``spp-l1c`` — georeference, orthorectify and co-register an L1B product.

Takes the acquisition package (for the platform telemetry and the calibration) and the
L1B radiance rasters (for the pixels), and writes projected products on a common grid.

The elevation model is fetched automatically over the network unless one is supplied.
Without it the stage still runs — georeferencing on the ellipsoid — and says clearly
that the result is not orthorectified.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

from spp.geometry import dem_source
from spp.geometry.telemetry import (
    attitude_from_ancillary,
    camera_from_package,
    ephemeris_from_ancillary,
    gnss_lock,
    line_timing_from_session,
)
from spp.geometry.terrain import DEMTerrain
from spp.pipeline.l1c_pipeline import L1CPipeline

logger = logging.getLogger("spp")


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="spp-l1c",
        description=(
            "Georeference, orthorectify and co-register an L1B product onto a common "
            "projected grid."
        ),
    )
    parser.add_argument("--input", required=True, type=Path, help="Acquisition package directory")
    parser.add_argument("--l1b", required=True, type=Path, help="Directory of L1B band rasters")
    parser.add_argument("--output", required=True, type=Path, help="Directory to write L1C into")
    parser.add_argument("--bands", nargs="+", help="Subset of bands (default: all)")
    parser.add_argument(
        "--dem",
        type=Path,
        help="Elevation raster or virtual mosaic. Fetched from Copernicus GLO-30 if omitted.",
    )
    parser.add_argument(
        "--no-dem",
        action="store_true",
        help="Skip terrain entirely: georeference on the ellipsoid, do NOT orthorectify.",
    )
    parser.add_argument("--gsd", type=float, help="Output resolution, metres (default: native)")
    parser.add_argument("--reference-band", default="PAN", help="Band the others co-register onto")
    parser.add_argument(
        "--no-refine",
        action="store_true",
        help="Skip the self-calibration of the missing line-of-sight term. The bands "
        "will not be co-registered.",
    )
    parser.add_argument("--step", type=int, default=8, help="Geolocation-lattice spacing, px")
    parser.add_argument("--quiet", action="store_true")
    return parser.parse_args(argv)


def _read_package(package: Path) -> dict:
    """Load the documents the geometric level needs, by discovery, not by fixed names."""
    def load(pattern: str) -> dict:
        matches = sorted(package.glob(pattern))
        if not matches:
            raise FileNotFoundError(f"No {pattern} in {package}")
        return json.loads(matches[0].read_text())

    ancillary = load("ancillary.json")
    metadata = load("metadata.json")
    cpf = load("calibration/CPF_*.json")
    filters = load("calibration/filters_*.json")

    stac = None
    for candidate in package.glob("*.json"):
        document = json.loads(candidate.read_text())
        if isinstance(document, dict) and "stac_version" in document:
            stac = document
            break
    if stac is None:
        raise FileNotFoundError(f"No STAC item in {package}")

    sessions = metadata["Sessions"]
    session = sessions[next(iter(sessions))]
    return {
        "ancillary": ancillary,
        "session": session,
        "cpf": cpf,
        "filters": filters,
        "stac": stac,
    }


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    logging.basicConfig(
        level=logging.WARNING if args.quiet else logging.INFO,
        format="%(levelname)s  %(message)s",
    )
    logging.getLogger("rasterio").setLevel(logging.ERROR)

    package = _read_package(args.input)
    ancillary, session = package["ancillary"], package["session"]
    cpf, filters, stac = package["cpf"], package["filters"], package["stac"]

    time_sync_offset = 0.0
    for key, value in stac.get("properties", {}).items():
        if key.endswith("time_sync_offset"):
            time_sync_offset = float(value)

    rows = {entry["band"]: int(entry["row"]) for entry in cpf["radiometric"]}
    tdi = {entry["band"]: int(entry["tdi"]) for entry in cpf["radiometric"]}

    camera = camera_from_package(ancillary, cpf.get("geometric"), rows, tdi)
    ephemeris = ephemeris_from_ancillary(ancillary)
    attitude = attitude_from_ancillary(ancillary)
    timing = line_timing_from_session(
        session,
        {entry["name"]: int(entry["id"]) for entry in filters},
        time_sync_offset_s=time_sync_offset,
    )

    if not gnss_lock(ancillary):
        logger.warning(
            "This acquisition has no satellite-navigation lock: its positions are "
            "propagated, not measured. The absolute geolocation is UNVALIDATED against "
            "anything independent of the telemetry."
        )

    terrain = None
    if not args.no_dem:
        dem_path = args.dem
        if dem_path is None:
            logger.info("Fetching an elevation model over the footprint...")
            dem_path = dem_source.build_mosaic(
                tuple(stac["bbox"]), args.output / "cache" / "dem.vrt"
            )
        if dem_path is not None:
            terrain = DEMTerrain.from_raster(str(dem_path), tuple(stac["bbox"]))
    else:
        logger.warning(
            "Running without an elevation model: the product will be georeferenced but "
            "NOT orthorectified. Relief will be displaced."
        )

    pipeline = L1CPipeline(
        camera,
        ephemeris,
        attitude,
        timing,
        terrain=terrain,
        reference_band=args.reference_band,
        refine=not args.no_refine,
        gsd_m=args.gsd,
        step=args.step,
    )
    result = pipeline.run(args.l1b, args.output, bands=args.bands)

    _summarise(result)
    return 0


def _summarise(result) -> None:
    """Print what was produced — and, as prominently, what was not established."""
    qa = result.qa
    grid = result.grid

    print()
    print(f"  L1C: {len(result.products)} band(s) on {grid.crs.to_string()} "
          f"@ {grid.gsd_m} m  ({grid.width} x {grid.height} px)")
    print(f"       native sampling {qa['grid']['native_gsd_m']} m, "
          f"lattice error {qa['interpolation']['max_error_px']} px")
    print(f"       terrain: {qa['terrain']['source']}")

    coregistration = qa.get("coregistration", {})
    if coregistration:
        fitted = [b for b, v in coregistration.items() if v["fitted"]]
        refused = {b: v for b, v in coregistration.items() if not v["fitted"]}
        print()
        print(f"  co-registration: {len(fitted)}/{len(coregistration)} band(s) calibrated")
        for band, value in coregistration.items():
            if value["fitted"]:
                print(
                    f"       {band:5s} systematic {value['systematic_before_px']:6.2f} -> "
                    f"{value['systematic_after_px']:.3f} px   scatter {value['scatter_px']:.2f} px"
                )
        for band, value in refused.items():
            print(f"       {band:5s} REFUSED — {value['reason']}")

    if qa["flags"]:
        print()
        print("  NOT ESTABLISHED BY THIS RUN:")
        explanations = {
            "no_gnss_lock": "the platform had no navigation lock; positions are propagated",
            "absolute_accuracy_unvalidated": (
                "absolute geolocation has not been checked against anything independent "
                "of the telemetry"
            ),
            "conventions_assumed": (
                "two telemetry conventions (scan direction, column sign) are mirrors that "
                "geometry cannot see; they are assumed"
            ),
            "coregistration_not_achieved": (
                "one or more bands could not be co-registered — see the refusals above"
            ),
            "zone_straddle": "the footprint crosses a projection-zone boundary",
        }
        for flag in qa["flags"]:
            print(f"       - {flag}: {explanations.get(flag, '')}")
    print()


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
