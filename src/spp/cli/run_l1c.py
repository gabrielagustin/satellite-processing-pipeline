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
import os
import sys
import warnings
from datetime import date
from pathlib import Path

from rasterio.errors import NotGeoreferencedWarning

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
    parser.add_argument(
        "--no-parity-check",
        action="store_true",
        help="Do not resolve the scan/column mirrors against the terrain; use the "
        "defaults. Only sensible if the scene has no land/water contrast.",
    )
    parser.add_argument(
        "--no-reference-image",
        action="store_true",
        help="Do not search for a reference orthoimage. Absolute geolocation then falls "
        "back to the elevation model's coastline, which is coarser (tens of metres).",
    )
    parser.add_argument(
        "--no-absolute-correction",
        action="store_true",
        help="Do not correct the absolute pointing bias against the terrain. The product "
        "keeps whatever geolocation bias the telemetry carries (~1 km without a "
        "navigation lock).",
    )
    parser.add_argument(
        "--no-scene-correction",
        action="store_true",
        help="Skip the scene-local (attitude) correction. The instrument calibration is "
        "still estimated, but the bands will not fully co-register.",
    )
    parser.add_argument(
        "--per-band",
        action="store_true",
        help="Also write one raster per band. Redundant with stack.tif, and it costs a "
        "compress and a set of overviews each.",
    )
    parser.add_argument(
        "--no-stack",
        action="store_true",
        help="Do not write the multi-band stack.tif (per-band rasters only).",
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
    # L1B rasters are in sensor coordinates and carry no CRS *by design*. GDAL says so
    # on every open; it is not news, and it drowns the log.
    warnings.filterwarnings("ignore", category=NotGeoreferencedWarning)

    # Remote Cloud-Optimized GeoTIFFs are read once (the elevation model, the reference
    # orthoimage) and then cached locally. These make that one read as fast as it can be.
    os.environ.setdefault("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
    os.environ.setdefault("GDAL_HTTP_MULTIPLEX", "YES")
    os.environ.setdefault("GDAL_HTTP_VERSION", "2")
    os.environ.setdefault("VSI_CACHE", "TRUE")
    os.environ.setdefault("GDAL_NUM_THREADS", "ALL_CPUS")

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

    # The reference is searched inside the pipeline, once it knows which band it will
    # match -- so that our band and the reference's are the same WAVELENGTH.
    reference_search = None
    if not args.no_absolute_correction and not args.no_reference_image and terrain is not None:
        reference_search = (
            tuple(stac["bbox"]),
            date.fromisoformat(stac["properties"]["datetime"][:10]),
            args.output / "cache" / "reference.vrt",
        )

    # Central wavelength per band, from the imager configuration. It is what lets the
    # reference image be requested in the SAME band as ours: matching a 665 nm red
    # against a reference's 842 nm near-infrared correlates two different pictures of the
    # same ground, and produces a confident, precise, wrong answer.
    band_ids = {entry["name"]: int(entry["id"]) for entry in filters}
    cwls = session["ImagerConfiguration"].get("BandCWL", [])
    wavelengths = {
        name: float(cwls[band_id])
        for name, band_id in band_ids.items()
        if band_id < len(cwls)
    }

    pipeline = L1CPipeline(
        camera,
        ephemeris,
        attitude,
        timing,
        terrain=terrain,
        reference_band=args.reference_band,
        refine=not args.no_refine,
        resolve_parities=not args.no_parity_check,
        correct_absolute=not args.no_absolute_correction,
        band_wavelengths=wavelengths,
        reference_search=reference_search,
        scene_correct=not args.no_scene_correction,
        stack=not args.no_stack,
        per_band=args.per_band,
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
    if result.stack:
        print(f"  L1C: stack.tif — {len(result.products)} band(s), co-registered, "
              f"{grid.crs.to_string()} @ {grid.gsd_m} m  ({grid.width} x {grid.height} px)")
    else:
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

    absolute_qa = qa.get("absolute")
    if absolute_qa and absolute_qa.get("fitted"):
        print()
        print("  absolute geolocation — the only reference in this run the telemetry")
        print("  did not produce:")
        print(f"       reference       {absolute_qa.get('reference', '?')}")
        print(f"       offset before   {absolute_qa['offset_before_m']:.0f} m")
        print(f"       boresight       roll {absolute_qa['roll_urad']:+.0f} urad, "
              f"pitch {absolute_qa['pitch_urad']:+.0f} urad")
        print(f"       NOTE: the along-track part is equally explained by a clock offset "
              f"of {absolute_qa['equivalent_clock_offset_s']:+.3f} s.")
        print( "             They are not separable from one strip; pitch is a convention.")

    drift_qa = qa.get("drift")
    if drift_qa and drift_qa.get("fitted"):
        print()
        print("  along-track drift (a rigid boresight cannot remove this — it grows):")
        print(f"       residual  {drift_qa['residual_before_m']:.0f} m -> "
              f"{drift_qa['residual_after_m']:.0f} m   over {drift_qa['n_windows']} windows")

    scene_correction = qa.get("scene_correction", {})
    fitted_scene = [b for b, v in scene_correction.items() if isinstance(v, dict) and v.get("fitted")]
    if fitted_scene:
        print()
        print("  scene correction (THIS acquisition's pointing, not the instrument —")
        print("  applied to the product, never written to the calibration):")
        for band in fitted_scene:
            value = scene_correction[band]
            print(f"       {band:5s} band-to-band RMS {value['rms_before_px']:6.2f} -> "
                  f"{value['rms_after_px']:.2f} px")

    timings = qa.get("timings_s", {})
    if timings:
        total = sum(timings.values())
        print()
        print(f"  where the time went ({total:.0f} s total):")
        for title, seconds in list(timings.items())[:5]:
            print(f"       {seconds:6.1f} s  ({100 * seconds / total:4.1f}%)  {title}")

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
                "geometry cannot see, and the terrain could not settle them either"
            ),
            "absolute_accuracy_terrain_only": (
                "the pointing bias was corrected against a coastline in the elevation "
                "model -- a real independent reference, but good only to the DEM's "
                "resolution and the shoreline's sharpness (tens of metres). It has NOT "
                "been checked against a reference orthoimage"
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
