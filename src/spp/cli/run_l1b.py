"""Command-line entry point: produce an L1B radiance product from a package.

Usage
-----
    spp-l1b --input <package_dir> --output <output_dir> [options]

Reads an acquisition package, calibrates every band to at-aperture spectral
radiance, writes one GeoTIFF per band plus a ``qa_report.json``, and prints a
short summary with the runtime.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

from spp.pipeline.l1b_pipeline import L1BPipeline
from spp.products.geotiff_writer import GeoTIFFWriter
from spp.qa.radiometric_validator import QAThresholds, RadiometricValidator
from spp.readers.package_reader import PackageReader


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="spp-l1b",
        description="Calibrate an acquisition package to L1B TOA radiance.",
    )
    parser.add_argument(
        "--input", "-i", required=True, type=Path,
        help="Path to the acquisition package directory.",
    )
    parser.add_argument(
        "--output", "-o", required=True, type=Path,
        help="Output directory for the L1B band rasters and QA report.",
    )
    parser.add_argument(
        "--bands", "-b", nargs="+", default=None,
        help="Subset of band names to process (default: all).",
    )
    parser.add_argument(
        "--window-lines", type=int, default=2048,
        help="Along-track lines processed per window (default: 2048).",
    )
    parser.add_argument(
        "--saturation-dn", type=int, default=None,
        help="DN at/above which a pixel is flagged as saturated (default: off).",
    )
    parser.add_argument(
        "--quicklook", action="store_true",
        help="Also write an RGB quicklook PNG (needs the R, G, B bands).",
    )
    parser.add_argument(
        "--stac", action=argparse.BooleanOptionalAction, default=True,
        help="Write a STAC item describing the product (default: on).",
    )
    parser.add_argument(
        "--quiet", "-q", action="store_true",
        help="Only print the final summary.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.WARNING if args.quiet else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )

    pipeline = L1BPipeline(
        reader=PackageReader(args.input),
        writer=GeoTIFFWriter(args.output),
        validator=RadiometricValidator(
            QAThresholds(saturation_dn=args.saturation_dn)
        ),
        bands=args.bands,
        window_lines=args.window_lines,
    )

    product = pipeline.run()

    qa = product.metadata["qa"]
    report_path = args.output / "qa_report.json"
    report_path.write_text(json.dumps(qa, indent=2))

    quicklook_path = _write_quicklook(product, args.output) if args.quicklook else None

    if args.stac:
        _write_stac(product, args.output, quicklook_path, report_path)

    _print_summary(product, report_path)
    return 0 if qa["passed"] else 1


def _write_quicklook(product, output_dir: Path) -> Path | None:
    from spp.qa.quicklook import write_quicklook

    required = ("R", "G", "B")
    if not all(b in product.bands for b in required):
        logging.getLogger(__name__).warning(
            "Skipping quicklook: needs bands %s", required
        )
        return None
    path = write_quicklook(product.bands, output_dir / "quicklook.png")
    logging.getLogger(__name__).info("Wrote quicklook %s", path)
    return path


def _write_stac(
    product, output_dir: Path, quicklook_path: Path | None, report_path: Path
) -> None:
    from spp.products.stac_writer import write_stac_item

    item_path = output_dir / f"{product.scene_id}_{product.level}.json"
    path = write_stac_item(
        product,
        item_path,
        quicklook_href=quicklook_path.name if quicklook_path else None,
        qa_href=report_path.name,
    )
    logging.getLogger(__name__).info("Wrote STAC item %s", path)


def _print_summary(product, report_path: Path) -> None:
    qa = product.metadata["qa"]
    runtime = product.metadata.get("runtime_seconds")
    print("\n=== L1B product ===")
    print(f"scene      : {product.scene_id}")
    print(f"level/units: {product.level}  [{product.units}]")
    print(f"bands      : {len(product.bands)}")
    print(f"runtime    : {runtime} s")
    print(f"QA passed  : {qa['passed']}")
    print(f"{'band':>5} {'mean':>10} {'min':>10} {'max':>10} {'nodata%':>9}  flags")
    for name, r in qa["bands"].items():
        print(
            f"{name:>5} {r['mean']:>10.3f} {r['min']:>10.3f} {r['max']:>10.3f} "
            f"{r['nodata_fraction'] * 100:>8.2f}%  {','.join(r['flags']) or '-'}"
        )
    print(f"\nwritten to : {report_path.parent}")
    print(f"qa report  : {report_path}")


if __name__ == "__main__":
    sys.exit(main())
