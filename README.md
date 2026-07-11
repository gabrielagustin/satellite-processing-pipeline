# Satellite Processing Pipeline

[![CI](https://github.com/gabrielagustin/satellite-processing-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/gabrielagustin/satellite-processing-pipeline/actions/workflows/ci.yml)

A modular, sensor-agnostic Earth Observation processing framework that turns raw
multispectral payload acquisitions into calibrated products through a transparent
L0 → L2 processing chain.

The framework implements the **L1B transition end-to-end** — radiometric
calibration of raw digital numbers (DN) to top-of-atmosphere (TOA) spectral
radiance — and specifies the remaining levels (L1A, L1C, L2A) in detail. It is
designed to read as a reusable EO tool: sensor-specific knowledge lives in the
data it discovers and in a swappable calibration model, not hard-coded in the
processing engine.

---

## Processing levels

| Level | Output | Status |
|------|--------|--------|
| L0 / L1A | Per-band DN rasters (sensor coords) | Ingested |
| **L1B** | **TOA spectral radiance** | **Implemented** |
| L1C | Georeferenced / orthorectified radiance | **Partly implemented** — georeferencing and orthorectification work; band co-registration does not yet. See [`docs/l1c_spec.md`](docs/l1c_spec.md) and [`docs/limitations.md`](docs/limitations.md) |
| L2A | Surface reflectance | Specified |

See [`docs/product_hierarchy.md`](docs/product_hierarchy.md) for the full
hierarchy and [`docs/remaining_levels.md`](docs/remaining_levels.md) for the
levels not yet implemented. L1C — georeferencing, orthorectification and band
co-registration — has a detailed implementable design in
[`docs/l1c_spec.md`](docs/l1c_spec.md).

---

## Installation

Requires Python 3.11+.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .            # runtime deps + the `spp-l1b` command
```

Core dependencies: `numpy`, `rasterio`, `matplotlib` (the RGB quicklook is on by
default). For the test suite, `pip install -e '.[test]'`.

Dependencies are declared in `pyproject.toml` (the packaging source of truth,
which also registers the `spp-l1b` command and the `test` extra). A
`requirements.txt` mirroring the runtime deps is also provided for convenience
and for environments that expect one (`pip install -r requirements.txt`);
installing the package with `pip install -e .` alone is sufficient.

---

## How to run

Calibrate an acquisition package to L1B radiance:

```bash
spp-l1b --input  /path/to/acquisition_package \
        --output /path/to/output_dir
```

This writes the band rasters, a QA report, a STAC item and an RGB quicklook
(disable the last two with `--no-stac` / `--no-quicklook`).

Or without installing:

```bash
PYTHONPATH=src python -m spp.cli.run_l1b --input <package> --output <out>
```

### L1C — georeference, orthorectify and co-register

Takes the acquisition package (for the platform telemetry and the calibration) and the
**L1B rasters** (for the pixels):

```bash
spp-l1c --input  /path/to/acquisition_package \
        --l1b    /path/to/l1b_output \
        --output /path/to/l1c_output
```

The elevation model is fetched automatically from Copernicus GLO-30 over the network —
no account, no manual download. Supply your own with `--dem <raster-or-vrt>`, or skip
terrain entirely with `--no-dem` (which georeferences but does **not** orthorectify, and
says so).

Writes one projected `float32` radiance raster per band on a common UTM grid, plus
`qa_report_l1c.json`. A full 8-band run takes roughly 7 minutes and needs ~2 GB of RAM.
Try two bands first:

```bash
spp-l1c --input <package> --l1b <l1b_out> --output <l1c_out> --bands PAN RE1
```

**Read the run's closing summary.** It ends with a section headed *NOT ESTABLISHED BY
THIS RUN*, which lists what the product does **not** demonstrate — the absolute
geolocation is unvalidated against anything independent of the telemetry, two telemetry
conventions are assumed rather than resolved, and some bands may have been refused by the
self-calibration. The rasters look finished whether or not any of that is true, which is
why the run says it out loud. See [`docs/limitations.md`](docs/limitations.md).

| Flag | Meaning |
|---|---|
| `--bands B G R ...` | Subset of bands (default: all) |
| `--dem <path>` | Elevation raster or virtual mosaic (default: fetch Copernicus GLO-30) |
| `--no-dem` | Georeference on the ellipsoid; do **not** orthorectify |
| `--gsd N` | Output resolution in metres (default: the native sampling, rounded up) |
| `--reference-band B` | Band the others co-register onto (default: `PAN`) |
| `--no-refine` | Skip the self-calibration; the bands will **not** be co-registered |
| `--step N` | Geolocation-lattice spacing in pixels (default 8, matched to the DEM) |

### Expected input

An acquisition package directory containing:

- one GeoTIFF per spectral band (`uint16` DN),
- `metadata.json` (imager configuration + detector telemetry),
- `ancillary.json` (platform telemetry — used by later levels),
- a STAC item describing the product,
- `calibration/` with the Calibration Parameter File, spectral filters and the
  solar reference.

See [`docs/input_package.md`](docs/input_package.md) for the layout. Band files
and calibration assets are **discovered** (via the STAC assets and glob
patterns), so exact file names are not hard-coded.

### Outputs

Written to the output directory:

- `<BAND>.tif` — one `float32` TOA radiance raster per band
  (`W / (m² · sr · µm)`), tiled with internal overviews, NoData = NaN, no CRS
  (L1B is still in sensor coordinates);
- `qa_report.json` — per-band statistics and quality flags;
- `<scene_id>_L1B.json` — a STAC item cataloguing the product (band assets,
  spectral/raster properties, processing lineage; `geometry` is `null` because
  L1B is not yet georeferenced). On by default; disable with `--no-stac`;
- `quicklook.png` — RGB preview. On by default; disable with `--no-quicklook`.

### Options

| Flag | Meaning |
|---|---|
| `--bands B G R ...` | Process a subset of bands (default: all) |
| `--window-lines N` | Along-track lines per processing window (default 2048) |
| `--saturation-dn N` | Flag DN ≥ N as saturated (default: off) |
| `--quicklook` / `--no-quicklook` | Write an RGB quicklook PNG (default: on) |
| `--stac` / `--no-stac` | Write a STAC item for the product (default: on) |
| `--quiet` | Only print the final summary |

### Approximate runtime

A full 8-band acquisition of ~4096 × 30948 pixels processes in **~2 minutes** on
a laptop (≈13–18 s per band), producing ~4 GB of `float32` output. Memory stays
flat (windowed streaming), independent of raster size.

---

## Tests

```bash
pip install -e '.[test]'
pytest
```

The suite is organised one file per module under test: `test_l1b_calibrator.py`
covers the L1B calibrator core (the `calibrate()` DN→radiance path, windowed
exactness, NoData masking, validation) and the temperature-interpolation model;
`test_package_reader.py` covers the reader timing/parsing helpers (timestamp→line
mapping, optional-field coercion). Synthetic fixtures only — no proprietary data
required.

---

## Repository structure

```text
satellite-processing-pipeline/
├── src/spp/
│   ├── core/          # domain entities (Band, Acquisition, Calibration, Product)
│   ├── readers/       # package discovery and parsing -> Acquisition
│   ├── calibration/   # radiometric calibration (L1B)
│   ├── geometry/      # sensor model, frames, terrain (L1C)
│   ├── resample/      # target grid + geolocation-array warp (L1C)
│   ├── qa/            # quality assessment + quicklook
│   ├── products/      # raster writers (GeoTIFF)
│   ├── pipeline/      # orchestration (L1B pipeline)
│   └── cli/           # command-line entry points
├── docs/              # product hierarchy, decisions, specs, validation
├── pyproject.toml
└── requirements.txt
```

---

## Documentation

| Document | Contents |
|---|---|
| [architecture.md](docs/architecture.md) | Component design and data flow |
| [product_hierarchy.md](docs/product_hierarchy.md) | L0→L2 levels: inputs, outputs, justification |
| [remaining_levels.md](docs/remaining_levels.md) | Specs for L1A, L1C, L2A |
| [l1c_spec.md](docs/l1c_spec.md) | L1C geometric design: sensor model, orthorectification, band co-registration |
| [l1c_findings.md](docs/l1c_findings.md) | What measurement overturned: five assumptions, and how each was caught |
| [decision_log.md](docs/decision_log.md) | Engineering decisions and trade-offs |
| [limitations.md](docs/limitations.md) | Failure modes and next steps |
| [validation_strategy.md](docs/validation_strategy.md) | How correctness is established |
| [qa_report.md](docs/qa_report.md) | `qa_report.json` schema: metrics, flags, thresholds |
| [generalisation.md](docs/generalisation.md) | Adapting to a new mission |
| [performance.md](docs/performance.md) | Processing large strips: memory & performance |
| [input_package.md](docs/input_package.md) | Expected input package layout |
| [references.md](docs/references.md) | Standards and tooling references |

---

## Design principles

- **Separation of concerns** — reader, calibrator, QA, writer, pipeline are
  independent and individually testable.
- **Sensor-agnostic core** — entities carry no mission/sensor constants; the
  radiometric model is swappable behind an interface.
- **Streaming I/O** — large rasters are processed in windows; memory stays flat.
- **Reproducibility** — deterministic outputs; every run emits a QA report.

---

## License

MIT License.
