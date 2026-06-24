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
| L1C | Georeferenced radiance | Specified |
| L2A | Surface reflectance | Specified |

See [`docs/product_hierarchy.md`](docs/product_hierarchy.md) for the full
hierarchy and [`docs/remaining_levels.md`](docs/remaining_levels.md) for the
levels not yet implemented.

---

## Installation

Requires Python 3.11+.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .            # runtime deps + the `spp-l1b` command
```

Core dependencies: `numpy`, `rasterio`. Quicklook generation additionally needs
`matplotlib` (`pip install -e '.[viz]'`); for the test suite,
`pip install -e '.[test]'`.

Dependencies are declared in `pyproject.toml` (the packaging source of truth,
which also registers the `spp-l1b` command and the `viz` / `test` extras). A
`requirements.txt` mirroring the runtime deps is also provided for convenience
and for environments that expect one (`pip install -r requirements.txt`);
installing the package with `pip install -e .` alone is sufficient.

---

## How to run

Calibrate an acquisition package to L1B radiance:

```bash
spp-l1b --input  /path/to/acquisition_package \
        --output /path/to/output_dir \
        --quicklook
```

Or without installing:

```bash
PYTHONPATH=src python -m spp.cli.run_l1b --input <package> --output <out>
```

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
- `quicklook.png` — RGB preview (with `--quicklook`).

### Options

| Flag | Meaning |
|---|---|
| `--bands B G R ...` | Process a subset of bands (default: all) |
| `--window-lines N` | Along-track lines per processing window (default 2048) |
| `--saturation-dn N` | Flag DN ≥ N as saturated (default: off) |
| `--quicklook` | Also write an RGB quicklook PNG |
| `--stac` / `--no-stac` | Write a STAC item for the product (default: on) |
| `--quiet` | Only print the final summary |

### Approximate runtime

A full 8-band acquisition of ~4096 × 31000 pixels processes in **~1–2 minutes** on
a laptop, producing ~4 GB of `float32` output. Memory stays flat (windowed
streaming), independent of raster size.

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
