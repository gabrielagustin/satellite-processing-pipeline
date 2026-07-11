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
| L1C | Georeferenced / orthorectified radiance | **Partly implemented** — georeferencing and orthorectification work; band co-registration does not yet. See [`docs/l1c_spec.md`](docs/l1c/spec.md) and [`docs/limitations.md`](docs/l1b/limitations.md) |
| L2A | Surface reflectance | Specified |

See [`docs/product_hierarchy.md`](docs/product_hierarchy.md) for the full
hierarchy and [`docs/remaining_levels.md`](docs/remaining_levels.md) for the
levels not yet implemented. L1C — georeferencing, orthorectification and band
co-registration — has a detailed implementable design in
[`docs/l1c_spec.md`](docs/l1c/spec.md).

---

## Installation

Requires Python 3.11+.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .            # runtime deps + the `spp-l1b` and `spp-l1c` commands
```

Core dependencies: `numpy`, `rasterio`, `matplotlib` (the RGB quicklook is on by
default) and `pyproj` (coordinate systems and the geoid, used by L1C). For the test
suite, `pip install -e '.[test]'`.

L1C fetches its elevation model and geoid grid over the network on first use; both are
cached, and `--dem` accepts a local one.

Dependencies are declared in `pyproject.toml` (the packaging source of truth,
which also registers the `spp-l1b` command and the `test` extra). A
`requirements.txt` mirroring the runtime deps is also provided for convenience
and for environments that expect one (`pip install -r requirements.txt`);
installing the package with `pip install -e .` alone is sufficient.

---

## L1B — calibrate to radiance

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

### L1B — expected input

An acquisition package directory containing:

- one GeoTIFF per spectral band (`uint16` DN),
- `metadata.json` (imager configuration + detector telemetry),
- `ancillary.json` (platform telemetry — used by later levels),
- a STAC item describing the product,
- `calibration/` with the Calibration Parameter File, spectral filters and the
  solar reference.

See [`docs/input_package.md`](docs/l1b/input_package.md) for the layout. Band files
and calibration assets are **discovered** (via the STAC assets and glob
patterns), so exact file names are not hard-coded.

### L1B — outputs

Written to the output directory:

- `<BAND>.tif` — one `float32` TOA radiance raster per band
  (`W / (m² · sr · µm)`), tiled with internal overviews, NoData = NaN, no CRS
  (L1B is still in sensor coordinates);
- `qa_report.json` — per-band statistics and quality flags;
- `<scene_id>_L1B.json` — a STAC item cataloguing the product (band assets,
  spectral/raster properties, processing lineage; `geometry` is `null` because
  L1B is not yet georeferenced). On by default; disable with `--no-stac`;
- `quicklook.png` — RGB preview. On by default; disable with `--no-quicklook`.

### L1B — options

| Flag | Meaning |
|---|---|
| `--bands B G R ...` | Process a subset of bands (default: all) |
| `--window-lines N` | Along-track lines per processing window (default 2048) |
| `--saturation-dn N` | Flag DN ≥ N as saturated (default: off) |
| `--quicklook` / `--no-quicklook` | Write an RGB quicklook PNG (default: on) |
| `--stac` / `--no-stac` | Write a STAC item for the product (default: on) |
| `--quiet` | Only print the final summary |

### L1B — approximate runtime

A full 8-band acquisition of ~4096 × 30948 pixels processes in **~2 minutes** on
a laptop (≈13–18 s per band), producing ~4 GB of `float32` output. Memory stays
flat (windowed streaming), independent of raster size.

---

## L1C — georeference, orthorectify, co-register

Takes the **acquisition package** (for the platform telemetry and the calibration) *and*
the **L1B rasters** (for the pixels). Both are required: the L1B GeoTIFFs carry only the
radiance, not the ephemeris, attitude, per-line timestamps or camera geometry the
geometric level needs.

```bash
spp-l1c --input  /path/to/acquisition_package \
        --l1b    /path/to/l1b_output \
        --output /path/to/l1c_output
```

The elevation model is fetched automatically from Copernicus GLO-30 over the network —
no account, no manual download — and cached in the output directory. Supply your own with
`--dem <raster-or-vrt>`, or skip terrain with `--no-dem` (which georeferences but does
**not** orthorectify, and says so).

### L1C — the reference band is mandatory

The bands are co-registered **onto a reference band**, `PAN` by default. If you process a
subset, **the reference band must be in it**:

```bash
spp-l1c ... --bands PAN R G B          # correct: PAN included
spp-l1c ... --bands R G B              # runs, but the bands will NOT co-register
```

Without it the run still produces georeferenced, orthorectified rasters — it simply skips
the self-calibration, raises `coregistration_not_achieved`, and says so in the summary.
Override the reference with `--reference-band <BAND>`.

### L1C — outputs

Written to the output directory:

- **`stack.tif`** — **the product.** One multi-band `float32` raster: every band on the
  same projected grid, co-registered, with band names set. This is what to open.
- `<BAND>.tif` — the same bands as individual rasters, for convenience.
- `qa_report_l1c.json` — conventions used, grid, terrain, interpolation error,
  per-band co-registration, the estimated line-of-sight calibration, and the flags.
- `cache/dem.vrt` — the elevation mosaic that was fetched. Reuse it across runs with
  `--dem`.

### L1C — options

| Flag | Meaning |
|---|---|
| `--bands B G R ...` | Subset of bands (default: all). **Must include the reference band** |
| `--reference-band B` | Band the others co-register onto (default: `PAN`) |
| `--dem <path>` | Elevation raster or virtual mosaic (default: fetch Copernicus GLO-30) |
| `--no-dem` | Georeference on the ellipsoid; do **not** orthorectify |
| `--gsd N` | Output resolution in metres (default: the native sampling, rounded up) |
| `--no-refine` | Skip the self-calibration of the missing line-of-sight term |
| `--no-scene-correction` | Skip the scene-local (attitude) correction |
| `--no-stack` | Do not write `stack.tif` (per-band rasters only) |
| `--step N` | Geolocation-lattice spacing in pixels (default 8, matched to the DEM) |
| `--quiet` | Only print the final summary |

### L1C — approximate runtime

A full 8-band run takes roughly **10 minutes** and needs **~2 GB of RAM**, producing a
~3.5 GB stack. Four bands take about 4 minutes. The warp is **not** memory-bounded (unlike
L1B) — see [`docs/performance.md`](docs/l1b/performance.md).

### L1C — read the closing summary

The run ends with a section headed **NOT ESTABLISHED BY THIS RUN**. It lists what the
product does *not* demonstrate: that the absolute geolocation has never been checked
against anything independent of the telemetry that produced it, that two telemetry
conventions are mirrors geometry cannot see and are therefore assumed, and which bands the
self-calibration refused.

The rasters open in any geographic information system and overlay a basemap plausibly
whether or not any of that is true, which is exactly why the run says it out loud. Band
co-registration currently reaches **2–3 px**, not the sub-pixel target; the floor is
attitude jitter. See [`docs/limitations.md`](docs/l1b/limitations.md).

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

Organised by processing level, with the cross-cutting documents at the root.

| | |
|---|---|
| **[docs/l1b/](docs/l1b/)** | **L1B — radiometric calibration.** Architecture, input package, validation, QA schema, performance, limitations |
| **[docs/l1c/](docs/l1c/)** | **L1C — geometric processing.** [Spec](docs/l1c/spec.md), [findings](docs/l1c/findings.md), architecture, validation, QA schema, performance, [limitations](docs/l1c/limitations.md) |

| Cross-cutting | Contents |
|---|---|
| [product_hierarchy.md](docs/product_hierarchy.md) | L0→L2 levels: inputs, outputs, justification |
| [decision_log.md](docs/decision_log.md) | Engineering decisions and trade-offs, in order |
| [remaining_levels.md](docs/remaining_levels.md) | Specs for the levels not yet built (L1A, L2A) |
| [generalisation.md](docs/generalisation.md) | Adapting the framework to a new mission |
| [references.md](docs/references.md) | Standards and tooling references |

**[`docs/l1c/findings.md`](docs/l1c/findings.md)** is the one to read if you read only one:
six of the geometric specification's load-bearing assumptions were wrong, four were
indistinguishable from an irreducible platform error, and every one was caught by a check
that could have failed.

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
