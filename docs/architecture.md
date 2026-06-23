# Architecture

The Satellite Processing Pipeline (SPP) is a modular, sensor-agnostic Earth
Observation framework. It transforms raw payload acquisitions into calibrated
products through a chain of small, single-responsibility stages. This document
describes the **implemented L1B architecture** and the interfaces it is built on.

For the product levels themselves see
[`product_hierarchy.md`](product_hierarchy.md); for how large strips are
processed efficiently see [`performance.md`](performance.md).

---

## Data flow

```text
acquisition package (filesystem)
        │
        ▼
   PackageReader.read()  ──►  Acquisition  (lazy bands + calibration + config)
        │
        ▼
   L1BPipeline.run()
        │   for each band, streaming windows of lines:
        │       read DN window
        │       ─► L1BCalibrator.calibrate(band, dn, line_start)  ─► radiance
        │       ─► GeoTIFFWriter.write_block(radiance, line_start)
        │       ─► RadiometricValidator accumulator.update(radiance, dn)
        ▼
   Product  +  qa_report.json  +  STAC item  (+ quicklook.png)
```

Stages communicate through the core domain entities, never through ad-hoc dicts.

---

## Layers and interfaces

Each stage is an abstract base class (the seam) with a concrete implementation.
The signatures below are the real ones.

### Reader — `readers/`

```python
class Reader(ABC):
    def read(self) -> Acquisition: ...
```

`PackageReader` discovers the band rasters (from the STAC assets) and the
calibration assets (by glob), reads only raster **metadata** (lazy bands), parses
the CPF and imager configuration, and assembles an `Acquisition`. No calibration,
no processing.

### Calibrator — `calibration/`

```python
class Calibrator(ABC):
    def calibrate(self, band_name: str, dn: np.ndarray,
                  line_start: int = 0) -> np.ndarray: ...
```

`L1BCalibrator` converts a **block** of DN to at-aperture radiance. It is pure
(NumPy in, NumPy out, no I/O); `line_start` lets along-track-varying terms (the
temperature profile) be indexed correctly, which is what makes windowed
processing exact.

### Quality validator — `qa/`

```python
class QAValidator(ABC):
    def validate(self, product) -> dict: ...
```

`RadiometricValidator` provides a streaming `accumulator(band)` fed per window,
plus a `validate(product)` that re-reads a finished product. It only measures
(non-destructive) and emits the per-band metrics/flags documented in
[`qa_report.md`](qa_report.md). `write_quicklook()` produces an RGB preview.

### Product writer — `products/`

```python
class ProductWriter(ABC):
    def open_band(self, name, *, width, height, dtype, nodata): ...  # context manager
```

`GeoTIFFWriter` opens a band as a context-managed handle with
`write_block(array, line_start)`, writing a tiled, compressed GeoTIFF with
overviews and no CRS (L1B is in sensor coordinates).

### Pipeline — `pipeline/`

```python
class ProcessingPipeline(ABC):
    def run(self) -> Product: ...
```

`L1BPipeline` wires reader → calibrator → QA → writer and owns **only**
orchestration and the windowing loop. Collaborators are injected, so each can be
swapped or tested in isolation.

### CLI — `cli/`

`run_l1b` (`spp-l1b`) parses arguments, builds the pipeline, writes the band
rasters + `qa_report.json` (+ optional quicklook) and returns an exit code driven
by the QA result.

---

## Core data model — `core/`

Pure dataclasses, no raster-library dependency, no sensor proper nouns.

| Entity | Role |
|---|---|
| `Band` | One spectral band — **lazy**: path + raster properties, never pixels. |
| `Acquisition` | The package as a domain object: bands + calibration + imager config + telemetry. |
| `ImagerConfiguration` | Per-acquisition imager settings (line period, per-band TDI/start row/CWL — central wavelength). |
| `RadiometricCalibration` | Coefficients for one `(band, start_row, tdi)`. |
| `CalibrationParameters` | Calibration collection with `lookup(band, start_row, tdi)`; retains the geometric block for L1C. |
| `Product` | Generated product described **by reference** (output paths) + metadata + provenance. |

The physics lives in the calibrator, not the entities — the core stays a pure
data layer.

---

## Design principles

- **Separation of concerns** — reader, calibrator, QA, writer, pipeline are
  independent and individually testable.
- **Sensor-agnostic core** — entities carry no mission/sensor constants; the
  radiometric model is swappable behind the `Calibrator` interface (see
  [`generalisation.md`](generalisation.md)).
- **Streaming I/O** — large rasters are processed in windows; memory stays flat
  (see [`performance.md`](performance.md)).
- **Reproducibility** — deterministic outputs; every run emits a QA report.
- **No premature abstraction** — only the seams that earn their keep today exist;
  future layers (geolocation, atmospheric, catalogue) are specified, not stubbed.

---

## Future stages

`L1C` (geolocation) and `L2A` (surface reflectance) are specified in
[`remaining_levels.md`](remaining_levels.md). They slot in as additional
`Calibrator`/pipeline implementations consuming the same `Acquisition` plus their
external data (DEM, atmospheric state); the reader already retains the geometric
calibration and ancillary telemetry they need.
