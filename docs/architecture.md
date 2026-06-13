# Architecture

## Overview

Satellite Processing Pipeline (SPP) is a modular Earth Observation processing framework designed to transform raw payload acquisitions into higher-level geospatial products.

The framework follows a processing-level architecture inspired by Committee on Earth Observation Satellites (CEOS) product levels and common EO processing systems.

```text
Raw Acquisition
      │
      ▼
     L0
      │
      ▼
     L1A
      │
      ▼
     L1B
      │
      ▼
     L1C
      │
      ▼
      L2
```

The current implementation focuses on the L0 → L1B transition while maintaining an architecture that can support future processing levels.

---

# Design Goals

## Sensor Agnostic

Sensor-specific information must be externalized whenever possible.

Examples:

* spectral bands
* calibration coefficients
* detector geometry
* metadata mappings
* quality thresholds

The processing engine should not require code modifications when introducing a new sensor.

---

## Modular Processing

Each processing step should be independently testable.

```text
Reader
   ↓
Calibrator
   ↓
QA Validator
   ↓
Product Writer
```

Modules communicate through well-defined data structures.

---

## Reproducibility

Given the same inputs and configuration, the pipeline must always produce identical outputs.

---

## Transparency

Every product level should clearly define:

* inputs
* outputs
* assumptions
* validation criteria

---

# High-Level Components

## Reader Layer

Responsible for loading acquisitions and calibration assets.

### Responsibilities

* Parse raw acquisition files
* Load metadata
* Validate file integrity
* Expose data in a common format

### Interface

```python
class Reader:

    def load_acquisition(self):
        pass

    def load_metadata(self):
        pass
```

### Examples

```text
Raw Binary Reader
GeoTIFF Reader
NetCDF Reader
```

---

# Calibration Layer

Transforms instrument measurements into physically meaningful quantities.

### Responsibilities

* Dark current correction
* Offset correction
* Gain correction
* Radiometric conversion

### Interface

```python
class Calibrator:

    def apply(self, image):
        pass
```

### Current Focus

L1B radiometric calibration.

---

# Quality Assurance Layer

Evaluates product validity.

### Responsibilities

* Detect invalid values
* Detect saturation
* Verify metadata consistency
* Generate QA reports

### Interface

```python
class QAValidator:

    def validate(self, product):
        pass
```

### Outputs

```text
QA Report
Statistics
Warnings
Quicklooks
```

---

# Product Layer

Responsible for generating output products.

### Responsibilities

* Write raster outputs
* Generate metadata
* Apply naming conventions

### Interface

```python
class ProductWriter:

    def write(self, product):
        pass
```

---

# Metadata Layer

Responsible for metadata generation and catalog integration.

### Responsibilities

* Product metadata
* Processing lineage
* STAC generation

### Interface

```python
class MetadataBuilder:

    def build(self):
        pass
```

---

# Pipeline Layer

Coordinates all processing steps.

### Responsibilities

* Execute processing workflow
* Track processing status
* Generate logs

### Interface

```python
class ProcessingPipeline:

    def run(self):
        pass
```

---

# Data Model

## Acquisition

Represents raw instrument data.

```python
@dataclass
class Acquisition:

    image: np.ndarray
    metadata: dict
```

---

## CalibrationModel

Represents sensor calibration information.

```python
@dataclass
class CalibrationModel:

    gain: float
    offset: float
```

Future implementations may support:

```python
gain_per_band
offset_per_band
dark_reference
flat_field
spectral_response
```

---

## Product

Represents a generated processing product.

```python
@dataclass
class Product:

    level: str
    image: np.ndarray
    metadata: dict
```

---

# Processing Flow

## L0 → L1A

Input:

```text
Raw acquisition
```

Output:

```text
Detector counts
```

Responsibilities:

* unpacking
* ordering
* integrity checks

---

## L1A → L1B

Input:

```text
Detector counts
Calibration assets
```

Output:

```text
Radiance product
```

Responsibilities:

1. Dark current correction
2. Offset correction
3. Gain correction
4. Radiance conversion
5. QA generation

---

## Future L1B → L1C

Input:

```text
Radiance product
DEM
Attitude data
Ephemeris
```

Output:

```text
Georeferenced radiance
```

Candidate libraries:

* GDAL
* Rasterio

---

## Future L1C → L2

Input:

```text
Georeferenced radiance
Atmospheric parameters
```

Output:

```text
Surface reflectance
```

Candidate libraries:

* Py6S
* Acolite

---

# Configuration Strategy

All mission-specific parameters should live outside the processing code.

Example:

```yaml
sensor:
  name: example_sensor

bands:
  - blue
  - green
  - red
  - nir

calibration:
  gain: 0.01
  offset: 12
```

This allows new missions to be supported through configuration rather than code changes.

---

# Testing Strategy

Unit Tests

* Reader tests
* Calibration tests
* QA tests
* Metadata tests

Integration Tests

* End-to-end L0 → L1B

Validation Tests

* Radiometric sanity checks
* Histogram verification
* Product consistency checks

---

# Future Extensions

Potential roadmap:

Phase 1

* L1B implementation

Phase 2

* STAC catalog support

Phase 3

* L1C orthorectification

Phase 4

* L2 atmospheric correction

Phase 5

* Multi-sensor support

Phase 6

* Distributed processing with Dask

---

# Architectural References

The architecture is inspired by concepts commonly used in operational EO processing systems, including:

* CEOS Product Levels
* STAC
* SatPy
* Pygac
* GDAL ecosystem

The implementation intentionally favors lightweight, transparent, and sensor-agnostic components over large mission-specific frameworks.
