# Satellite Processing Pipeline

A modular Earth Observation (EO) processing framework designed to transform raw satellite payload acquisitions into analysis-ready products through configurable Level-0 to Level-2 processing chains.

The project focuses on building transparent, reproducible, and sensor-agnostic processing workflows that can be adapted to different satellite missions, calibration models, spectral configurations, and geospatial reference systems.

---

## Features

* Modular processing architecture
* Sensor-agnostic design
* Radiometric calibration workflows
* Product-level abstraction (L0 → L2)
* Quality Assurance (QA) framework
* STAC-compatible metadata generation
* Extensible processing pipeline
* Reproducible scientific workflows

---

## Processing Levels

### Level 0 (Raw Acquisition)

Raw payload data as received from the instrument.

Typical operations:

* Data ingestion
* Integrity checks
* Metadata extraction
* Packet reconstruction

Output:

```text
Raw acquisition package
```

---

### Level 1A

Instrument data organized into detector counts and acquisition structures.

Typical operations:

* Decompression
* Line reconstruction
* Detector organization
* Missing data detection

Output:

```text
Detector counts
```

---

### Level 1B

Radiometrically calibrated imagery.

Typical operations:

* Dark current correction
* Offset correction
* Gain correction
* Conversion from digital numbers to radiance
* QA flag generation

Output:

```text
Top-of-Atmosphere Radiance
```

---

### Level 1C

Geometrically corrected imagery.

Typical operations:

* Georeferencing
* Orthorectification
* Terrain correction
* Reprojection

Output:

```text
Georeferenced radiance product
```

---

### Level 2

Surface-derived geophysical products.

Typical operations:

* Atmospheric correction
* Reflectance generation
* Environmental retrievals

Output:

```text
Surface Reflectance
Derived Products
```

---

## Repository Structure

```text
satellite-processing-pipeline/

├── README.md
├── requirements.txt
├── pyproject.toml
│
├── configs/
│   ├── sensor_template.yaml
│   └── processing_levels.yaml
│
├── data/
│   ├── raw/
│   ├── calibration/
│   └── external/
│
├── docs/
│   ├── architecture.md
│   ├── product_hierarchy.md
│   ├── decision_log.md
│   ├── validation_strategy.md
│   └── references.md
│
├── outputs/
│   ├── l1b/
│   ├── qa/
│   └── quicklooks/
│
├── src/
│   └── spp/
│
│       ├── readers/
│       ├── calibration/
│       ├── geolocation/
│       ├── atmospheric/
│       ├── products/
│       ├── metadata/
│       ├── qa/
│       ├── pipeline/
│       └── utils/
│
├── scripts/
│
└── tests/
```

---

## Installation

Create a virtual environment:

```bash
python -m venv .venv
source .venv/bin/activate
```

Install dependencies:

```bash
pip install -r requirements.txt
```

---

## Core Dependencies

```text
numpy
rasterio
xarray
dask
pystac
matplotlib
pyyaml
pytest
```

Optional:

```text
gdal
py6s
```

---

## Example Workflow

Run a Level-1B calibration:

```bash
python scripts/run_l1b.py \
    --input data/raw/acquisition \
    --calibration data/calibration \
    --output outputs/l1b
```

Generate QA products:

```bash
python scripts/generate_quicklook.py
```

---

## Design Principles

### Sensor Agnostic

Mission-specific parameters should be externalized into configuration files.

Examples:

* Spectral bands
* Calibration coefficients
* Detector geometry
* Metadata mappings
* Coordinate reference systems

---

### Reproducibility

All processing steps should be:

* deterministic
* documented
* version controlled

---

### Transparency

Each processing level should clearly define:

* inputs
* outputs
* assumptions
* limitations
* validation strategy

---

## Quality Assurance

Recommended QA checks:

* Missing calibration coefficients
* Invalid radiance values
* Saturated pixels
* Histogram anomalies
* Metadata consistency checks

Outputs should include:

```text
QA report
Quicklook imagery
Processing log
```

---

## Future Work

Potential future extensions include:

* Orthorectification workflows
* Atmospheric correction
* BRDF correction
* STAC catalog generation
* Cloud-native GeoTIFF support
* Distributed processing with Dask
* Multi-sensor support

---

## References

* CEOS Product Levels
* STAC Specification
* GDAL
* Rasterio
* PySTAC
* SatPy
* Pygac
* Copernicus DEM
* Py6S

---

## License

MIT License

```
```
