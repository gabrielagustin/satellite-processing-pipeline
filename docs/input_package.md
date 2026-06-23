# Input Package Specification

The processing pipeline expects an acquisition package containing:

```text
acquisition_package/

├── B.tiff
├── G.tiff
├── R.tiff
├── RE1.tiff
├── RE2.tiff
├── RE3.tiff
├── NIR.tiff
├── PAN.tiff
│
├── <scene_id>.json          # STAC item (band asset table + scene identity)
├── metadata.json
├── ancillary.json
│
├── calibration/
│   ├── CPF_*.json              # calibration parameter file (matched by CPF*.json)
│   ├── filters_*.json          # spectral filter definitions (matched by filters*.json)
│   └── solar/
│       ├── spectral_solar.json
│       └── temporal_solar.json
│
├── README.md                # package documentation (optional)
└── thumbnail.webp           # preview image (optional)
```

Notes

* File names shown above are illustrative.
* Proprietary datasets are intentionally excluded from the repository.
* Example files are synthetic and only intended to document expected interfaces.

Band data are delivered as individual TIFF files.

Calibration and ancillary information are provided as JSON documents. The solar
reference (`spectral_solar.json`, `temporal_solar.json`) lives under
`calibration/solar/`; it is consumed by reflectance-level processing, not by L1B.

Inside `calibration/`, the reader locates the two documents by **glob pattern**,
not by an exact name, so payload-specific suffixes do not need to be hard-coded:

* the calibration parameter file (CPF) is matched by `CPF*.json`;
* the spectral filter definitions are matched by `filters*.json`.

The first match for each pattern is used (and a missing match raises a clear
error). The filter file provides the authoritative band-name↔detector-id map.

The STAC item at the package root carries the scene identity and the band asset
table, which is the reader's primary source for discovering the band rasters; it
falls back to glob patterns when the STAC item is absent.

The package reader is responsible for converting the package contents into an Acquisition domain object that can be consumed by downstream processing stages.
