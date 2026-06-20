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
├── metadata.json
├── ancillary.json
│
├── calibration/
│   ├── calibration_parameters.json
│   └── filter_definitions.json
│
└── solar/
    ├── spectral_solar.json
    └── temporal_solar.json
```

Notes

* File names shown above are illustrative.
* Proprietary datasets are intentionally excluded from the repository.
* Example files are synthetic and only intended to document expected interfaces.

Band data are delivered as individual TIFF files.

Calibration and ancillary information are provided as JSON documents.

The package reader is responsible for converting the package contents into an Acquisition domain object that can be consumed by downstream processing stages.
