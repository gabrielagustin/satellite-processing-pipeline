# Product Hierarchy

This document specifies the Earth Observation product levels this framework
targets, from raw acquisition (L0) up to surface products (L2). The hierarchy
follows the CEOS (Committee on Earth Observation Satellites) processing-level
conventions and is written for a generic
multispectral **pushbroom (linescan)** imager with the following characteristics
(the values are those of the reference acquisition used to develop the
framework):

- **8 spectral bands** — one panchromatic plus blue, green, red, three red-edge
  and a near-infrared band (VNIR — Visible and Near Infrared, ~490–842 nm).
- A single-frame **CIS (CMOS Image Sensor) detector** read out in linescan mode, with per-band
  **Time-Delay Integration (TDI)** and a per-band detector **start row**.
- Cross-track width of 4096 detector columns; along-track dimension built up over
  time (one raster row per readout line).
- Delivered quantisation: `uint16` digital numbers (DN).

The framework currently implements the **L1B** transition end-to-end
(radiometric calibration to top-of-atmosphere radiance). The remaining levels are
specified in [`remaining_levels.md`](remaining_levels.md).

---

## Overview

| Level | Name | Input | Output | Domain |
|------|------|-------|--------|--------|
| **L0** | Raw acquisition | Session binary (packets) | Raw DN, unstructured | Instrument |
| **L1A** | Reconstructed counts | L0 + session metadata | Per-band DN rasters (sensor coords) | Instrument |
| **L1B** | TOA radiance | L1A + radiometric calibration | Per-band radiance (sensor coords) | Physical |
| **L1C** | Georeferenced radiance | L1B + ephemeris/attitude + DEM | Radiance on a map grid (CRS) | Geometric |
| **L2A** | Surface reflectance | L1C + solar + atmospheric model | Surface reflectance | Geophysical |

> **A note on the delivered package.** The acquisition used here is labelled
> "L0", but it has already been unpacked into per-band DN rasters — i.e. it is
> effectively delivered at the **L0/L1A boundary**. True packet-level
> reconstruction (rebuilding the per-band rasters from the raw session binary) is
> therefore out of scope; the framework starts from the per-band rasters and the
> calibration assets. See the L0/L1A sections below for what that step entails.

---

## L0 — Raw acquisition

**Definition.** The data exactly as downlinked from the payload: a session
binary containing the detector readout stream plus housekeeping/telemetry, with
no image structure imposed.

**Input.** Raw session binary (`raw.bin`) and the downlink metadata.

**Output.** Byte-faithful raw packets in archival storage; no per-band imagery.

**Processing steps.**
1. Demultiplex the downlink stream; verify frame/packet integrity (CRC — cyclic
   redundancy check).
2. Extract session metadata (imager configuration, telemetry, timing).
3. Persist the raw binary and metadata unchanged for traceability.

**Justification.** L0 is the archival, loss-free record. Keeping it untouched
guarantees every higher product can be regenerated and audited from source.

---

## L1A — Reconstructed detector counts

**Definition.** The raw stream organised into one DN raster per spectral band, in
**sensor (detector) coordinates** — still raw counts, no radiometry applied.

**Input.** L0 session binary + `ImagerConfiguration` (line period, per-band TDI,
per-band start row, scan direction, binning).

**Output.** One `uint16` raster per band, `columns ∈ [0, 4096)` (cross-track) ×
`lines` (along-track); NoData = 0; no CRS.

**Processing steps.**
1. Decompress the session payload (if compressed).
2. Reconstruct readout lines in acquisition order using the line period and scan
   direction.
3. For each band, extract its detector window (start row, TDI height) and
   accumulate TDI lines into one output line per band.
4. Flag missing/dropped lines and dead pixels.

**Justification.** Separating "image assembly" from "radiometry" keeps the
instrument-geometry concerns (TDI, start rows, line timing) isolated from the
physics of calibration. The per-band raster is the natural unit every downstream
level consumes.

> The provided package is delivered **at this stage** (per-band DN rasters), so
> the framework's reader ingests L1A-equivalent rasters directly.

---

## L1B — Top-of-atmosphere radiance  *(implemented)*

**Definition.** Per-band, radiometrically calibrated **at-aperture spectral
radiance**, still in sensor coordinates.

**Input.**
- L1A per-band DN rasters.
- Radiometric calibration (CPF — Calibration Parameter File): per
  `(band, start_row, tdi)` an absolute scale,
  an optional offset, and three per-column arrays — non-uniformity,
  thermal-intercept and thermal-gradient.
- Per-line detector temperature telemetry.

**Output.** One `float32` radiance raster per band (units `W / (m² · sr · µm)`),
NoData = NaN, no CRS; a per-band QA report; and a **STAC item** cataloguing the
product (band assets, spectral/raster properties, processing lineage). The STAC
item's geometry is `null` — L1B is not yet georeferenced.

**Processing steps** (per band, applied per detector column `c` and line `ℓ`):
1. **Temperature-dependent dark subtraction.** Build a per-line temperature
   profile `T(ℓ)` by interpolating the telemetry, then
   `darkfield(c, ℓ) = thermal_intercept[c] + thermal_gradient[c] · T(ℓ)` and
   subtract it. This removes the dark-signal pedestal, which is present because
   on-detector electronic-black correction is disabled, and which drifts with
   detector temperature over the acquisition.
2. **Non-uniformity (relative) correction.** Multiply by `non_uniformity[c]` to
   equalise the pixel-to-pixel sensitivity (flat-field) across the detector.
3. **Absolute calibration.** `radiance = absolute · (relative_corrected) +
   absolute_offset`, converting corrected counts to physical radiance.

**Justification.** Radiance is the first physically meaningful, sensor-model-free
quantity: it expresses what the aperture measured in SI units, independent of
detector idiosyncrasies. It is the prerequisite for both geometric processing
(L1C) and reflectance retrieval (L2). Dark/flat/absolute is the standard
radiometric chain; doing it before geometry avoids resampling artefacts
contaminating the calibration.

**Validation.** Radiance magnitudes and inter-band ordering compared against
expectation (decreasing from visible to NIR over the scene); per-band statistics
and quality flags (NoData fraction, negative-radiance fraction, saturation,
dynamic range) recorded in the QA report. See
[`validation_strategy.md`](validation_strategy.md).

---

## L1C — Georeferenced (orthorectified) radiance

**Definition.** L1B radiance resampled onto a map grid in a defined CRS, with
bands co-registered to a common ground geometry.

**Input.** L1B radiance; platform ephemeris (position/velocity) and attitude
(quaternions); camera intrinsics; boresight alignment and per-band line-of-sight
arrays; a Digital Elevation Model (DEM).

**Output.** Per-band radiance on a projected grid (e.g. UTM/EPSG), with a valid
geotransform and CRS; band co-registration to sub-pixel accuracy.

**Processing steps (summary).** Build the per-pixel viewing rays from the
line-of-sight model + boresight + attitude; intersect them with the DEM-bearing
ellipsoid to get ground coordinates; construct the sensor-to-ground model and
resample each band onto the output grid. Detailed in
[`remaining_levels.md`](remaining_levels.md).

**Justification.** Geometric correction is what makes the imagery measurable in
ground space and stackable with other geodata. It depends on **external data**
(DEM) and the full platform geometry, so it is kept separate from radiometry.

---

## L2A — Surface reflectance

**Definition.** Atmospherically corrected **surface reflectance** (bottom-of-
atmosphere), the analysis-ready geophysical product.

**Input.** L1C radiance; band-averaged solar irradiance (ESUN — exo-atmospheric
solar irradiance) from the solar
reference and the per-band spectral response; sun/view geometry; an atmospheric
model (aerosol optical depth, water vapour, ozone); DEM.

**Output.** Per-band surface reflectance (unitless, 0–1) on the L1C grid; optional
water-leaving reflectance variant for aquatic scenes.

**Processing steps (summary).** Convert radiance to TOA reflectance using ESUN
and sun geometry; run an atmospheric radiative-transfer correction (e.g. 6S or an
aquatic-specific scheme) to remove scattering/absorption and retrieve surface
reflectance. Detailed in [`remaining_levels.md`](remaining_levels.md).

**Justification.** Surface reflectance removes illumination and atmosphere so
scenes are comparable across time and sensors — the basis for indices and
quantitative retrievals. It is the most external-data-intensive level.

---

## References

- CEOS — Earth Observation product level definitions: <https://ceos.org>
- STAC Specification — product/catalogue metadata: <https://stacspec.org>
- Rasterio / GDAL — raster I/O and warping: <https://rasterio.readthedocs.io>, <https://gdal.org>
- Copernicus DEM — elevation reference for L1C: <https://spacedata.copernicus.eu>
- 6S / Py6S — atmospheric radiative transfer for L2: <https://py6s.readthedocs.io>

See [`references.md`](references.md) for the full list.
