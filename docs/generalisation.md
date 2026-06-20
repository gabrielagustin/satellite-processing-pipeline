# Pipeline Generalisation

How the processing chain adapts to a **future mission** with different spectral
bands, a different calibration model, and a different native CRS — and which
parts of the design are sensor-specific versus generic.

The architecture already separates the two: the **core entities, reader
interface, windowing pipeline, QA and writer are sensor-agnostic**; the
sensor-specific knowledge is concentrated in (a) the data the reader discovers
and (b) the radiometric model. The sections below map each axis of change to the
seam that absorbs it.

---

## What is generic vs sensor-specific

| Element | Sensor-specific? | Where it lives | How it generalises |
|---|---|---|---|
| Band set, names, CWL, order, PAN-or-not | Data | STAC assets, filter file, metadata | Discovered at read time — no code change |
| Detector geometry (width, start rows, TDI, line period) | Data | `ImagerConfiguration` (metadata) | Parsed generically per acquisition |
| File layout & naming | Format | `PackageReader` | New layout → new `Reader` implementation |
| **Radiometric model** (the equation + coefficient semantics) | **Code** | `L1BCalibrator`, CPF parsing | Swap via `Calibrator` ABC + model registry |
| NoData, dtype, units | Param | Calibrator / writer args | Configuration values |
| Native CRS & sensor→ground model | Code (L1C) | (future) geometric stage | Target CRS as a parameter; LoS/ephemeris are generic |

The core entities (`Band`, `Acquisition`, `CalibrationParameters`, `Product`)
carry **no** mission proper nouns or sensor constants today.

---

## Axis 1 — different spectral bands

Already handled by data, not code. Bands are discovered from the STAC assets,
the band name↔id map from the filter file, and central wavelengths from the
metadata. Adding, removing or reordering bands (more red-edge bands, no PAN, a
SWIR band) requires **no code change**: the reader builds whatever band set the
package declares, and the pipeline iterates it.

**Caveat.** A band outside the VNIR range (e.g. SWIR/thermal) may need a
different radiometric or atmospheric treatment — that is an Axis-2 concern, not a
band-count one.

---

## Axis 2 — different calibration model

This is the **main code-touch point**. The current model is the dark/thermal +
non-uniformity + absolute chain, and `CalibrationParameters` models exactly the
three named coefficient arrays of this CPF. A future sensor might use a
polynomial gain, a per-module correction, a look-up table, or no thermal term.

**How to abstract.**
1. The `Calibrator` ABC already defines the seam (`calibrate(band, dn,
   line_start)`). Implement the new physics as another subclass (e.g.
   `PolynomialCalibrator`) — the pipeline, QA and writer are unchanged because
   they only see "a calibrator that turns a DN block into a value block".
2. Select the calibrator via a **model registry / factory** keyed by a model id
   declared in a sensor profile, instead of importing `L1BCalibrator` directly in
   the pipeline.
3. Generalise `CalibrationParameters` to hold a **named dict of coefficient
   arrays** plus scalars, rather than three fixed fields, so different CPF schemas
   parse into the same entity. The reader's CPF parsing becomes a small
   per-schema adapter.

After that, a new radiometric model = one new `Calibrator` subclass + one
parsing adapter; nothing else moves.

---

## Axis 3 — different native CRS

L1B is produced in **sensor coordinates with no CRS**, so the implemented level
is already CRS-agnostic — a different native CRS changes nothing at L1B.

CRS becomes relevant at **L1C**, where it is simply a parameter of the geometric
stage: the sensor→ground model (line-of-sight + ephemeris + attitude + DEM) is
generic, and the **target CRS/GSD is chosen at warp time**. The writer already
accepts a profile; for L1C it would additionally take a `crs` and `transform`.
No part of the radiometric chain depends on the CRS.

---

## Worked example

> *A future 12-band mission: a polynomial radiometric model, native UTM output.*

Changes required:
1. **Sensor profile** (config): band definitions / file patterns, the
   radiometric `model: polynomial` with its parameters, quantisation, NoData,
   units, default CRS.
2. **`PolynomialCalibrator(Calibrator)`** implementing the new equation, plus a
   CPF parsing adapter for its coefficient schema.
3. Point the reader at the new file patterns (or a new `Reader` if the container
   differs).

Unchanged: the windowing pipeline, QA validator, GeoTIFF writer, quicklook and
CLI. Choosing UTM is an L1C parameter, not an L1B concern.

---

## Recommended abstraction work (to make the above turnkey)

- **Sensor profile files** (YAML) read by the reader and pipeline, so a new
  mission is onboarded by adding data, not editing code.
- **Calibrator registry** mapping a model id → `Calibrator` implementation.
- **Generic coefficient container** in `CalibrationParameters` (named arrays).
- **Reader registry** keyed by package/container type.

These keep the single most sensor-specific concern — the radiometric model —
behind a stable interface, while the rest of the chain stays untouched across
missions.
