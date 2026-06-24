# Limitations, Failure Modes & Next Steps

## Failure modes

For each: what it is, how to **detect** it, and how to **handle** it in
production.

### 1. Temperature sample → line mapping (timestamped, with uniform fallback)

The per-line darkfield needs each detector-temperature sample placed at its
correct along-track line. The framework now does this **from the sample
timestamps** (`ImagerTime`), anchoring the line clock to the imager clock via the
`TimeSync` block and the acquisition start time — exact even when the telemetry
is irregularly spaced or **brackets** the imaging window (as it does in the
reference scene: telemetry spans ~19 s around a 16 s acquisition). The
provenance records `temperature_model: per_line_timestamped_interpolation`.

The residual risk is the **fallback**: when the timing inputs are absent (no
`TimeSync` anchor, no start time, or samples without `ImagerTime`), the code
spreads the samples **uniformly** across the lines
(`per_line_uniform_interpolation`). That assumes equal spacing *and* that the
telemetry exactly spans the acquisition — on the reference scene this uniform
model mis-assigns temperature by up to ~15 °C in places, biasing per-line
radiance by up to **~3.5 %** along-track (a banding gradient), while leaving the
band **mean** almost unchanged — which is why scalar QA does not catch it.

- **Detect.** Check that the timestamped model was used (provenance flag); if on
  the fallback, flag it. Monitor along-track mean radiance for a slow gradient or
  steps correlated with temperature transitions.
- **Handle.** Prefer the timestamped model (now the default). On the fallback,
  raise a QA flag and record it in provenance; interpolate across telemetry gaps
  and flag affected lines; fall back to the scene-mean temperature only if
  telemetry is missing entirely.

### 2. Missing or mismatched calibration entry

The calibrator selects a radiometric entry by `(band, start_row, tdi)`. If the
acquisition runs in a configuration absent from the CPF — Calibration Parameter File — (e.g. a different TDI (Time-Delay Integration) or
start row), the lookup raises and the run fails; a wrong-but-present entry would
silently miscalibrate.

- **Detect.** At read time, resolve every band's entry and assert it exists and
  that its coefficient-array length equals the detector width (the calibrator
  already guards the width). Cross-check `band`/`row`/`tdi` against the imager
  configuration.
- **Handle.** Fail fast with an explicit, actionable error listing the missing
  tuple; optionally fall back to the nearest available configuration with a loud
  warning and a QA flag, never silently.

### 3. Saturation passes through unflagged

By default the saturation check is disabled (`saturation_dn = None`) because the
saturation level depends on the quantisation mode. Saturated pixels then become
plausible-looking high radiance.

- **Detect.** Set `saturation_dn` from the quantisation/ADC (Analogue-to-Digital Converter) configuration (e.g.
  the 12-bit maximum) and inspect the saturated-pixel and high-radiance fractions
  in the QA report.
- **Handle.** Enable the threshold per acquisition mode; raise the `saturation`
  QA flag and, downstream, mask or down-weight saturated pixels.

### 4. NaN NoData propagation downstream

NaN is the natural float NoData, but tools that do not treat NaN as NoData (or
naive statistics) can propagate it or skew results.

- **Detect.** Validate that NaN locations match the input NoData mask; check
  downstream stats are computed NaN-aware.
- **Handle.** Document the convention; offer a sentinel-NoData output mode for
  consumers that need it; always write the NoData tag in the raster profile
  (done).

---

## Other limitations (scope)

- **No geometric correction.** Output is in sensor coordinates; no CRS. L1C is
  specified, not implemented.
- **No reflectance / atmospheric correction.** L2A is specified, not implemented.
- **Single-acquisition oriented.** No batch/catalogue orchestration yet.
- **Test coverage.** Validation so far is via end-to-end and numerical-reference
  checks (see [`validation_strategy.md`](validation_strategy.md)); a formal unit
  test suite + CI is not yet in place.

---

## Reflection — what I would do next

- **Test suite + CI.** A focused unit suite exists for the temperature-timing
  model (`tests/`); extend it per module (reader parsing, calibration math, QA
  flags, writer round-trip) plus an integration test on a small window, wired
  into CI.
- **L1C prototype.** Implement the line-of-sight + ephemeris/attitude + DEM
  geolocation, with GCP (Ground Control Point) refinement to compensate for the absent GNSS (Global Navigation Satellite System) lock.
- **TOA reflectance.** A quick, high-value step toward L2: compute per-band ESUN
  (exo-atmospheric solar irradiance) from the provided solar + filter assets and
  emit TOA reflectance.
- **Sensor profiles via config.** Externalise per-sensor parameters into config
  files so a new mission is onboarded without code changes.
- **Parallelism.** Bands are independent — process them concurrently
  (multiprocessing/Dask) to cut wall-clock time.
- **Detector artefacts.** Add dead/hot-pixel interpolation and de-striping.
