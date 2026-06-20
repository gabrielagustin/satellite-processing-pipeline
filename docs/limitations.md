# Limitations, Failure Modes & Next Steps

## Failure modes

For each: what it is, how to **detect** it, and how to **handle** it in
production.

### 1. Temperature model assumes uniformly-spaced telemetry

The per-line darkfield uses a temperature profile built by spreading the detector
telemetry samples **uniformly** across the lines. The samples actually carry
their own timestamps (`ImagerTime`); if they are irregularly spaced or have gaps,
the uniform assumption mis-assigns temperature to lines, causing along-track
**banding** from over/under-subtracted dark signal.

- **Detect.** Compare the uniform mapping against the true line↔time mapping from
  `ImagerTime`; flag if the residual exceeds a fraction of a sample interval.
  Monitor along-track mean radiance for steps correlated with temperature
  transitions.
- **Handle.** Map samples to line times using `ImagerTime` and the per-line
  timing (exact interpolation); interpolate across gaps and flag affected lines;
  fall back to the scene-mean temperature only if telemetry is missing entirely,
  recording the fallback in the QA report.

### 2. Missing or mismatched calibration entry

The calibrator selects a radiometric entry by `(band, start_row, tdi)`. If the
acquisition runs in a configuration absent from the CPF (e.g. a different TDI or
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

- **Detect.** Set `saturation_dn` from the quantisation/ADC configuration (e.g.
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

- **Exact temperature timing.** Use `ImagerTime` to map telemetry to line times,
  removing the uniform-spacing assumption (failure mode #1).
- **Test suite + CI.** Unit tests per module (synthetic fixtures) and an
  integration test on a small window, wired into CI.
- **TOA reflectance.** A quick, high-value step toward L2: compute per-band ESUN
  from the provided solar + filter assets and emit TOA reflectance.
- **L1C prototype.** Implement the line-of-sight + ephemeris/attitude + DEM
  geolocation, with GCP refinement to compensate for the absent GNSS lock.
- **STAC output.** Emit a STAC item for the L1B product (metadata/catalogue
  layer) for interoperability.
- **Sensor profiles via config.** Externalise per-sensor parameters into config
  files so a new mission is onboarded without code changes.
- **Parallelism.** Bands are independent — process them concurrently
  (multiprocessing/Dask) to cut wall-clock time.
- **Detector artefacts.** Add dead/hot-pixel interpolation and de-striping.
