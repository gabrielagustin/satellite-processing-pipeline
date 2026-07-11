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
- **Test coverage.** A unit suite (`tests/`) runs on every push/PR via GitHub
  Actions (see [`validation_strategy.md`](validation_strategy.md)), but it is
  still narrow — it covers the L1B calibrator core (the DN→radiance path), the
  temperature-timing model and a couple of reader helpers, not yet every module.

---

## Reflection — what I would do next

- **Broaden test coverage.** CI (GitHub Actions) and a focused unit suite
  (calibrator core + temperature timing) are in place; extend per module (reader
  parsing, QA flags, writer round-trip) plus an integration test on a small
  window.
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
- **Structured logging.** The stages emit `logging` warnings (e.g. dropped
  temperature timing, skipped bands), but there is no centralised setup yet —
  add a configured root logger with a consistent format, wired to logging levels
  (extending the existing `--quiet` and adding a `--verbose`), so
  provenance-relevant events are captured rather than relying on the default
  last-resort handler.

---

## L1C — the geometric level

The level is **partly implemented**, and the parts are not equally finished. What
follows is the product's real state, stated plainly, because every gap below produces a
file that *looks* correct: it opens in any geographic information system, it overlays a
basemap, the bands stack without complaint.

| Capability | State |
|---|---|
| Georeferencing | **Works.** Reproduces the delivered footprint to 208 m |
| Orthorectification | **Works.** Terrain intersection converges everywhere; the geoid datum is handled and checked |
| Band co-registration | **Does not work.** The bands are still ~8 px apart |
| Absolute accuracy | **Unvalidated** against anything independent of the telemetry |

### Band co-registration does not yet work

Resampling through the physical model left the band residual where it started: 30.2 m
before, 33.2 m after. The model carries a per-band error of its own — a constant
(−5.1, +6.5) px — the same size as the misregistration it was meant to remove.

The cause is diagnosed, not merely suspected: the residual is **constant across the
whole strip and the whole swath** (correlations with along-track and cross-track
position: −0.09, −0.05, +0.01, −0.17). A residual that does not vary with position is a
fixed angular offset per band — the **interior orientation** — which is precisely the
per-band line-of-sight term the Calibration Parameter File ships unpopulated. Eight
detector pixels of offset is an ordinary amount of optical distortion for an instrument
whose geometric calibration was never filled in.

**Next step.** Estimate the missing calibration from the imagery (the self-calibration
of [`l1c_spec.md`](l1c_spec.md) §8.2). It is not a refinement; it is the prerequisite.

### Absolute accuracy is unvalidated

The model reproduces the delivered footprint to 208 m — but that footprint was almost
certainly derived by the provider from the **same telemetry** the model consumes, so the
agreement validates the implementation and says nothing about the orbit. Nothing in the
package is independent of the ephemeris it would have to check.

The acquisition also has **no satellite-navigation lock**: its positions are propagated,
not measured.

**Next step.** Match against an external reference orthoimage
([`l1c_spec.md`](l1c_spec.md) §8.1). Until then the product carries the
`absolute_accuracy_unvalidated` flag, and any absolute figure quoted from it would be
false precision.

### Two telemetry conventions are assumed, not resolved

The scan direction (a north–south mirror) and the detector column sign (an east–west
mirror) cannot be determined from geometry. Both map the footprint's corners onto each
other *and* displace every band identically, so every geometric probe returns
bit-identical scores. They are **provably invisible** to the model and need image
content to settle. They are carried as declared assumptions and flagged.

A product that is silently mirrored is much harder to catch than one that says it might
be.

### The warp is not memory-bounded

L1B streams in windows and its peak memory is independent of strip height. **The L1C
warper is not**: it holds the whole source band (~507 MB) and the whole destination
(~1.15 GB) in memory, so it uses ~1.7 GB per band and grows with the strip. Workable at
this scale, not at the next. See [`performance.md`](performance.md).

### Not yet built

- A `spp-l1c` command-line entry point (the level currently runs as a library).
- The geometric section of `qa_report.json` is **specified** ([`qa_report.md`](qa_report.md))
  but not yet emitted.
- A STAC item and a map-projected quicklook for the L1C product.
