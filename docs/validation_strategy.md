# Validation Strategy

How correctness is established for the implemented L1B level, and how the
remaining levels should be validated. The framework is **deterministic** (same
inputs + configuration → identical outputs), which makes every check below
repeatable.

---

## Layers

1. **Numerical correctness** — the algorithm matches an independent reference.
2. **Streaming correctness** — windowed processing equals whole-array processing.
3. **Physical plausibility** — outputs make radiometric/physical sense.
4. **Quality assessment** — per-band metrics and flags on every run.

---

## L1B — what was validated

**Numerical reference.** The calibrator output was compared, on a real window,
against an independent `float64` NumPy implementation of the calibration chain.
Maximum difference: **2.3 × 10⁻⁵** (floating-point rounding only) — confirming
the production path computes the intended physics.

**Window invariance (streaming).** Calibrating lines `[r0:r1]` gives an
**identical** result whether done as one block or split into sub-blocks. This
guarantees the windowed pipeline produces the same product as an in-memory one
and that `line_start`-dependent terms (the temperature profile) are indexed
correctly.

**NoData masking.** The number of NaN output pixels exactly equals the number of
input NoData (DN = 0) pixels; all non-NoData outputs are finite.

**QA accumulator.** Incremental accumulation across windows yields the same mean
and standard deviation as a single-shot computation, so the streamed QA is exact.

**QA flag logic.** Synthetic cases exercise each flag: `all_invalid`,
`high_nodata`, `negative_radiance`, `saturation`, `zero_dynamic_range`. The full
per-band metric and flag schema is documented in [`qa_report.md`](qa_report.md).

**L1B calibration core.** Unit tests (`tests/test_l1b_calibrator.py`) pin the
public `calibrate()` path: the documented DN→radiance formula (checked against an
independent reference and hand-computed pixels), the windowed-`line_start`
equivalence to the full frame that makes streaming exact, NoData masking, and the
shape/bounds guards.

**Temperature line-timing.** Unit tests cover the timestamp→line mapping
(`tests/test_package_reader.py`), its uniform fallback and the edge cases
(missing timing, single sample) (`tests/test_l1b_calibrator.py`). The impact was quantified on the real scene by
comparing, for the R band, the **per-line mean TOA spectral radiance** (the L1B
output quantity, `W·m⁻²·sr⁻¹·µm⁻¹`) under the timestamped vs uniform models: the
band-mean radiance barely moves (39.38 vs 39.41) but the **along-track per-line
mean** shifts by up to **~3.5 %**, confirming the timestamped model corrects a
banding gradient that scalar (whole-band) QA cannot see.

**Physical plausibility (full scene).** Across all eight bands the run passes QA,
with radiance magnitudes in the **physically expected range** and the band-mean
radiance **decreasing from the visible to the NIR**.

- *Magnitude.* Red ≈ 40 W·m⁻²·sr⁻¹·µm⁻¹ (≈ 0.04 W·m⁻²·sr⁻¹·nm⁻¹). Cross-check
  against the package solar spectrum: at 665 nm `E_sun ≈ 1.55 W·m⁻²·nm⁻¹`
  (`spectral_solar.json`), so for a surface reflectance ρ and high Sun,
  `L ≈ ρ·E_sun·cos θ_s / π ≈ 0.04 W·m⁻²·sr⁻¹·nm⁻¹` at ρ ≈ 0.1 — matching the
  product. (A red radiance of 40 *per nm* would exceed the solar irradiance and
  is unphysical; the unit is **per µm**.)
- *Spectral order.* The visible→NIR decrease follows from the solar spectral
  irradiance peaking in the visible and falling toward the NIR, combined with low
  NIR reflectance over the predominantly water/coastal surface.

The RGB quicklook reproduces recognisable surface features (coastline, water,
land), an end-to-end visual sanity check. See
[`references.md`](references.md) for the solar-spectrum references underpinning
the expectation.

**Reproducibility.** Per-band statistics and flags are written to
`qa_report.json` on every run for audit and regression comparison.

---

## Ongoing / recommended

- **Unit tests** — a focused suite covers the L1B calibrator core and the
  temperature-timing model; extend per module with small synthetic fixtures
  (reader parsing, QA flags, writer round-trip).
- **Integration test** on a small window, asserting against a stored reference.
- **Regression fixtures** — pin `qa_report.json` summary metrics for the
  reference scene and diff on change.
- **CI** — GitHub Actions runs `pytest -v` on every push/PR on Python 3.11
  (`.github/workflows/ci.yml`); extend it as the suite grows.

---

## Remaining levels (see [`remaining_levels.md`](remaining_levels.md))

- **L1A** — round-trip against a reference L1A product; line count vs timing.
- **L1C** — GCP (Ground Control Point)/tie-point RMSE against a reference orthoimage; band
  co-registration error; computed footprint vs the STAC `bbox`.
- **L2A** — comparison to a coincident reference surface-reflectance product over
  stable targets; AERONET-based checks; physical sanity (NIR over deep water ≈ 0).

---

## L1C — what was validated

The geometric level is validated against **things independent of the model**, which is
harder than it sounds: almost everything in the delivered package derives from the same
telemetry the model consumes, so agreeing with it proves less than it appears to. Where
an independent check exists it is used; where it does not, that is stated rather than
papered over.

### Numerical correctness

Every primitive is checked against a closed form or an analytic inverse, never against
the model's own output — a rotation library that is self-consistently wrong produces a
plausible product on the wrong piece of ground.

- **Quaternion algebra** — rotation matrices are orthonormal and right-handed;
  composition matches matrix composition; conjugation inverts; the exponential and
  logarithmic maps round-trip, including in the small-angle branch that production
  actually runs.
- **Geodetic conversion** — round-trips to millimetres, including at the poles and the
  equator, where a wrong formula fails.
- **Ray–ellipsoid intersection** — matches the analytic solution; returns the *near*
  side; and returns `NaN` for a ray pointing at the sky rather than a spurious hit
  behind the platform (a bug the tests caught).
- **Interpolation** — reproduces the sampled states *exactly* at the sample times.
  Hermite is checked against an exact circular orbit and must beat linear interpolation
  by a factor of 100, or the extra machinery earns nothing.
- **Relief displacement** — matches `height × tan(incidence angle)`, with the
  incidence angle at the ground, not the view angle at the platform (they differ by the
  Earth-curvature factor `(R + altitude)/R` = 1.063 at 400 km).

### Correctness against ground truth (synthetic)

The convention harness and the sensor model are tested on a **synthetic acquisition
whose geometry is known by construction** — a circular orbit, a yaw-steered platform, and
bands staggered in time exactly as a real product staggers them. This is what makes the
tests test the *harness* rather than the scene.

- The harness **recovers the five resolvable conventions** from the synthetic data.
- It **admits it cannot see** the two parities — reversing the scan mirrors the strip
  north–south and flipping the column sign mirrors it east–west, and both leave the
  footprint's corners and the band coherence bit-identical. A harness that claimed
  otherwise would be manufacturing a result, so the test pins the *admission*.
- Reading the ephemeris as Earth-fixed when it is inertial must fail **grossly**, not
  subtly — a close call is what lets an error like that slip through.

### Correctness on the real acquisition

| Check | Result |
|---|---|
| Footprint vs the delivered catalogue geometry | **208 m** corner RMSE |
| Strip dimensions vs the delivered footprint | 117.02 × 15.48 km against 117.03 × 15.44 km (0.03%) |
| Timing model vs the delivered acquisition window | reproduces 08:15:02.257 → 08:15:18.257 UTC to sub-second |
| Terrain intersection | 100% converged, within 2 iterations |
| Geoid undulation | −32.3 to −24.6 m, and **checked against being the silent zero** |
| Band-pair stereo coefficient | 0.01524 measured against 0.0152 predicted |
| Geolocation-lattice interpolation error | **0.026 px** (budget: 0.1) |

**The 208 m is not a measurement of absolute accuracy**, and treating it as one would be
the single easiest mistake to make here. The delivered footprint was almost certainly
derived by the provider from the *same* telemetry, so the agreement demonstrates that an
independent implementation reproduces theirs. It validates the **model**, not the orbit.
Genuine absolute accuracy requires matching against an external reference image, and it
is **not yet done**.

### The check that failed

Band co-registration was given an exit criterion phrased as a number it had to beat, not
an artefact it had to produce: *the residual must fall below the 1–7 px it starts at.*
**It did not** — 30.2 m before, 33.2 m after — and the failure is the most useful result
the level produced. The follow-up test (is the residual constant across the strip, or
does it vary?) distinguished a missing calibration from a modelling bug, and pinned the
cause to the interior orientation. See [`l1c_findings.md`](l1c_findings.md).

A criterion that only ever passes is not a criterion.

### What is not validated

- **Absolute geolocation**, against anything independent of the telemetry.
- **Band co-registration** — it does not yet work; see above.
- **The two parities** — carried as declared assumptions, flagged in the quality report.
