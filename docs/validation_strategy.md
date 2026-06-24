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

**L1B calibration core.** Unit tests (`tests/test_temperature_timing.py`) pin the
public `calibrate()` path: the documented DN→radiance formula (checked against an
independent reference and hand-computed pixels), the windowed-`line_start`
equivalence to the full frame that makes streaming exact, NoData masking, and the
shape/bounds guards.

**Temperature line-timing.** Unit tests (`tests/test_temperature_timing.py`)
cover the timestamp→line mapping, its uniform fallback, and the edge cases
(missing timing, single sample). The impact was quantified on the real scene by
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
