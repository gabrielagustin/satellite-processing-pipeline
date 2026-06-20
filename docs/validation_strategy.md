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
`high_nodata`, `negative_radiance`, `saturation`, `zero_dynamic_range`.

**Physical plausibility (full scene).** Across all eight bands the run passes QA,
with radiance magnitudes in the expected range (e.g. red ≈ 40 W·m⁻²·sr⁻¹·nm⁻¹)
and the band-mean radiance **decreasing from the visible to the NIR**, as
expected for the scene. The RGB quicklook reproduces recognisable surface
features (coastline, water, land), an end-to-end visual sanity check.

**Reproducibility.** Per-band statistics and flags are written to
`qa_report.json` on every run for audit and regression comparison.

---

## Ongoing / recommended

- **Unit tests** per module with small synthetic fixtures (reader parsing,
  calibration math, QA flags, writer round-trip).
- **Integration test** on a small window, asserting against a stored reference.
- **Regression fixtures** — pin `qa_report.json` summary metrics for the
  reference scene and diff on change.
- **CI** to run the above on every commit.

---

## Remaining levels (see [`remaining_levels.md`](remaining_levels.md))

- **L1A** — round-trip against a reference L1A product; line count vs timing.
- **L1C** — GCP/tie-point RMSE against a reference orthoimage; band
  co-registration error; computed footprint vs the STAC `bbox`.
- **L2A** — comparison to a coincident reference surface-reflectance product over
  stable targets; AERONET-based checks; physical sanity (NIR over deep water ≈ 0).
