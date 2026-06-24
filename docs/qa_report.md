# QA Report Reference

Every L1B run writes a `qa_report.json` to the output directory. This document
defines its schema field by field, the quality flags and their thresholds, and
shows a real example.

## What this QA is (and is not)

The report is a **radiometric quality assessment** of the calibrated radiance
product: per-band descriptive statistics plus threshold-based quality flags. It
measures **internal consistency and physical plausibility** — it is
**non-destructive** (it never modifies the product) and computed in a single
streaming pass.

It is **not** an absolute radiometric accuracy / cal-val assessment: it does not
compare against reference targets, vicarious calibration or another sensor. That
is out of scope (see [`remaining_levels.md`](remaining_levels.md)). A passing QA
means "the product is well-formed and plausible", not "the radiance is absolutely
accurate".

---

## Structure

```jsonc
{
  "scene_id": "<source scene id>",
  "level": "L1B",
  "passed": true,            // true only if every band passed
  "bands": {
    "<BAND>": { ...per-band metrics... },
    ...
  }
}
```

---

## Per-band fields

Statistics are computed over **valid** pixels only (finite radiance, i.e.
excluding NoData/NaN and any non-finite values). Radiance units are
`W / (m² · sr · µm)`.

| Field | Type | Definition |
|---|---|---|
| `band` | str | Band name. |
| `pixels` | int | Total pixels assessed (whole band). |
| `valid` | int | Finite radiance pixels (not NaN, not ±inf). |
| `nodata` | int | NoData pixels (NaN) — pixels masked from invalid input DN. |
| `nodata_fraction` | float | `nodata / pixels`. |
| `min` | float | Minimum valid radiance (NaN if no valid pixels). |
| `max` | float | Maximum valid radiance (NaN if no valid pixels). |
| `mean` | float | Mean valid radiance. |
| `std` | float | Standard deviation of valid radiance. |
| `negative` | int | Valid pixels with radiance `< 0` (near-dark pixels after dark subtraction). |
| `negative_fraction` | float | `negative / valid`. |
| `non_finite` | int | Count of `±inf` outputs (should be 0; non-zero signals a numerical problem). |
| `saturated` | int | Pixels with `DN ≥ saturation_dn`. **Only present** when `--saturation-dn` is set. |
| `flags` | list[str] | Quality flags raised (see below); empty means clean. |
| `passed` | bool | `true` when `flags` is empty. |

> **Why negatives are kept.** Radiance is not clipped: slightly-negative values at
> near-dark pixels are preserved so QA can *measure* dark-subtraction bias rather
> than hide it. See [`decision_log.md`](decision_log.md) §4.

---

## Quality flags

| Flag | Raised when | Controlled by |
|---|---|---|
| `all_invalid` | `valid == 0` (no finite pixels). | — |
| `zero_dynamic_range` | `valid > 0` and `max == min` (flat band). | — |
| `high_nodata` | `nodata_fraction > max_nodata_fraction`. | `max_nodata_fraction` (default 0.50) |
| `negative_radiance` | `negative_fraction > max_negative_fraction`. | `max_negative_fraction` (default 0.05) |
| `non_finite_values` | `non_finite > 0`. | — |
| `saturation` | saturation checked and `saturated > 0`. | `--saturation-dn` (default: off) |

A band `passed` only if it raises **no** flags; the product `passed` only if
**every** band passed. The `run_l1b` CLI exit code is `0` when the product passed
and `1` otherwise, so it can gate automated pipelines.

### Thresholds / configuration

Defined by `QAThresholds` in
[`../src/spp/qa/radiometric_validator.py`](../src/spp/qa/radiometric_validator.py):

| Parameter | Default | Meaning |
|---|---|---|
| `saturation_dn` | `None` (off) | DN at/above which a pixel counts as saturated. |
| `max_negative_fraction` | `0.05` | Negative-radiance fraction that raises `negative_radiance`. |
| `max_nodata_fraction` | `0.50` | NoData fraction that raises `high_nodata`. |

---

## Example (excerpt from a real 8-band run)

```jsonc
{
  "scene_id": "...",
  "level": "L1B",
  "passed": true,
  "bands": {
    "R": {
      "band": "R", "pixels": 126763008, "valid": 126763008,
      "nodata": 0, "nodata_fraction": 0.0,
      "min": 6.534, "max": 784.204, "mean": 39.411, "std": 23.247,
      "negative": 0, "negative_fraction": 0.0,
      "non_finite": 0, "flags": [], "passed": true
    },
    "RE3": {
      "band": "RE3", "pixels": 126763008, "valid": 126763008,
      "nodata": 0, "nodata_fraction": 0.0,
      "min": -5.833, "max": 467.729, "mean": 23.767, "std": 17.859,
      "negative": 76211, "negative_fraction": 0.0006,
      "non_finite": 0, "flags": [], "passed": true
    }
  }
}
```

Reading it: both bands are fully valid (no NoData), with plausible radiance
ranges. `RE3` has a small fraction of negative pixels (0.06 %) — well under the
5 % threshold, so no `negative_radiance` flag — consistent with near-dark pixels
after dark subtraction.
