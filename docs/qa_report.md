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

---

## The `geometric` section (L1C)

The geometric level extends `qa_report.json` with a `geometric` block. Its design
principle is the same one that governs the level: **report what was measured, and say
plainly what was not.** A metric that is absent because a check could not be run is not
the same as a metric that passed, and the schema keeps them distinguishable.

```jsonc
"geometric": {
  "conventions": {
    "resolved":   { "ephemeris_frame": "eci", "attitude_frame": "lvlh", ... },
    "assumed":    { "scan_direction": 1, "column_sign": 1 },   // NOT resolvable by geometry
    "band_coherence_m": 18.0
  },
  "footprint": {
    "corner_rmse_m": 208.0,        // vs the delivered catalogue geometry
    "centroid_offset_m": 130.0,
    "is_absolute_accuracy": false  // see below -- this is model-vs-model agreement
  },
  "grid": {
    "crs": "EPSG:32640",
    "gsd_m": 4.0,
    "native_gsd_m": 3.77,          // derived from the ephemeris, never a data sheet
    "shape": [29403, 9774],
    "data_coverage_pct": 40.0
  },
  "terrain": {
    "source": "Copernicus GLO-30",
    "geoid": "EGM2008",
    "undulation_range_m": [-32.3, -24.6],
    "converged_pct": 100.0,
    "void_fraction": 0.43          // open water, filled with sea level
  },
  "interpolation": { "lattice_step_px": 8, "max_error_px": 0.026 },
  "coregistration": {
    "reference_band": "PAN",
    "bands": { "RE3": { "residual_px": [-5.1, 6.5], "before_px": [3.9, -6.6] }, ... },
    "max_residual_px": 8.2,
    "status": "NOT_ACHIEVED"       // see the flag below
  },
  "los_correction": null           // populated by the self-calibration
}
```

### Flags

| Flag | Meaning |
|---|---|
| `no_gnss_lock` | The platform had no satellite fix; its positions are propagated, not measured. Always raised on the reference acquisition |
| `absolute_accuracy_unvalidated` | The footprint check compares against a geometry the provider derived from the *same* telemetry. It validates the model, **not the orbit** |
| `conventions_assumed` | One or more telemetry conventions could not be resolved from the data and are carried as assumptions |
| `coregistration_not_achieved` | The band residual did not fall below its starting value: the interior orientation is not calibrated |
| `telemetry_gap` | A line's exposure interval departs from the modal cadence — dropped or duplicated lines |
| `dem_voids` | Elevation-model voids over **land** (voids over water are expected and filled with sea level) |
| `intersection_not_converged` | Rays whose terrain intersection did not settle within the iteration cap |
| `zone_straddle` | The footprint crosses a UTM zone boundary; the centroid's zone was used |
| `partial_band_coverage` | Grid cells covered by some bands but not all — the strip ends, where the stagger runs off the raster |

### Two flags that exist to stop a lie

`absolute_accuracy_unvalidated` and `coregistration_not_achieved` both mark things that
**look fine in the product**. The rasters are georeferenced, they open in any geographic
information system, they overlay a basemap plausibly, and the bands stack without
complaint. Nothing about the file says that its absolute position has never been checked
against anything independent of the telemetry that produced it, or that its bands are
still eight pixels apart.

A quality report that only reported successes would be silent on both.
