# L1C — Quality Report

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
