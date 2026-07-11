# L1C — Limitations

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

### Band co-registration is partly solved, and the rest is attitude

The missing interior orientation is now **estimated from the imagery**: the sensor model
inverts the measured band-to-band displacement into a per-band line-of-sight offset, and
for **five of seven bands** that drives the systematic error from 1.4–5.3 px to
**~0.005 px**. Those offsets are written out in the same schema the Calibration Parameter
File ships empty.

**Two bands are refused**, and correctly. The estimator declines to fit a band whose
residual varies with position, because a constant angular offset cannot represent one.
The two refused bands (B and RE3) are the furthest from the reference on the detector,
and the position-dependence grows with each band's **time separation** from the reference
(correlation +0.82). That is not optics — it is **attitude**: the bands are up to 0.47 s
apart, the attitude jitters by ~0.005° about a smooth fit (about 9 px on the ground), and
that jitter does not cancel across half a second.

**The attitude part is corrected in the product, not in the calibration.** The
along-track-varying residual is fitted as a low-degree polynomial and applied to the
band's geolocation grid (`spp.refine.scene`). It is recorded as a **scene** property and
is never written back to the instrument model — doing so would produce a per-band
"optical" constant that is really a snapshot of one pass's pointing noise, wrong for
every other acquisition while looking, on this one, like a triumph. Keeping the two
corrections in separate places is the only thing that distinguishes a co-registered image
from a corrupted instrument model.

**Measured in the delivered `stack.tif`** (bands against PAN, 4.0 m pixels):

| Band | Median | 90th percentile |
|---|---:|---:|
| G | 1.87 px | 3.15 px |
| R | 2.28 px | 3.42 px |
| B | 2.80 px | 4.72 px |

Against ~8 px in the delivered L1B, and ~8.3 px in L1C without any calibration. So the
bands **are co-registered to 2–3 px, and are not co-registered to the sub-pixel target**.
The floor is attitude jitter (~0.005° root-mean-square about a smooth fit, about 9 px on
the ground) plus matching scatter, and neither a constant offset nor a smooth polynomial
removes it.

**Next step for the last few pixels.** Estimate the attitude itself — a filtered or
smoothed attitude, or a per-line correction driven by dense matching — rather than
absorbing its effect into a per-band displacement field.

### Absolute accuracy is unvalidated

The model reproduces the delivered footprint to 208 m — but that footprint was almost
certainly derived by the provider from the **same telemetry** the model consumes, so the
agreement validates the implementation and says nothing about the orbit. Nothing in the
package is independent of the ephemeris it would have to check.

The acquisition also has **no satellite-navigation lock**: its positions are propagated,
not measured.

**Next step.** Match against an external reference orthoimage
([`l1c_spec.md`](spec.md) §8.1). Until then the product carries the
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

- A STAC item and a map-projected quicklook for the L1C product.
- The `qa_report_l1c.json` the pipeline emits covers the conventions, grid, terrain,
  interpolation, co-registration and flags, but is not yet merged into the single
  `qa_report.json` the L1B stage writes.
