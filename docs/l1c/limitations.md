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
R 2.3 px, G 2.0 px, B 2.8 px median — against ~8 px in the delivered L1B.

**The residual is real, and it is not correctable by any displacement field.** Three
measurements say so, and they rule out the obvious explanations:

- It does **not shrink with window size** (128 → 1024 px leaves it at 2.0–2.9 px), so it
  is not matching noise. Averaging over 64× the area would have hidden noise; it does not
  move.
- A **dense two-dimensional displacement field does not beat a degree-3 polynomial** on
  held-out windows, so it is not a smooth spatial function.
- Fitting a **per-line profile and testing it on withheld columns makes it worse** for
  some bands (2.89 → 3.83 px). A dense fit is fitting noise.

A geometric error from attitude would be **spatially coherent** — the attitude is common
across the swath at any instant. This residual is not. What varies rapidly in both line
and column, and is stable in time, is **the scene itself**: phase correlation between two
*spectrally different* bands finds the shift that best aligns their gradients, and two
bands genuinely see different edges. A large part of the measured "misregistration" is
therefore **spectral, not geometric**, and no resampling corrects it because it is not a
position error.

**Two things that look like defects and are not:**

- **Ships appear displaced between bands.** They are. The bands image the same ground up
  to 0.78 s apart, and a vessel at 10 m/s moves 8 m — two pixels — in that time. The
  product is rendering real motion.
- **Water rendering as magenta in a red-green-blue composite** is a display stretch, not
  the data. With a per-band 2–98% stretch the water is the blue-cyan it should be.

**Open question, and the test that settles it.** How much of the 2–3 px is geometry and
how much is spectral? Measure the residual over spectrally *flat* terrain (bare desert)
against spectrally *varied* terrain (vegetation, urban). If it collapses over the desert
it is spectral, and there is nothing to correct.

### Absolute accuracy: corrected to ~20 m, against an external reference

The acquisition has **no satellite-navigation lock** — its positions are propagated, not
measured — and the resulting pointing bias was **919 m**, almost entirely along-track.
Nothing internal to the level could see it: the bands are co-registered onto *each other*,
and the delivered footprint was derived by the provider from the *same telemetry* the
model consumes. A bias shared by the whole product is invisible to both.

It is now corrected against a **reference orthoimage** (searched by catalogue, ranked by
closeness in time; for the reference acquisition, a same-day scene at 0.9% cloud), with
the correction applied to the **boresight** rather than as a shift of the output image.

Verified against two references that share no code path with the estimator:

| | |
|---|---:|
| uncorrected | 919 m |
| corrected — elevation-model coastline | **19 m** |
| corrected — Sentinel-2, same day | **29 m** |

**Two caveats that are not rounding errors.**

The reference's *own* absolute accuracy is specified at ~11 m. We are measuring against
it, so its error is inside ours. At ~20 m we are within a factor of two of the ruler's own
precision, and pressing further without a better reference would be false precision.

And the along-track correction is **degenerate with a clock offset** (+0.126 s here). From
a single strip the two are not separable. Attributing it to pitch is a **convention**, not
a measurement, and the quality report says so.

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
