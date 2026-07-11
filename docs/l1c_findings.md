# L1C — What Measurement Overturned

The design of the geometric level is in [`l1c_spec.md`](l1c_spec.md). This document
is about how it got there, and it exists because the honest answer is: **five of its
load-bearing assumptions were wrong, and measurement found every one of them.**

None of the five announced itself. Each produced a plausible product, a finite
number, a model that converged. Four of them looked exactly like an *irreducible
platform error* — the kind of thing you shrug at, attribute to the hardware, and
route into a refinement stage. That is what makes them worth writing down: not the
errors themselves, but the fact that only a check *that could have failed* told them
apart from the truth.

---

## The pattern

Every one of the five followed the same shape:

1. An assumption was made about the input, reasonably, from domain knowledge.
2. The assumption was never stated as a *question*, so nothing tested it.
3. The resulting error appeared as a plausible residual — not a crash, not a NaN
   (Not a Number), not an obviously wrong footprint.
4. It was found by asking **where the residual pointed**, rather than accepting that
   it was irreducible.

The framework's defence against this is not care. It is that the undocumented
properties of the input are made into **explicit parameters** and **measured**
(`geometry.Convention`), that derived quantities are **validated against the data
that produced them** (`Attitude.rates_are_consistent`, `terrain.check_undulation`),
and that every phase has an **exit criterion that can fail**.

---

## 1. The ephemeris is not Earth-fixed

**Assumed.** The platform position and velocity carry a numeric frame tag. The
package never defines it. The original specification read it as ECEF
(Earth-Centred, Earth-Fixed), which is the conventional choice and which the earlier
outline had also assumed.

**Cost.** The sub-satellite point landed at longitude 5.4° while the delivered
footprint sits at 56.1° — **5,000 km**. Every candidate convention failed by about
the same amount, which is what exposed it: the harness refused to name a winner
rather than crowning the least-bad of a bad field.

**Found by.** Two independent measurements, neither of which required trusting the
other:

- Rotating the ephemeris by the Earth-rotation angle put the sub-satellite point at
  55.8°, tens of kilometres from the footprint centroid instead of thousands.
- The delivered speed is **7,673.0 m/s**. The circular *inertial* orbital speed at
  that radius is **7,672.1 m/s** — agreement to 0.01%. An Earth-fixed speed would
  differ by about 440 m/s.

**Changed.** The ephemeris frame became the eighth axis of the convention harness,
and the largest lever in the whole model.

---

## 2. The bands were already co-registered

**Assumed.** The bands sit on different detector rows, so they look at different
along-track angles and image a given ground point up to 0.78 s and 5.7 km apart. The
specification concluded that the delivered rasters were therefore offset by hundreds
of lines, and that "a naive cross-correlation would have to search ~1500 lines".

**Cost.** None yet — but the entire band-co-registration section was written around
a problem that did not exist, and one interface (`LineTiming`) was designed to be
shared across bands when it must be per-band.

**Found by.** Cross-correlating the delivered L1B bands against each other. The
residuals were **a few pixels, not hundreds of lines**. Inspecting the session
metadata then showed why: each band carries its **own** per-line timestamps, offset
from its neighbours' by exactly `Δrow × line_period` (1504 rows → 777,568 µs, to the
microsecond). The product assembler had already staggered them.

**Changed.** Line timing is per band. And the real problem was identified: what
remains is a ~0.5% scale error in that stagger (the assembler used an integer line
count, which assumes the *nominal* altitude-to-ground-speed ratio, and this
acquisition is not on the nominal orbit), plus yaw coupling and terrain parallax.

**Worth noting:** the wrong assumption was *more alarming* than the truth. It would
have led to building a large-search-range matcher for a problem that needed a 0.5%
correction.

---

## 3. Precession, mistaken for a platform bias

**Assumed.** The Earth-rotation model applied only Greenwich Mean Sidereal Time. Its
own docstring justified the omission: *"precession, nutation and polar motion are
neglected — they are arcsecond-level effects, far below the pointing uncertainty
that dominates this level's error budget."*

**Cost.** **32.8 km** of footprint error — and this is the dangerous one, because the
specification *predicted* a large platform bias. The acquisition has no GNSS (Global
Navigation Satellite System) lock, so its positions are propagated rather than
measured, and a 32.8 km error looked exactly like the consequence. It would have been
absorbed into the absolute-refinement stage as if it were real, and the refinement
would have "worked".

**Found by.** Asking which *direction* the residual pointed, instead of accepting it.
The offset was −33.9 km east and −16.5 km north. Precession in right ascension over
the elapsed 26.2 years predicts **−33.5 km**; in declination, **−16.2 km**. Both
matched to 2%.

The justification was wrong by four orders of magnitude. Precession accumulates at
about 50 arcseconds **per year**, so a quarter-century past J2000 it is a third of a
degree — tens of kilometres at orbital radius, not arcseconds.

**Changed.** The full IAU-76/FK5 chain — precession → nutation → apparent sidereal
time. Footprint error fell to **208 m**.

---

## 4. The angular rates are not the attitude's derivative

**Assumed.** The telemetry delivers angular rates alongside the attitude quaternions.
The specification argued for rate-aware interpolation on the grounds that the
derivatives are delivered, so using them is free accuracy.

**Cost.** The band-to-band error rose from **18 m to 123 m** — a "refinement" that
made the model **seven times worse**, and which presented as an irreducible residual.

**Found by.** Comparing the delivered rates against the rates implied by the
quaternion sequence itself. They are **~55× larger**. They describe the body's motion
in some other frame — most likely an inertial one, where an Earth-pointing platform's
rate is dominated by the orbital rate — and not the drift of the small attitude offset
the quaternions encode.

**The mechanism is not the obvious one**, which is why it is worth recording. A slope
that is merely the wrong *magnitude* largely cancels: the two Hermite endpoint weights
sum to `s(2s−1)(s−1)`, which vanishes at the midpoint of every segment. The damage
comes from a rate pointing along a **different axis**, which injects rotation the
quaternions never contained.

**Changed.** Not a flag — a **check**. `Attitude.rates_are_consistent` compares the
rates against the quaternion sequence, and interpolation falls back to SLERP
(Spherical Linear intERPolation) when they disagree. Asking for rate-awareness is a
request, not an instruction: the data gets a veto, because the caller cannot know what
a given payload contains. Position keeps Hermite — its velocity *does* agree with its
positions, and the same check passes.

---

## 5. The model alone cannot co-register the bands

**Assumed.** The implementation plan had Phase 3 (resampling through the physical
model) already improving band co-registration, with Phase 4 (self-calibrating the
missing line-of-sight term) polishing it from good to sub-pixel.

**Cost.** Nothing — because the exit criterion caught it. This is the one finding that
was *designed for*.

**Found by.** The Phase 3 exit criterion, written to be able to fail: *"the band
residual must fall below the 1–7 px it starts at. If it does not shrink, stop — the
model is wrong."* It did not shrink: **30.2 m before, 33.2 m after**.

The model carries a per-band error of its own, and it is the *same size* as the
misregistration it was meant to remove: it places the RE3 band wrong by a constant
**(−5.1, +6.5) px**.

**The discriminating question** was whether that residual is a missing *calibration*
or a modelling *bug* — because feeding a residual back into the model would paper over
either one indiscriminately. The discriminator is whether it **varies with position**:

| | along-track correlation | cross-track correlation |
|---|---:|---:|
| line residual | −0.09 | −0.05 |
| column residual | +0.01 | −0.17 |

It does not vary. The residual is constant across the whole strip and the whole swath
— the signature of a fixed angular offset per band, i.e. the **interior orientation**.
That is exactly the per-band line-of-sight term the Calibration Parameter File ships
unpopulated, and ~8 detector pixels of offset is an ordinary amount of optical
distortion for an instrument whose geometric calibration was never filled in.

**Changed.** The phasing. Phase 4 is not a refinement; it is the **prerequisite** for
co-registration existing at all. Phases 1–3 deliver a georeferenced, orthorectified
product with its absolute error declared — which is real and usable — but the bands do
not align until the calibration the instrument never shipped is estimated from the
imagery.

Had the exit criterion been written as *"the products are produced"*, Phase 3 would
have passed.

---

## What was measured, and what it cost

| # | Assumption | Error it produced | Looked like | Found by |
|---|---|---:|---|---|
| 1 | Ephemeris is Earth-fixed | 5,000 km | nothing worked | harness refusing to pick a winner |
| 2 | Bands offset by ~1500 lines | — | a harder problem | cross-correlating the delivered bands |
| 3 | Precession is negligible | 32.8 km | the predicted platform bias | asking where the residual pointed |
| 4 | Rates are the attitude's derivative | 105 m | an irreducible residual | checking them against the quaternions |
| 5 | The model can co-register the bands | — | a phase that would have passed | an exit criterion that could fail |

Four of the five were **indistinguishable from a real, irreducible platform error**
until they were measured against something independent.

---

## The three defences

**Make undocumented properties explicit parameters, and measure them.** The
telemetry's frames, quaternion conventions and detector axes are not documented. They
are collected in `geometry.Convention`, all 384 combinations are enumerated, and the
data chooses. Six of eight resolve decisively; the harness **reports the other two as
unresolved** rather than crowning the least-bad — see
[`l1c_spec.md`](l1c_spec.md) §4.3.

**Validate a derived quantity against the data that produced it.** The angular rates
are checked against the quaternion sequence. The geoid undulation is checked against
being identically zero — because PROJ, asked to convert orthometric to ellipsoidal
height without its grid, does not raise: it performs a "ballpark" transformation that
returns the height unchanged, and every elevation would be wrong by ~31 m with no
other symptom.

**Write exit criteria that can fail.** Phase 3's criterion failed, and that was the
most valuable single outcome of the phase. A criterion that only ever passes is not a
criterion.

---

## Still open

The findings above are settled. These are not, and they are recorded as assumptions
rather than results:

- **Two parities remain unresolved** — the scan direction (a north–south mirror) and
  the detector column sign (an east–west mirror). Both map the footprint's corners
  onto each other *and* displace every band identically, so both probes return
  bit-identical scores. They are not merely unmeasured; they are **provably invisible
  to geometry**, and only image content settles them.
- **Absolute accuracy is unvalidated.** The model reproduces the delivered footprint
  to 208 m — but that footprint was almost certainly derived by the provider from the
  *same* telemetry, so the agreement validates the implementation, not the orbit.
  Nothing in the package is independent of the ephemeris it would have to check.
- **The bands are not co-registered.** See finding 5.

See [`limitations.md`](limitations.md) for the product's honest state.
