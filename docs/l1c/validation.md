# L1C — Validation Strategy

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
cause to the interior orientation. See [`l1c_findings.md`](findings.md).

A criterion that only ever passes is not a criterion.

### What is not validated

- **Absolute geolocation**, against anything independent of the telemetry.
- **Band co-registration** — it does not yet work; see above.
- **The two parities** — carried as declared assumptions, flagged in the quality report.
