# L1C — Geometric Processing Specification

Georeferencing, orthorectification and band co-registration: the transition from
L1B top-of-atmosphere (TOA) spectral radiance in **sensor coordinates** to
radiance on a **projected map grid**, with all spectral bands resampled onto a
single common grid.

This document supersedes the L1C sketch in
[`remaining_levels.md`](remaining_levels.md) with an implementable design. It is
written for a generic multispectral **pushbroom (linescan)** imager; every figure
quoted as "measured" was derived from the reference acquisition used to develop
the framework, and the procedure that derived it is reproducible from the
delivered package alone.

Terminology follows the rest of the docs: **pipeline** = the sequence of product
levels (L0 → L1A → L1B → L1C → L2A); L1C is the **geometric** level.

---

## 1. Goal and scope

**In scope.**

1. **Georeferencing** — assign a ground coordinate to every detector sample via a
   rigorous physical sensor model (ephemeris + attitude + camera geometry).
2. **Orthorectification** — intersect each viewing ray with a terrain surface
   (a digital elevation model, DEM), removing relief displacement.
3. **Band co-registration** — resample all bands onto one common grid so that a
   given output pixel refers to the same ground location in every band.
4. **Geometric self-calibration** — estimate the pointing corrections the
   delivered calibration does not provide (see §3.3), from the imagery itself.

**Out of scope.** Pan-sharpening, mosaicking across acquisitions, and any
radiometric change: L1C resamples radiance, it does not alter its physics. The
resampling kernel is the only radiometric operator applied, and it is chosen to
be radiometry-preserving (§7.3).

**Output.** Per-band Cloud-Optimized GeoTIFF (COG) radiance rasters on a common
projected grid, plus a SpatioTemporal Asset Catalog (STAC) item, a geometric
quality-assessment (QA) report, and a map-projected quicklook.

---

## 2. What the input actually provides

Everything below was verified against the delivered package before writing this
spec. The verification matters: two of these findings overturn assumptions in the
earlier sketch.

### 2.1 Inventory

| Input | Content | Verified state |
|---|---|---|
| Per-line timing | `ExposureTimestamp` per raster line | **All 30 948 lines present**, tick = 1 µs, uniform step of exactly 517.0 ticks — no dropped lines |
| Time anchor | `TimeSync` pairs (imager clock ↔ platform clock) | Present; maps imager ticks to Unix time |
| Clock correction | `time_sync_offset` = −6.357 ms | Present in the STAC item |
| Ephemeris | Position + velocity, Earth-Centred Earth-Fixed (ECEF) | 229 + 229 samples @ ~4.1 Hz, none interpolated |
| Attitude | Quaternion + angular rates | 223 + 223 samples @ ~4.3 Hz, none interpolated |
| Telemetry coverage | vs the 16.00 s image window | **Brackets it** with 16–20 s of margin on both sides |
| Camera intrinsics | focal length 580 mm, pixel pitch 5.5 µm, principal point (2048, 1536), detector 4096 × 3072, detector→body reference frame | Present |
| Geometric calibration | boresight alignment (3×3), per-band line-of-sight (LoS) coefficients | **Placeholder — see §3.3** |
| Band detector rows | Per-band detector start row, 628 … 2132 | Present; this is what separates the bands |
| GNSS lock | `gnss_lock` | **false** — see §3.2 |
| Footprint | STAC `bbox` / `geometry` | Present; used as an independent check |

### 2.2 The timing model closes

Reconstructing each line's acquisition time from the per-line imager timestamps,
the `TimeSync` anchor and the clock offset yields an image window of
**08:15:02.257 → 08:15:18.257 UTC**, against a delivered STAC window of
**08:15:02 → 08:15:18**. The timing model therefore reproduces the delivered
metadata to sub-second agreement, and needs no calibration of its own.

Consequence: **use the per-line timestamps, not a nominal line period.** They are
complete, so the usual pushbroom hazard — dropped lines silently compressing the
along-track scale — does not arise here, but the code must still detect it
(§11, `telemetry_gap`).

### 2.3 The ephemeris is trustworthy in *shape*, and it is not the nominal orbit

Deriving the acquisition geometry from the ephemeris and cross-checking it
against the delivered footprint:

| Quantity | From ephemeris + camera | From the STAC footprint | Agreement |
|---|---|---|---|
| Ground sample distance (GSD) at nadir | 3.77 m | 3.77 m (across-track) / 3.76 m (along-track) | **~0.5 %** |
| Swath width | 15.52 km | 15.45 km | ~0.5 % |
| Strip length | — | 116.4 km | — |
| Platform altitude | 397.9 km | — | — |

Two things follow.

**The nominal sensor specification does not describe this acquisition.** The
data sheet quotes ~4.75 m GSD and a 19 km swath, which correspond to a ~501 km
reference orbit. This acquisition was taken from **~398 km**, giving **3.77 m**
GSD and a **15.5 km** swath. Any code or documentation that hard-codes the
nominal figures will be wrong by 25 %. The framework must derive GSD and swath
from the ephemeris, never from a constant.

**The ephemeris and the delivered footprint agree to ~0.5 %.** Despite
`gnss_lock = false`, the platform state is geometrically self-consistent with the
footprint the provider shipped. This de-risks the whole level: a rigorous sensor
model built on this ephemeris should land close, and the footprint becomes a
usable end-to-end check (§12). It says nothing about *absolute* accuracy — see
§3.2.

### 2.4 The bands are separated in time and in view angle

Each band is read from a different detector row, so each looks along a different
along-track angle and images a given ground point at a different instant. From
the delivered start rows and intrinsics:

| Band | Detector row | Along-track view angle | Ground offset | Time offset | Line offset |
|---|---:|---:|---:|---:|---:|
| RE3 | 628 | −0.493° | −3.43 km | −0.47 s | −917 |
| RE2 | 864 | −0.365° | −2.54 km | −0.35 s | −679 |
| RE1 | 1112 | −0.230° | −1.60 km | −0.22 s | −428 |
| NIR | 1312 | −0.122° | −0.85 km | −0.12 s | −226 |
| PAN | 1532 | −0.002° | −0.02 km | −0.00 s | −4 |
| R | 1724 | +0.102° | +0.71 km | +0.10 s | +190 |
| G | 1916 | +0.207° | +1.43 km | +0.20 s | +384 |
| B | 2132 | +0.324° | +2.25 km | +0.31 s | +602 |

The extremes (RE3 and B) are **5.7 km and 0.78 s apart** on the ground. This is
the entire band-co-registration problem, and it is why co-registration is a
*geometric* operation and not an image-alignment afterthought: a naive
cross-correlation would have to search ~1500 lines.

It also drives the DEM requirement, in a way that is easy to get backwards —
see §3.4.

---

## 3. The four things that make this level hard

### 3.1 The attitude reference frame is not documented

The attitude quaternions carry a numeric `frame` tag whose meaning is not
defined in the package. The measured behaviour settles it:

- The quaternion is **near-identity** (scalar part ≈ −0.9996) and drifts by only
  **0.078° over the full 56 s** of telemetry.
- A nadir-pointing platform in a ~398 km orbit rotates ~3.4° in 56 s relative to
  an inertial frame. If the frame were Earth-Centred Inertial (ECI) or ECEF, the
  quaternion would sweep by that much. It does not, by a factor of ~44.
- The vector part is dominated by a single component of ≈ −0.0275, i.e. a
  rotation of **≈ 3.15°** about one axis — the expected magnitude of the **yaw
  steering** a sun-synchronous platform applies to compensate Earth rotation.

**Conclusion: `frame 1` is a local orbital (LVLH — Local Vertical, Local
Horizontal) frame, and the quaternion is the small body-pointing offset from
nadir.** The model must therefore *construct* the LVLH frame from the ECEF
position and velocity, and compose the quaternion onto it.

This is a hypothesis backed by two independent numerical arguments, not a
certainty. It is resolved for good by the disambiguation harness in §4.

### 3.2 There is no GNSS lock

`gnss_lock = false`: the platform positions are propagated, not measured.
Internally consistent (§2.3) does not mean absolutely accurate — a propagated
orbit can be self-consistent while being displaced bodily by hundreds of metres
to kilometres.

Combined with the error budget (§9), where attitude knowledge dominates by an
order of magnitude, this means **the physical model alone cannot deliver a
defensible absolute accuracy figure.** Hence the image-based refinement in §8.1,
and hence the decision *not* to quietly snap the product to the delivered
footprint: doing so would make the footprint check circular and hide the real
error.

### 3.3 The geometric calibration shipped is a placeholder

The calibration parameter file's geometric section contains:

- `boresightAlignment` = the 3×3 **identity** matrix.
- `lineOfSight` = 6 coefficients per band, per axis — **identical for all eight
  bands**, and equal to `[1, 0, 0, 0, 0, 0]`.

Identical LoS coefficients across bands cannot be physically true: the bands sit
on different detector rows and demonstrably look in different directions (§2.4).
These are unpopulated defaults, not a calibration.

**Two consequences, and they are the crux of the design.**

1. The interior orientation must be built from **first principles** — the pinhole
   model implied by focal length, pixel pitch, principal point, per-band detector
   row, and the detector→body reference-frame matrix (§5.3). The calibration's
   LoS then enters as a *correction layer* on top, which happens to be identity
   here.
2. The framework should **estimate** the correction it was not given. The
   band-to-band residuals measured after physical orthorectification are exactly
   the per-band LoS offsets the file leaves empty. Feeding them back (§8.2) turns
   the missing calibration from a blocker into an output: the pipeline
   self-calibrates its interior orientation from the imagery.

### 3.4 The DEM matters for co-registration, not for absolute accuracy

The instinct is that a DEM matters because terrain displaces pixels. For this
long-focal-length, narrow-field instrument, that effect is *small*: the
cross-track half-field is only **1.11°**, so a 10 m height error displaces a
pixel by 10 · tan(1.11°) ≈ **0.19 m** — a twentieth of a pixel. Absolute
orthorectification is barely DEM-sensitive.

But the **bands are separated by 0.82° of along-track view angle** (§2.4) — 37 %
of the entire cross-track field of view. Two bands viewing the same terrain from
0.82° apart are a **stereo pair with a short baseline**, and terrain at height
`h` displaces them *relative to each other* by ≈ `h · 0.0142`:

| Terrain height | Relative band displacement | In pixels (3.77 m) |
|---:|---:|---:|
| 200 m | 2.8 m | 0.8 px |
| 500 m | 7.1 m | 1.9 px |
| 1500 m | 21.3 m | 5.6 px |

So the DEM is what holds the bands together over relief. Skipping it would leave
a **terrain-correlated, several-pixel band misregistration** that no global shift
can remove — and which would be invisible over water and glaring over hills.
This is the argument for orthorectifying rather than merely georeferencing.

---

## 4. Resolving the undocumented conventions

Four conventions in the input are ambiguous. Guessing wrong is not subtle — the
footprint lands hundreds or thousands of kilometres away — so they are resolved
empirically, once, by a cheap harness.

| Ambiguity | Candidates |
|---|---|
| Attitude reference frame | LVLH (§3.1, expected) · ECI · ECEF |
| Quaternion direction | body→reference · reference→body |
| Quaternion component order | scalar-first `(w, q0, q1, q2)` (expected) · scalar-last |
| Scan direction | raster row 0 = first line acquired · = last |

**Harness.** Forward-project only the **four image corners** under each
combination (a few dozen rays in total, milliseconds), and score each by the
distance between the computed footprint and the delivered STAC footprint. The
correct combination matches to within the model's error (kilometres at worst);
every incorrect one is off by tens to thousands of kilometres. The margin is
enormous, so the test is decisive rather than a fit.

The scan-direction test is the sharpest: reversing the line order flips the
footprint north↔south, a ~116 km error.

**Deliverable.** `tests/test_frame_conventions.py` pins the resolved combination,
and the resolution is recorded in [`decision_log.md`](decision_log.md). If a
future acquisition contradicts it, the test fails loudly rather than producing a
plausible, wrong product.

---

## 5. The forward sensor model

The core object: given a band and a detector sample `(line ℓ, column c)`, return
the ground point `(longitude, latitude, height)` it observed.

### 5.1 Line timing

For raster line `ℓ` of band `b`:

```
t(ℓ) = unix(TimeSync anchor) + (ExposureTimestamp[ℓ] − anchor_imager_ticks) · 1e-6
       + time_sync_offset
```

Per-line timestamps are authoritative; the nominal line period is used only to
*detect* gaps (a step deviating from the modal step flags `telemetry_gap`), never
to generate times. The band's detector row is **not** a time offset to be applied
here — it is a view-angle offset (§5.3), and the time difference between bands
emerges from the geometry rather than being imposed on it.

### 5.2 State interpolation

Interpolate platform state to `t(ℓ)`:

- **Position** — cubic **Hermite** interpolation using the sampled positions
  *and velocities*. Velocity is delivered, so it should be used: Hermite is exact
  to well under a millimetre here, against ~6 cm for linear interpolation. It is
  free, so there is no reason to accept the error.
- **Attitude** — the sampled quaternions *and* angular rates are both delivered,
  so interpolate with rate-aware cubic interpolation (Hermite in rotation-vector
  space, or SQUAD — Spherical and Quadrangle interpolation). Plain SLERP
  (Spherical Linear intERPolation) ignores the rates and is the fallback.

At ~4 Hz sampling with a body rate of ~7 mrad/s, the attitude changes ~0.094°
between samples — roughly **650 m on the ground**. Interpolation quality is
therefore not a detail; it is a first-order term in the error budget (§9), and it
is the reason for using the angular rates rather than discarding them.

### 5.3 The viewing ray

For band `b`, detector column `c`, with focal length `f`, pixel pitch `p`, and
principal point `(c₀, r₀)`:

1. **Detector → camera.** The band's samples lie on detector row
   `r_b` (its start row; with Time-Delay Integration the effective row is the
   centre of the accumulated window). The pinhole direction is

   ```
   d_cam = normalize( ( (c − c₀)·p ,  (r_b − r₀)·p ,  f ) )
   ```

   The `(c − c₀)` term sweeps the cross-track field; the `(r_b − r₀)` term is
   **constant per band** and is precisely the along-track view angle of Table
   §2.4. This single term is the mechanism that both separates the bands and, via
   the model, re-registers them.

2. **Calibration correction.** Apply the per-band LoS correction — a polynomial
   in normalized detector column giving `(Δalong, Δacross)` in radians. Identity
   as delivered (§3.3); populated by self-calibration (§8.2). This is the seam
   where a real geometric calibration would drop in unchanged.

3. **Camera → body.** Apply the detector→body reference-frame matrix from the
   intrinsics, then the boresight alignment matrix.

4. **Body → LVLH → ECEF.** Apply the interpolated attitude quaternion (body
   offset from nadir), then the LVLH→ECEF rotation built from the interpolated
   ECEF position and velocity:

   ```
   ẑ = −normalize(r)                     (nadir)
   ŷ = −normalize(r × v)                 (negative orbit normal)
   x̂ = ŷ × ẑ                             (completes the triad, ~along-track)
   ```

   The exact axis convention is one of the items pinned by §4.

The result is a unit ray direction in ECEF, originating at the interpolated
platform position.

### 5.4 Ground intersection

Iterative ray–terrain intersection:

1. Intersect the ray with the WGS84 (World Geodetic System 1984) ellipsoid →
   first guess `P₀`.
2. Sample the DEM at `P₀` → orthometric height `H`. Convert to an **ellipsoidal**
   height via the geoid undulation: `h = H + N(φ, λ)`.
3. Re-intersect the ray with the surface of constant ellipsoidal height `h` →
   `P₁`.
4. Repeat until the ground displacement between iterations is below a tolerance
   (default 0.1 px ≈ 0.4 m).

**The height datum is a real trap.** The DEM is referenced to the EGM2008
(Earth Gravitational Model 2008) geoid; the ephemeris is referenced to the WGS84
ellipsoid. In many regions the undulation is tens of metres. Mixing them is a
silent, systematic height error. It is small in *ground* displacement here (§3.4)
but it biases band-to-band registration over terrain, so it is applied properly
rather than neglected.

Convergence is fast — the narrow field means the ellipsoid guess is already
within a metre or two of the answer. Two iterations will typically suffice;
the loop is capped and non-convergence is flagged.

---

## 6. External data

| Data | Source | Purpose | Notes |
|---|---|---|---|
| DEM, 30 m | Copernicus DEM GLO-30 (COGs on AWS Open Data, no authentication) | Terrain intersection | Global; EGM2008-referenced; has voids over water — treated as height 0 (see below) |
| Geoid | EGM2008 grid, via PROJ's data package (`projsync`) | Orthometric → ellipsoidal height | Cached locally |
| Reference orthoimage | Sentinel-2 L2A (AWS Open Data / Copernicus Data Space / Planetary Computer) | Absolute refinement (§8.1) | Match by wavelength: red 665 nm ↔ S2 B4; near-infrared 842 nm ↔ S2 B8. Prefer a low-cloud, same-season scene |

**Water.** DEM voids over water are filled with the geoid (orthometric height 0),
which is the physically correct sea surface — not the ellipsoid. Voids over land
are interpolated and flagged (`dem_voids`).

**Offline operation.** All three sources are cached to a local directory and the
pipeline runs offline once populated. A `--dem` / `--reference` override accepts
any local raster, so no run depends on network availability at execution time.

---

## 7. Resampling

### 7.1 Mechanism: geolocation arrays, not GCPs

The forward model maps sensor → ground. Resampling needs the inverse (for each
output pixel, which detector sample?), and for a pushbroom sensor the inverse has
no closed form — each output pixel requires solving for the line time that saw
it. Doing that per output pixel is ~10⁹ solves across eight bands.

The standard solution, and the one specified here: evaluate the forward model on
a **subsampled grid**, write it as a **GDAL geolocation array**, and let GDAL
perform the inverse warp. This is the same mechanism used to grid swath products
such as those from MODIS and Sentinel-3.

- **Node spacing: 8 × 8 detector samples** (configurable). For the reference
  acquisition that is 512 × 3 869 nodes ≈ **32 MB** of `float64` longitude and
  latitude per band — trivial to hold in memory.
- **8 px ≈ 30 m on the ground, matching the DEM resolution.** This is the reason
  for 8 and not 64. The geolocation field has a *smooth* component (orbit and
  attitude, which a coarse grid would capture perfectly) plus a **terrain
  component that varies at DEM scale**. Subsampling at 64 px would low-pass the
  terrain correction and quietly undo part of the orthorectification. Node
  spacing must track the DEM, not the smoothness of the orbit.
- `float64` is mandatory: `float32` resolves longitude to only ~0.7 m at this
  latitude, which would inject noise comparable to the effects being corrected.

GCP (Ground Control Point) grids with a polynomial or thin-plate-spline fit were
considered and rejected: a low-order polynomial cannot represent
terrain-dependent displacement at all, and a spline over ~10⁴ GCPs is
prohibitively slow.

**Validation.** On one test patch, evaluate the model at full resolution and
compare against the interpolated geolocation grid. Require max error < 0.1 px,
and report it (`interp.max_error_px`).

### 7.2 Target grid

- **CRS (Coordinate Reference System)** — Universal Transverse Mercator (UTM)
  zone derived automatically from the footprint centroid (EPSG:326xx north /
  327xx south), with a polar-stereographic fallback beyond ±80°. Never hard-coded;
  overridable with `--crs`.
  - A 116 km strip can straddle a UTM zone boundary. The centroid's zone is used
    and the `zone_straddle` flag is raised when the footprint crosses one; the
    distortion is acceptable, and `--crs` is the escape hatch.
- **GSD** — default 4.0 m: the native 3.77 m rounded up to a clean grid. Rounding
  *up* slightly oversamples rather than decimating, avoiding aliasing.
  Configurable via `--gsd`; the native value is derived from the ephemeris and
  reported, never assumed.
- **Extent** — the **union** of all eight band footprints, snapped to GSD
  multiples. Because the bands are offset by up to 5.7 km along-track (§2.4), the
  strip ends have partial band coverage. The union preserves all data; the
  all-band intersection is computed and reported so a consumer can crop to fully
  co-registered coverage.

All bands warp to this **one** grid. That is what co-registers them.

### 7.3 Kernel and NoData

- **Bilinear** by default; cubic available. **Lanczos is deliberately not the
  default**: its negative lobes ring around saturated pixels and NoData edges,
  which corrupts a physical radiance quantity.
- Radiance is `float32` with **NaN** NoData, consistent with the L1B products.
  NaN must not be allowed to bleed through interpolation — invalid inputs are
  masked before resampling, not averaged into their neighbours.
- Per-band validity masks are written alongside the radiance.

---

## 8. Refinement

Both refinements share one component — a phase-correlation matcher — and both
apply their correction **inside the physical model**, then re-run geolocation.
This is a deliberate choice: shifting the output image instead would correct the
average error while leaving the terrain-dependent part wrong, and would leave the
sensor model knowingly incorrect. Corrections belong where the error physically
originates.

### 8.1 Absolute geolocation, against a reference orthoimage

With no GNSS lock (§3.2), absolute accuracy must be measured against external
truth, not asserted.

1. Produce a preliminary model-only L1C for one band with a good reference match
   (red or near-infrared).
2. Select match windows: on **land** (near-infrared / water-index threshold),
   with sufficient **texture** (local gradient energy above a threshold), spread
   across the strip.
3. **Phase-correlate** each window against the reference orthoimage, upsampled
   for sub-pixel precision (~10×). Implemented directly on NumPy FFTs — about
   twenty lines — to avoid adding a dependency.
4. Aggregate **robustly**: reject outliers by median-absolute-deviation, then fit
   either a constant bias or a bias plus an along-track linear drift (which is
   what an attitude drift or a clock-rate error looks like).
5. **Back-propagate into the model**, as pointing corrections, because pointing
   dominates the error budget (§9):
   - across-track ground bias → boresight **roll** bias ≈ `dx / H`
   - along-track ground bias → boresight **pitch** bias ≈ `dy / H`
6. Re-run geolocation with the corrected boresight.

**A degeneracy that must be stated rather than hidden.** From a single strip, an
along-track ground shift is produced *identically* by a pitch bias and by a clock
offset (`Δt = dy / v_ground`); a constant across-track shift is likewise
degenerate with a cross-track position error. They are not separable from this
data. **Convention: absorb along-track bias into pitch, hold the delivered
`time_sync_offset` fixed, and record the choice in the QA report.** The product is
correct either way; the attribution is a convention, and pretending otherwise
would be false precision.

**Fallback.** The reference acquisition is largely maritime, and open water has
no matchable texture. If fewer than a configurable minimum of windows match
(default 10), refinement is **skipped**, the model-only geolocation is kept, the
`refinement_skipped_insufficient_matches` flag is raised, and the model-only bias
against the STAC footprint is reported regardless. A degraded product that says
so is worth more than a confident one that is wrong.

### 8.2 Band co-registration, and self-calibrating the missing LoS

After physical orthorectification the bands should already be aligned — the
per-band view angle and the DEM have removed both the 5.7 km offset and the
terrain parallax. What remains is the residual error of the interior orientation:
exactly the per-band LoS calibration the file ships as identity (§3.3).

1. Phase-correlate each band against a reference band (**PAN** — broadest
   spectral coverage, best signal-to-noise ratio) over textured windows.
2. **Match on gradient magnitude, not raw radiance.** Cross-band correlation of
   radiance is unreliable: contrast inverts between, say, blue and
   near-infrared over vegetation and water, which biases or breaks the match.
   Edge structure is largely spectrally invariant; gradients are the robust
   feature.
3. Aggregate robustly per band → residual `(dx, dy)` in metres.
4. **Back-propagate to per-band LoS offsets** `(Δalong, Δacross)` in radians —
   populating the calibration's empty correction term.
5. Re-run the ortho with the corrected per-band LoS.

The estimated LoS offsets are written out as a **calibration artefact**: a real
geometric calibration, estimated from imagery, in the same schema the input file
left unpopulated. It can be fed back into subsequent runs of the same instrument.

**Target:** residual band-to-band registration **< 0.3 px** after correction.
Measured before and after, and reported.

---

## 9. Error budget

For the reference acquisition (altitude `H` ≈ 398 km, ground speed ≈ 7.2 km/s,
GSD 3.77 m):

| Source | Magnitude | Ground error | In pixels |
|---|---:|---:|---:|
| **Attitude knowledge** | 0.01° | **69 m** | **18 px** |
| Attitude knowledge (poor) | 0.1° | 695 m | 184 px |
| Attitude interpolation, 4 Hz, rate-aware | ~0.001° | ~7 m | ~2 px |
| **Position (no GNSS lock)** | 10²–10³ m | 10²–10³ m | 26–265 px |
| Clock offset, if not applied | 6.36 ms | 46 m | 12 px |
| Clock offset, applied | ≪1 ms | <7 m | <2 px |
| Ephemeris interpolation (Hermite) | — | <0.001 m | ~0 |
| DEM height error → absolute position | 10 m | 0.19 m | 0.05 px |
| DEM height → *band-to-band* (h = 500 m) | — | 7.1 m | **1.9 px** |
| Geolocation-grid interpolation (8 px nodes) | — | <0.4 m | <0.1 px |

**What this table says.** Attitude and position dominate by two orders of
magnitude over everything the model itself controls. Perfecting the interpolation,
the DEM or the resampling cannot fix a product whose pointing is unknown at the
0.01°–0.1° level. That is the entire justification for §8.1: **absolute accuracy
comes from the reference image, not from the telemetry.** Everything else in this
spec exists to make sure that once the pointing is corrected, nothing *else* is
wrong — and, per §3.4, to keep the bands mutually registered, which the telemetry
*can* do.

---

## 10. Architecture

Mirrors the existing layering — an abstract base class per seam, concrete
implementation behind it, domain entities between stages
(see [`architecture.md`](architecture.md)).

```text
L1B radiance (sensor coords)  +  ephemeris/attitude  +  DEM  +  reference image
        │
        ▼
  SensorModel        (line, col, band) ──► (lon, lat, h)     [§5]
        │              frames · time · ephemeris · camera · terrain
        ▼
  GeolocationGrid    subsampled lon/lat arrays per band      [§7.1]
        │
        ▼
  Warper             ──► preliminary L1C on the target grid  [§7.2]
        │
        ▼
  Refiner            relative (band → LoS)  +  absolute (→ boresight)   [§8]
        │              └── corrections fed back into SensorModel, re-run
        ▼
  L1C products  +  geometric QA report  +  STAC item  +  quicklook
```

```text
src/spp/
  geometry/
    frames.py         # WGS84, ECEF/ECI/LVLH, quaternion algebra
    timing.py         # imager ticks → UTC; per-line times; gap detection
    ephemeris.py      # Hermite (position+velocity), SQUAD (attitude+rates)
    camera.py         # detector → line-of-sight rays (intrinsics + LoS correction)
    terrain.py        # DEM + geoid abstraction; iterative ray intersection
    sensor_model.py   # composes the above: (band, line, col) → (lon, lat, h)
    geoloc_grid.py    # subsampled geolocation arrays
  resample/
    grid.py           # target grid: CRS/GSD/extent selection
    warper.py         # geolocation-array warp (GDAL)
  refine/
    matcher.py        # phase correlation on textured windows (shared)
    absolute.py       # vs reference orthoimage → boresight correction
    relative.py       # band → PAN → per-band LoS correction
  qa/
    geometric_validator.py
  pipeline/
    l1c_pipeline.py
  cli/
    run_l1c.py        # `spp-l1c`
```

New seams:

```python
class SensorModel(ABC):
    def locate(self, band: str, lines: np.ndarray, cols: np.ndarray) -> np.ndarray:
        """(N,) lines, (N,) cols -> (N, 3) lon/lat/height."""

class TerrainModel(ABC):
    def intersect(self, origin_ecef: np.ndarray, direction_ecef: np.ndarray) -> np.ndarray:
        """Iterative ray-terrain intersection -> (N, 3) ECEF ground points."""

class Refiner(ABC):
    def estimate(self, product, reference) -> Correction:
        """Measure residual displacement -> a correction to the sensor model."""
```

**New dependencies:** `pyproj` (CRS transforms, geoid), `scipy` (interpolation).
GDAL arrives with the existing `rasterio`. Phase correlation is implemented on
NumPy FFTs rather than pulling in `scikit-image`.

**Streaming.** The warp is block-wise (GDAL), and the geolocation grid is ~32 MB
per band, so memory stays bounded — consistent with the windowed L1B design in
[`performance.md`](performance.md). The refinement passes read only small windows.

---

## 11. Quality assessment

Extends `qa_report.json` with a `geometric` section:

| Metric | Meaning |
|---|---|
| `footprint.bias_m`, `footprint.corner_rmse_m` | Model-only footprint vs the delivered STAC footprint — reported **before** any refinement, so it stays an independent check |
| `absolute.n_matches`, `.bias_x_m`, `.bias_y_m`, `.residual_rmse_m` | Reference-image refinement (§8.1) |
| `absolute.boresight_roll_deg`, `.boresight_pitch_deg` | Corrections applied |
| `coregistration.<band>.residual_px` (before / after) | Band-to-band residual vs PAN (§8.2) |
| `coregistration.max_residual_px` | Headline co-registration figure; target < 0.3 px |
| `los_correction.<band>` | The self-calibrated LoS offsets (radians) |
| `dem.source`, `.coverage_pct`, `.voids_filled` | Terrain data actually used |
| `interp.max_error_px` | Geolocation-grid interpolation error (§7.1) |
| `grid.crs`, `.gsd_m`, `.native_gsd_m`, `.shape` | Output grid, and the ephemeris-derived native GSD |

**Flags:** `no_gnss_lock` (always, for this acquisition),
`refinement_skipped_insufficient_matches`, `dem_voids`, `telemetry_gap`,
`zone_straddle`, `intersection_not_converged`, `partial_band_coverage`.

---

## 12. Validation

Layered as in [`validation_strategy.md`](validation_strategy.md):

**Numerical correctness.**
- Quaternion/rotation algebra against an independent implementation (`scipy.spatial.transform`).
- Hermite and SQUAD interpolators reproduce the sampled states **exactly at the sample times**.
- Ray–ellipsoid intersection against a closed-form analytic solution.
- Round-trip: `locate(ℓ, c)` → ground → inverse-solve → recovers `(ℓ, c)` to < 0.01 px.

**Geometric correctness (end-to-end).**
- **Footprint vs the delivered STAC geometry** — the headline model check. The
  ephemeris already reproduces the footprint's scale to 0.5 % (§2.3), so a
  correct model must land within a few kilometres; anything worse means a
  convention is wrong.
- **Convention harness** (§4) — pinned by a test.
- **Independent scale check** — model-derived GSD and swath vs those implied by
  the footprint (the check of §2.3, as a regression test).

**Registration.**
- Band-to-band residuals by cross-correlation, before and after refinement;
  target < 0.3 px.
- Residuals stratified by terrain height — the discriminating test for §3.4. If
  the DEM is working, residuals must show **no correlation with elevation**. A
  residual that grows with height is the signature of a broken terrain
  correction, and a scene-average figure would hide it.

**Absolute accuracy.**
- Root-mean-square error against independent check points held out from the
  refinement fit — never against the points used to estimate it.
- Visual overlay on a basemap; coastlines are the honest test.

**Streaming correctness.** Warping in blocks equals warping whole — the same
invariance property already established for L1B.

---

## 13. Risks

| Risk | Impact | Mitigation |
|---|---|---|
| Attitude frame / quaternion convention wrong | Fatal — a plausible-looking product, hundreds of km off | §4 harness; a pinned test; the error is far too large to slip through |
| No GNSS lock | Absolute accuracy unknown | §8.1; report honestly; never snap to the footprint |
| Largely maritime scene → too few match windows | Refinement unavailable | Documented fallback + flag; PAN/near-infrared over the coastline and islands is the best available texture |
| Reference imagery unavailable, cloudy, or seasonally mismatched | Refinement degraded | Band matched by wavelength; several archives; gradient-based matching is robust to radiometric differences |
| LoS polynomial convention undefined (file is identity) | Cannot validate the correction path against the delivered file | We define the convention and document it; self-calibration (§8.2) exercises it end-to-end |
| Long strip straddles a UTM zone | Projection distortion | Centroid zone + flag + `--crs` override |
| DEM voids over water | Wrong heights | Fill with the geoid (height 0), which is the correct sea surface |

---

## 14. Implementation plan

Incremental and reviewable — one merge request per phase, with a pause for review
between them.

| Phase | Deliverable | Exit criterion |
|---|---|---|
| **0** | This spec + the convention harness (§4) | Frame, quaternion and scan conventions resolved and pinned by a test |
| **1** | Frames, timing, ephemeris, camera; sensor model on the **ellipsoid** | Computed footprint matches the STAC geometry within a few km |
| **2** | DEM + geoid + iterative terrain intersection | Intersection converges; heights validated against the geoid |
| **3** | Geolocation grid + warp → **first L1C products** (model-only) | Eight bands on one grid; band residuals *measured* |
| **4** | Relative refinement — self-calibrated per-band LoS (§8.2) | Band-to-band residual < 0.3 px, uncorrelated with terrain height |
| **5** | Absolute refinement vs reference orthoimage (§8.1) | Absolute root-mean-square error reported against held-out check points |
| **6** | Geometric QA report, STAC item, quicklook, tests, docs | QA schema complete; `spp-l1c` documented end-to-end |

Phases 1–3 produce a usable, honest product on their own: georeferenced,
orthorectified, physically co-registered, with its absolute error *measured and
declared* rather than corrected. Phases 4–5 improve the accuracy; they are not
prerequisites for a product.

---

## 15. Open questions

1. **Is `frame 1` really LVLH?** Two independent numerical arguments say yes
   (§3.1); §4 settles it. Highest-impact unknown in the level.
2. **Does the LoS polynomial run in normalized or raw detector column, and to
   what degree?** Six coefficients per axis suggests degree 5. The delivered file
   is identity, so the convention cannot be inferred from it — we define one and
   document it.
3. **Should the along-track bias be attributed to pitch or to the clock?** Not
   separable from one strip (§8.1). Convention chosen; revisit with a
   cross-strip or a manoeuvre.
4. **Is `PAN` the right co-registration reference?** Best signal-to-noise ratio
   and broadest spectral overlap, so yes by default — but its 625 nm centre is
   spectrally distant from blue, which is where the match is weakest. A red or
   green reference for the short-wavelength bands may prove more robust; decide
   with the Phase 3 residuals.
```

