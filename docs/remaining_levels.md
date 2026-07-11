# Specification of the Remaining Levels

This document specifies the processing levels **not** implemented in this
framework, in enough detail for a colleague to implement them: algorithm,
required inputs (with emphasis on **external data**), key challenges and a
validation strategy. The implemented level (L1B) is documented in
[`product_hierarchy.md`](product_hierarchy.md) and the code under
[`../src/spp/`](../src/spp/).

Levels covered:

- [L1A — detector-count reconstruction](#l1a--detector-count-reconstruction)
- [L1C — georeferencing / orthorectification](#l1c--georeferencing--orthorectification) *(external data: DEM, ephemeris, attitude)*
- [L2A — surface reflectance](#l2a--surface-reflectance) *(external data: solar, atmospheric state)*

---

## L1A — detector-count reconstruction

**Goal.** Turn the raw session binary into one `uint16` DN raster per band in
sensor coordinates. *(In the provided package this is already done — the bands
ship as per-band rasters — so this is specified for completeness and for a true
L0 ingest.)*

**Inputs.**
- `raw.bin` — the raw session binary (detector readout stream + housekeeping).
- `ImagerConfiguration` — line period, `SpectralBands`, per-band TDI (Time-Delay Integration)
  (`BandSetup`), per-band detector `BandStartRow`, `ScanDirection`,
  `BinningFactor`.
- `SensorConfiguration` — ADC (Analogue-to-Digital Converter) range, readout offsets, gain, e-black flag.

**Algorithm.**
1. Parse the binary container; validate frame/packet CRCs (cyclic redundancy checks); recover the readout
   line sequence with timestamps.
2. Decompress the payload if compressed (a compression flag in the product
   metadata indicates the codec; `0` = none).
3. For each readout line, slice the detector frame; for each band, take its
   window `[start_row : start_row + tdi)` of the `4096`-wide detector axis.
4. Accumulate the `tdi` lines into a single output line per band (TDI summation),
   respecting `ScanDirection`.
5. Apply binning if `BinningFactor > 0`.
6. Assemble per-band rasters; mark dropped/duplicated lines and known dead pixels
   as NoData (`0`).

**Key challenges.**
- Binary format details (endianness, packet boundaries, sync markers).
- TDI summation can overflow 16-bit intermediates — accumulate in wider integers.
- Timing gaps / dropped lines create along-track discontinuities that must be
  flagged, not silently closed.
- Dead/hot pixel maps must be applied consistently with the calibration's pixel
  ordering.

**Validation.**
- Round-trip: regenerate the per-band rasters and compare to a reference L1A
  product (bit-exact where lossless).
- Line count vs expected `acquisition_duration / line_period`.
- Per-band histograms within the ADC dynamic range; NoData only where expected.

---

## L1C — georeferencing / orthorectification

> **Superseded.** This section is the original outline. L1C now has a full
> implementable design in [`l1c_spec.md`](l1c_spec.md) — methods, error budget,
> QA schema and a phased implementation plan. Where the two differ, the spec
> wins: it was written against the *verified* contents of the package, and it
> corrects two assumptions made below (the geometric calibration is a placeholder,
> not a populated per-column line-of-sight table; and the nominal ground sample
> distance and swath do not describe the reference acquisition).

**Goal.** Resample each L1B radiance band onto a projected map grid (defined CRS,
e.g. UTM), co-registering bands to a common ground geometry. **This is the
external-data-dependent geometric level.**

**Inputs.**
- L1B radiance bands (sensor coordinates).
- **Ephemeris** — platform position/velocity in ECEF — Earth-Centred, Earth-Fixed — (`ancillary.extrinsics.hist`,
  frame 3), time-tagged.
- **Attitude** — body-frame quaternions (`ancillary.extrinsics.hist`, frame 1),
  time-tagged.
- **Camera intrinsics** — focal length, pixel pitch, principal point, sensor
  size, body→detector reference frame (`ancillary.intrinsics`).
- **Geometric calibration (CPF — Calibration Parameter File)** — `boreSightAlignment` (3×3) and per-band
  `lineOfSight` arrays (`along`, `across`, radians per detector column).
- **DEM** — external elevation model (e.g. Copernicus GLO-30) over the footprint.
- Earth model (WGS84) and a time system (UTC ↔ GPS/leap seconds) for ephemeris
  interpolation.

**Algorithm.**
1. **Line timing.** Assign each raster line `ℓ` a UTC time from the acquisition
   start, line period and `ScanDirection` (or per-line timing if available).
2. **State interpolation.** Interpolate ephemeris (position/velocity) and
   attitude to each line time — linear/Hermite for position, SLERP (spherical linear interpolation) for
   quaternions.
3. **Viewing ray.** For detector column `c`, form the line-of-sight unit vector
   from `(along[c], across[c])`; rotate by `boreSightAlignment`, then by the
   body→ECEF attitude rotation, to get the ray direction in ECEF.
4. **Ground intersection.** Intersect the ray from the platform position with the
   DEM-bearing Earth surface. Iterate: intersect the ellipsoid, sample the DEM at
   that point, re-intersect at the new height until convergence (terrain
   correction).
5. **Sensor model.** The set of per-`(ℓ, c)` ground coordinates defines a
   sensor-to-ground mapping. Densify it into a GCP (Ground Control Point) grid or fit RPCs (Rational Polynomial Coefficients).
6. **Resampling.** Warp each band onto the target grid/CRS at a chosen GSD (Ground Sample Distance)
   (inverse mapping: for each output pixel, find the source `(ℓ, c)` and
   interpolate). Each band uses its **own** LoS (and start-row time offset), which
   is what co-registers the bands.

**Key challenges.**
- **No GNSS (Global Navigation Satellite System) lock** in this acquisition (`extrinsics.gnss_lock = false`) →
  degraded absolute geolocation; expect a bulk offset and plan for GCP/tie-point
  refinement against a reference image.
- **Time synchronisation** — a small clock offset (`time_sync_offset` in STAC)
  shifts the whole strip along-track; it must be applied.
- **Band-to-band co-registration** — bands are read from different detector start
  rows, so they image the same ground point at slightly different times/angles;
  parallax over terrain must be handled via the per-band LoS + DEM.
- **Attitude quality** dominates the error budget for a long focal length; small
  pointing errors map to large ground displacements.
- **Tall, narrow strip** (≈4096 × 30948) — resampling and DEM tiling must stream;
  projection choice should minimise distortion over the long along-track extent.

**Validation.**
- Tie points / GCPs against a reference orthoimage; report RMSE in metres.
- Band co-registration error (cross-correlation between bands) at sub-pixel level.
- Computed footprint vs the STAC `bbox`/`geometry`.
- Visual overlay on a basemap; check coastlines/known features.

---

## L2A — surface reflectance

**Goal.** Convert L1C radiance to **surface reflectance** (bottom-of-atmosphere)
by removing illumination and atmospheric effects. **This is the most
external-data-intensive level.**

**Inputs.**
- L1C georeferenced radiance.
- **Solar reference** — `spectral_solar.json` (reference solar spectral
  irradiance) and `temporal_solar.json` (annual Earth–Sun distance model).
- Per-band **spectral response** — `filters_payload_0.json` (`wavelength`,
  `relative`).
- Sun/view geometry — solar zenith/azimuth and view zenith/azimuth per pixel
  (from acquisition time + L1C geometry).
- **Atmospheric state** — aerosol optical depth, water vapour, ozone, surface
  pressure (from ancillary sources such as CAMS/MODIS, or retrieved from the
  scene); an aerosol model.
- A radiative-transfer engine / LUT (e.g. 6S via Py6S; ACOLITE for water).
- DEM (surface pressure / altitude correction).

**Algorithm.**
1. **Band solar irradiance (ESUN — exo-atmospheric solar irradiance).** For each band, weight the reference solar
   spectrum by the band's relative spectral response and integrate:
   `ESUN_b = Σ(E_sun(λ)·rsr_b(λ)) / Σ(rsr_b(λ))`. Scale by the temporal solar
   model evaluated at the acquisition date (Earth–Sun distance `d`).
2. **TOA reflectance.** `ρ_TOA = π · L · d² / (ESUN_b · cos θ_s)`, with `L` the
   L1C radiance and `θ_s` the solar zenith angle.
3. **Atmospheric correction.** Run the radiative-transfer (RT) model with the
   atmospheric state and sun/view geometry to obtain path radiance, atmospheric
   transmittance (up/down) and spherical albedo; invert to surface reflectance
   `ρ_s = f(ρ_TOA; xa, xb, xc)`. For aquatic scenes, use a water-optimised scheme
   (SWIR (Short-Wave Infrared)/NIR-based aerosol estimation; retrieve water-leaving reflectance).
4. Mask clouds/cloud-shadow/sun-glint; write per-band surface reflectance.

**Symbols (every term in the equations above).**
- `b` — spectral band index; `λ` — wavelength.
- `E_sun(λ)` — reference solar spectral irradiance at `λ` (from
  `spectral_solar.json`).
- `rsr_b(λ)` — relative spectral response of band `b` at `λ` (from
  `filters_payload_0.json`).
- `ESUN_b` — band-integrated exo-atmospheric solar irradiance for band `b`; the
  response-weighted average of `E_sun` over the band. Express it in the **same
  spectral unit as `L`** (per µm — i.e. the per-nm solar integral × 10³) so
  `ρ_TOA` comes out unitless.
- `d` — Earth–Sun distance (astronomical units) at the acquisition date (from
  `temporal_solar.json`); the `d²` term normalises for the seasonal Sun distance.
- `L` — at-sensor (L1C) spectral radiance, `W·m⁻²·sr⁻¹·µm⁻¹` (the L1B output unit).
- `θ_s` — solar zenith angle; `cos θ_s` corrects for the illumination geometry.
- `π` — converts radiance (per steradian) to a (Lambertian) reflectance.
- `ρ_TOA` — top-of-atmosphere (apparent) reflectance, unitless in `[0, 1]`.
- `ρ_s` — surface (bottom-of-atmosphere) reflectance — the L2A output, unitless.
- `xa, xb, xc` — the atmospheric-correction coefficients from the RT model
  (6S convention): `xa` the inverse total transmittance/gain, `xb` the
  path-radiance term, `xc` the atmospheric spherical albedo — combined as
  `y = xa·L − xb`, then `ρ_s = y / (1 + xc·y)`.

**Key challenges.**
- **Aerosol retrieval** is the dominant uncertainty, especially over bright land
  or sun-glint-prone water.
- **Sun-glint and adjacency** effects over water; thin cirrus contamination.
- **Gas absorption** in red-edge/NIR (water vapour) requires accurate column
  amounts.
- **Cloud masking** without a thermal band — rely on spectral tests + thresholds.
- **Validation data scarcity** — coincident in-situ/AERONET matchups are rare.

**Validation.**
- Compare to a coincident reference surface-reflectance product (e.g. Sentinel-2
  L2A) over stable targets; report per-band bias/RMSE.
- AERONET-based atmospheric correction validation where available.
- Physical sanity: near-IR surface reflectance over deep clear water ≈ 0;
  spectral shape of known targets (vegetation red-edge, soil).

---

## References

See [`references.md`](references.md). Key external resources: Copernicus DEM,
GDAL warp/VRT, Py6S / 6S, ACOLITE, STAC, CEOS (Committee on Earth Observation Satellites) product-level definitions.
