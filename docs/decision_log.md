# Decision Log

Key engineering decisions taken while building this framework. For each:
alternatives considered, what was chosen, and the trade-offs.

---

## 1. Which level to implement end-to-end → **L1B (DN → TOA radiance)**

**Alternatives.**
- **L1A** (count reconstruction) — but the package already ships per-band DN
  rasters and the raw session binary is not included, so there is nothing to
  reconstruct.
- **L1C** (georeferencing) — depends on external data (DEM) and the full
  ephemeris/attitude chain; it also cannot be *absolutely* validated from the
  package alone (no GNSS lock, no reference orthoimage shipped).
- **L2A** (surface reflectance) — the most external-data-intensive level
  (atmospheric state, aerosol model), with the least in-package validation.

**Choice.** L1B. It is **self-contained** (every required asset — CPF, filters,
solar, telemetry — is in the package), produces the **first physically
meaningful** product (radiance in SI units), is the **foundation** every higher
level builds on, **exercises the provided calibration assets**, and is **fully
verifiable** from the delivered data.

**Trade-offs.** It does not tackle the hardest geometric/atmospheric problems —
those are instead specified in detail in [`remaining_levels.md`](remaining_levels.md).

---

## 2. Treat the "L0" package as image-assembled (L0/L1A boundary)

**Alternatives.** Attempt true packet-level reconstruction from `raw.bin`.

**Choice.** Start from the per-band rasters. The package is a hybrid: it delivers
assembled per-band DN rasters plus the calibration assets normally consumed by
the next stage. `raw.bin` is intentionally not shipped.

**Trade-offs.** True L0→L1A is out of scope; documented as such so the scope
boundary is explicit.

---

## 3. Per-line temperature interpolation (vs a single mean)

**Alternatives.** Use one mean detector temperature for the whole scene.

**Choice.** Interpolate the telemetry samples to a per-line temperature profile
and evaluate the darkfield per line.

**Trade-offs.** The detector temperature swings ~17 °C across the acquisition;
with a thermal gradient of ~0.45 DN/°C that is ~7.6 DN of dark drift — a
systematic, along-track error a single mean would bake in. Interpolation removes
it. It assumes the telemetry samples are uniformly spaced in time (see
[`limitations.md`](limitations.md) for the exact-timing refinement).

---

## 4. NaN NoData, and **no clipping** of negative radiance

**Alternatives.** Clip negative radiance to 0; or use a sentinel (e.g. -9999).

**Choice.** Mask NoData as NaN and preserve slightly-negative radiance from
near-dark pixels.

**Trade-offs.** Clipping would look cleaner but hides dark-subtraction bias;
preserving negatives lets QA *measure* it. NaN is the natural float NoData but
can be awkward for some GIS tools — a sentinel option is noted as future work.

---

## 5. Windowed streaming I/O (vs whole-array processing)

**Alternatives.** Read each band fully into memory, calibrate, write.

**Choice.** Process line-windows; read DN, calibrate, write and QA per block.

**Trade-offs.** More code (block loop, `line_start` plumbing, accumulators) in
exchange for bounded memory: a band is ~250 MB as `uint16` and ~1 GB as
`float64`; eight at once is not viable. Streaming keeps RAM flat and was verified
to be window-invariant (identical result regardless of window size).

---

## 6. Separation of concerns; pure-numpy calibrator; lazy bands

**Alternatives.** A single calibration script doing I/O + math + writing.

**Choice.** Distinct reader / calibrator / QA / writer / pipeline modules. The
calibrator is pure (numpy in, numpy out, no I/O); `Band` is lazy (path + raster
properties, never pixels); the pipeline owns orchestration only.

**Trade-offs.** More files and indirection, but each stage is independently
testable, the calibration math is reusable outside any I/O context, and the core
entity layer carries no raster-library dependency.

---

## 7. `float32` radiance output

**Alternatives.** `float64`.

**Choice.** `float32`. Verified against an independent `float64` reference: max
difference 2.3 × 10⁻⁵ (rounding only) — negligible for radiance, and it halves
output size and memory.

---

## 8. GeoTIFF output in sensor coordinates (no CRS), tiled + overviews + predictor

**Alternatives.** NetCDF; a single multiband file; no overviews.

**Choice.** One COG-friendly (Cloud-Optimized GeoTIFF) GeoTIFF per band: tiled
(512), compressed with a floating-point predictor, internal overviews, **no CRS**
(L1B is still in sensor coordinates).

**Trade-offs.** GeoTIFF + GDAL ecosystem is the most interoperable choice for
downstream geometric processing; per-band files keep the model simple and match
the input layout. NetCDF/Zarr would suit a multi-dimensional cube better and is
noted as future work.

---

## 9. Mission-agnostic framework with externalised, discovered configuration

**Alternatives.** Hard-code band names/files and sensor constants.

**Choice.** Discover bands from the STAC assets, calibration files by glob, and
the band name↔id map from the filter file; keep the code free of mission/sensor
proper nouns.

**Trade-offs.** Slightly more discovery logic, but the framework reads as a
reusable EO tool and adapts to a new sensor through data, not code edits (see
[`generalisation.md`](generalisation.md)).

---

## 10. Emit a STAC item for the L1B product, with `null` geometry

**Alternatives.** Ship only the rasters + a QA report; or fabricate a footprint
on the output item by carrying the source scene's `bbox`/`geometry`.

**Choice.** Write a STAC 1.0.0 item per product (eo/raster/processing
extensions): band assets with spectral/raster properties, processing lineage
(software version, `derived_from` the source scene), quicklook and QA as assets.
The geometry is `null` because L1B is still in sensor coordinates.

**Trade-offs.** A real footprint would make the item map-discoverable, but L1B is
**not georeferenced** — asserting a geometry would be misleading; an honest
`null` geometry defers that to L1C. The item adds negligible runtime and makes
the product catalogue-ready and interoperable; it is on by default (`--no-stac`
to skip).

> **Documented alternative (discoverability).** The source package's STAC item
> *does* carry a **nominal scene footprint** (the planned `bbox`/`geometry`).
> That outline could be propagated to the L1B item to make the product findable
> on a map, **provided it is explicitly labelled as nominal** — a planning-grade
> outline, not a pixel-accurate footprint — since the assets remain in sensor
> coordinates and this scene has no GNSS lock (so even the nominal outline is
> approximate). This is a legitimate choice when catalogue discoverability is
> prioritised over strict rigour; it would be a deliberate, documented toggle,
> never a silent default. We keep `null` here as the more rigorous default.
