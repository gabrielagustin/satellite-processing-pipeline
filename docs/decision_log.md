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
  package alone (no GNSS — Global Navigation Satellite System — lock, no reference orthoimage shipped).
- **L2A** (surface reflectance) — the most external-data-intensive level
  (atmospheric state, aerosol model), with the least in-package validation.

**Choice.** L1B. It is **self-contained** (every required asset — CPF (Calibration Parameter File), filters,
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
it. How each sample maps to a line is itself a decision — see §11.

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

---

## 11. Map temperature samples to lines by timestamp (not uniform spacing)

**Alternatives.** Spread the telemetry samples **uniformly** across the lines
(assume equal spacing and that the telemetry exactly spans the acquisition).

**Choice.** Place each sample at its true line from its `ImagerTime` timestamp,
anchoring the line clock to the imager clock via the `TimeSync` block (an
`ImagerTime`↔platform-epoch tie, confirmed by 1 Hz PPS (Pulse Per Second) pulses) and the
acquisition start time; line `ℓ` is then at `t_line0 + ℓ · line_period`. The
uniform model remains as a **fallback** when timing inputs are missing, and the
chosen model is recorded in provenance.

**Trade-offs.** More parsing (TimeSync, epoch conversion) and a sensor-specific
timing assumption in the reader — but it removes a real, measured error. On the
reference scene the telemetry **brackets** the imaging (≈19 s of telemetry around
a 16 s acquisition) and the temperature profile is non-monotonic, so the uniform
model mis-assigns temperature by up to ~15 °C, biasing per-line radiance by up to
**~3.5 %** along-track (a banding gradient) while leaving the band **mean**
almost unchanged — invisible to scalar QA. The timestamped model corrects it.
Verified on the real scene (per-line radiance comparison) and covered by unit
tests (`tests/test_package_reader.py` for the timestamp→line mapping,
`tests/test_l1b_calibrator.py` for the interpolation).

---

## 12. Keep both `requirements.txt` and `pyproject.toml`

**Alternatives.** Ship only one: a bare `requirements.txt` (as the challenge
asked for), or only `pyproject.toml`.

**Choice.** Keep both, with `pyproject.toml` as the **packaging source of
truth** — it declares the runtime dependencies, the `test` extra and the
`spp-l1b` console script, and enables an editable install (`pip install -e .`).
`requirements.txt` is retained because the challenge requested it, kept as a thin
**mirror of the runtime dependencies** for convenience and for tooling that
expects one.

**Trade-offs.** Two dependency lists to keep in sync, but the duplication is
small (three runtime packages) and bounded: the `test` extra and the entry point
live only in `pyproject.toml`, so `pip install -e .` alone is sufficient to
install and run. `requirements.txt` mirrors the full runtime set (`numpy`,
`rasterio`, `matplotlib`), so either install path yields the same working tool.

---

## 13. Quicklook on by default; `matplotlib` a core dependency

**Alternatives.** Keep the RGB quicklook opt-in (`--quicklook`) with `matplotlib`
as an optional `viz` extra — lighter base install, but a default run produces no
visual output.

**Choice.** Produce the quicklook **by default** (`--quicklook` /
`--no-quicklook`, mirroring `--stac`), and promote `matplotlib` to a **core
dependency**. A first run out of the box then yields a human-readable RGB PNG
alongside the rasters, QA report and STAC item — what a reviewer most wants to
see — with no extra install step or flag.

**Trade-offs.** Every install now pulls `matplotlib` (a non-trivial dependency)
even for purely programmatic use, and it reverses the earlier "matplotlib is
optional" framing. Judged worth it: the quicklook is part of how the product
demonstrates correctness, the cost is one well-established package, and the
default stays overridable (`--no-quicklook`). The quicklook still skips
gracefully (with a warning) when the R/G/B bands are absent.

---

## 14. Resolve the telemetry's undocumented conventions by measurement, not by reading

**Context.** The acquisition package delivers platform position, velocity and
attitude tagged with numeric frame identifiers it never defines, and a
detector-to-body matrix that does not say which way round the detector's two axes
feed into it. **Eight** properties of the telemetry are therefore ambiguous:
the ephemeris frame, the attitude frame, the quaternion's component order and
rotation direction, the scan direction, and the detector column axis and the sign
of each detector axis.

**Alternatives.** Read the conventions off the vendor documentation and hard-code
them. Rejected: the documentation does not state them, and an assumption made in
their place is invisible in the code and catastrophic in the product. Our own
first reading assumed the ephemeris was Earth-fixed and put the sub-satellite
point **5,000 km** from the delivered footprint — a wrong answer that no amount of
downstream care would have recovered from, and one that looked entirely plausible
until it was measured.

**Choice.** Make every ambiguous convention an explicit parameter
(`geometry.Convention`), enumerate all 384 combinations, and **measure** which one
the data supports, using two probes that answer different questions:

- **Footprint match** against the delivered catalogue geometry. Resolves the
  frames: a wrong frame throws the footprint 700 km or more.
- **Band coherence** — locate the same raster sample through two bands on widely
  separated detector rows and measure how far apart they land. This is the sharp
  probe, because an absolute bias displaces both bands identically and **cancels**.

The two are deliberately **not** combined into one score. Without a GNSS lock even
the correct model misses the delivered footprint by ~30 km (see
[`l1c_spec.md`](l1c_spec.md) §3.2), so summing them would let that irreducible
bias swamp the hundreds-of-metres signal that separates the fine conventions. The
footprint rejects gross failures; coherence ranks the survivors. Coherence fell
from 8,449 m under the initial reading to **123 m** under the resolved one.

**Resolved:** ephemeris frame = **inertial** (not Earth-fixed); attitude frame =
**LVLH**; quaternion = **scalar-first**; detector column axis = **y**; detector row
sign = **negative**.

**Trade-offs.** Three conventions remain **unresolved**, and the harness reports
them as such rather than picking the least-bad. Scan direction and detector column
sign are *parities* — they mirror the strip north–south and east–west, mapping the
footprint's corners onto each other and displacing both bands equally, so both
probes are blind to them by construction. Quaternion direction is preferred 5:1 by
coherence, short of the 10× margin required. All three need **image content** — the
reference-image matching of `l1c_spec.md` §8.1 — and until then they are carried as
declared assumptions, flagged in the quality report.

Reporting an unresolved axis as resolved would be manufacturing a result, and
would be worse than the honest gap: a product that is silently mirrored is harder
to catch than one that says it might be.

---

## 15. Check the telemetry's derivatives before using them

**Context.** The platform telemetry delivers derivatives alongside the states —
velocity with position, angular rates with attitude. Cubic Hermite interpolation
can use them, and at a few hertz the accuracy is worth having, so
[`l1c_spec.md`](l1c_spec.md) originally specified rate-aware interpolation for
both on the grounds that "the rates are delivered, so using them is free".

**What happened.** They are not free. The delivered angular rates are **~55x
larger** than the rotation the quaternion sequence itself undergoes between
samples: they describe the body's motion in some other frame, not the drift of the
small attitude offset the quaternions encode. Used as Hermite endpoint slopes they
inject rotation the attitude never contained, and the band-to-band error rose from
**18 m to 123 m** — a "refinement" that made the model seven times worse and looked,
until it was measured, like an irreducible platform bias.

The failure mode is not the obvious one. A slope of merely the wrong *magnitude*
largely cancels: the two Hermite endpoint weights sum to `s(2s-1)(s-1)`, which
vanishes at the midpoint. The damage comes from a rate pointing along a **different
axis**, which is exactly what is delivered.

**Choice.** Do not add a flag; add a **check**. `Attitude.rates_are_consistent`
compares the delivered rates against the rates implied by the quaternion sequence,
and interpolation falls back to SLERP — automatically, with a warning — whenever
they disagree. Asking for rate-awareness is a request, not an instruction: **the
data gets a veto**, because the caller cannot know what a given payload contains.
Position keeps Hermite: its velocity *does* agree with its positions, and the check
passes.

**Trade-offs.** SLERP discards curvature that a correct rate would capture, so a
payload whose rates are honest is interpolated slightly less well than it could be
— but only until the check passes, at which point Hermite is used. The cost is one
cheap comparison per acquisition; the alternative is a silently wrong model.

**Generalisation.** The same discipline caught the neglected precession
(entry 16). An unvalidated assumption about the input does not announce itself: it
presents as an irreducible error, and it will be attributed to the platform rather
than to the model unless someone asks where the residual points.

---

## 16. Model Earth rotation with precession and nutation, not sidereal time alone

**Context.** The sensor model must rotate the platform state from its inertial
frame into an Earth-fixed one. The first implementation applied only the Greenwich
Mean Sidereal Time angle, documenting the omission as "precession, nutation and
polar motion are neglected — they are arcsecond-level effects, far below the
pointing uncertainty that dominates this level's error budget".

**What happened.** That justification was wrong by four orders of magnitude.
Precession accumulates at ~50 arcseconds **per year**, so a quarter-century after
J2000 it is a third of a degree — **tens of kilometres** at orbital radius, not
arcseconds. It left a 32.8 km footprint error that was indistinguishable from the
platform bias the spec predicted for a GNSS-less acquisition, and would have been
absorbed into the "absolute refinement" stage as if it were real.

It was caught by asking which *direction* the residual pointed. The offset was
−33.9 km east and −16.5 km north; precession in right ascension over the elapsed
26.2 years predicts −33.5 km, and in declination −16.2 km. Both matched to 2%.

**Choice.** Implement the full IAU-76/FK5 chain: precession → nutation → apparent
sidereal time. Footprint error fell from **32.8 km to 208 m**. Polar motion (metres)
and the truncated nutation terms (hundreds of metres) are still neglected — but now
that omission is *justified* rather than asserted, being far below the geolocation
error the level carries without a GNSS lock.

**Trade-offs.** ~60 lines of standard astrodynamics and a handful of magic
polynomial coefficients, against a 150x reduction in geolocation error. No contest.
The residual 208 m is model-versus-model agreement, not absolute accuracy: the
delivered footprint was almost certainly derived from the same telemetry, so it
validates the implementation, not the orbit.

---

## 17. Demand a real geoid; refuse PROJ's silent "ballpark" fallback

**Context.** A digital elevation model reports **orthometric** height (above the
geoid); the platform ephemeris is referenced to the **ellipsoid**. Over the
reference acquisition those differ by about **−31 m**. The conversion is a standard
PROJ vertical transformation.

**The trap.** Asked to convert without the geoid grid installed, PROJ does **not**
raise. It performs what it calls a "ballpark vertical transformation", which returns
the height unchanged — an undulation of exactly zero, everywhere. Nothing downstream
notices: heights are finite, rays converge, the product looks correct, and every
elevation is wrong by the local undulation. A silent, uniform, tens-of-metres error
is exactly the kind that survives to production.

**Choice.** Enable PROJ's network grid access, and **validate the result**:
`terrain.check_undulation` rejects an undulation that is identically zero or
non-finite, with an error message naming the cause and the fix. The geoid is not
trusted because it was requested; it is trusted because it was checked.

**Trade-offs.** A first run needs network access to fetch the grid (or a one-off
`projsync`), and the check would false-positive on a scene where the true undulation
really is zero everywhere — a set of measure zero on the actual geoid, and a case
that would still deserve a second look.

**Related.** This is the same discipline as entries 15 and 16: an unvalidated
assumption about the input does not announce itself. It presents as a plausible
result, and only a check that could have failed distinguishes the two.

---

## 18. Write exit criteria that can fail

**Context.** The L1C level was planned in phases, each with an exit criterion. The
Phase 3 criterion was written as: *"eight bands on one grid, and the band residual
has fallen below the 1–7 px it starts at. This is the test that the model is real:
if the residual does not shrink, stop — the model is wrong, and no amount of
refinement will save it."*

**Alternatives.** The obvious criterion — *"the products are produced"* — would have
been met. The geolocation grid built, the interpolation error came in at 0.026 px
against a 0.1 px budget, every node found the ground, and eight bands landed on one
projected grid. By any check of the machinery, the phase succeeded.

**What happened.** The criterion failed. Warping through the physical model left the
band residual exactly where it started: 30.2 m before, 33.2 m after. The model turned
out to carry a per-band error of its own — a constant (−5.1, +6.5) px — which is the
*same size* as the misregistration it was meant to remove.

**Why the criterion earned its place.** Its real job was to distinguish a **missing
calibration** from a **modelling bug**, because feeding the residual back into the
model (the Phase 4 self-calibration) would paper over either one indiscriminately. The
discriminator is whether the residual varies with position: correlations with
along-track and cross-track position came out at −0.09, −0.05, +0.01 and −0.17 — none.
A residual constant across the whole strip and the whole swath is a fixed angular
offset per band, which is the interior orientation, which is precisely the line-of-sight
term the calibration file ships unpopulated.

So the self-calibration is legitimate — but the **phasing was wrong**. Phase 4 is not a
refinement that polishes good co-registration into sub-pixel co-registration; it is the
prerequisite for co-registration existing at all. That is now recorded in
[`l1c_spec.md`](l1c_spec.md) §14.1.

**Choice.** Keep writing criteria that can fail, and state the *number* they must beat
rather than the artefact they must produce. A criterion phrased as an output ("the file
exists", "the pipeline runs", "the products are produced") cannot distinguish a working
model from a broken one — it can only distinguish a working model from a crash.

**Trade-offs.** A failing criterion costs a phase's worth of schedule and forces the
plan to be rewritten in public. That is the price, and it is small: the alternative was
shipping a product whose bands do not align, discovering it somewhere downstream, and
having no record of which stage was responsible.

**Related.** The same discipline, applied to inputs rather than phases, is entries 15,
16 and 17. The full narrative is in [`l1c_findings.md`](l1c_findings.md).
