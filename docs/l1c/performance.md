# L1C — Performance

Measured on the reference acquisition (4096 × 30948 per band, eight bands), on the same
machine as the L1B figures above.

| Stage | Cost | Notes |
|---|---|---|
| Geolocation grid | **~3.8 s** per band | 513 × 3870 lattice at step 8; includes the terrain intersection at every node (39% of that time) |
| Lattice in memory | **32 MB** per band | `float64` longitude and latitude; `float32` would resolve longitude to only ~0.7 m |
| Warp | **~14 s** per band | Bilinear, through the geolocation array — **was 46 s** |

### Where the time actually went

Profiled rather than guessed, because the answer was not where the interesting code is.
**The warp was 93% of the per-band cost** (46.5 s of 50 s); the geolocation grid — every
line of physics in this level — was 3.8 s. Two defaults were doing the damage, and neither
involved a trade-off:

| | before | after |
|---|---:|---:|
| Resampling (`reproject`) | 24.6 s, **one core** | **6.5 s**, ten cores |
| Write + overviews | 20.4 s (Deflate + predictor) | **5.7 s** (Zstandard level 1) |

GDAL's warper is single-threaded unless told otherwise, and it had not been told. And
Zstandard produces a file the **same size** as Deflate (432 MB against 434 MB) in a quarter
of the time — Deflate was costing 15 seconds a band for nothing at all, because the
floating-point predictor is what is doing the compressing, before the codec ever sees the
data.

**Per band: 50 s → 18 s.**

### And then the profile moved

Optimising the warp was right and it was **incomplete**, because the profile above was
taken with every refinement stage switched off. With them on — which is how the pipeline
actually runs — the picture is different:

| Stage | Time (4 bands) |
|---|---:|
| **Absolute pointing correction** | **57 s** |
| Resampling (4 bands) | ~56 s |
| Geolocation grids (4 bands) | 14 s |
| Self-calibration | 8 s |
| Along-track drift | 4 s |
| Telemetry parities | 2 s |

The absolute correction now costs as much as the entire warp. It had never been profiled,
because it had been disabled in the run that was.

Most of that was **not arithmetic — it was the network**. The reference is a virtual mosaic
over *remote* Cloud-Optimized GeoTIFFs, and the coarse-to-fine refinement calls the
estimator six times (three lattices, each re-measured to prove it improved). Every one of
those was re-fetching the same few hundred megabytes of Sentinel-2. The patch depends only
on the footprint and the footprint does not move by more than metres during the
refinement, so it is now read **once** and kept. That halved the stage.

The lesson is the one the whole profile teaches: **a benchmark that disables the expensive
parts measures a program you are not running.**

### What was tried and reverted

The per-band rasters are an intermediate — the stack is assembled from them and they are a
duplicate of it. Writing the stack **directly from memory**, band by band as each comes off
the resampler, would skip a compress, a decompress and two passes over 1.7 GB.

It was tried. It produced a **stack that was entirely NaN**: every band resampled correctly
(the log said 39% coverage), the file was 652 MB, it had the right georeferencing and the
right band names, and every pixel read back as NoData. The failure does not reproduce in
isolation — the same open-stack, band-interleaved, write-band-by-band pattern works at full
scale in a test — and I could not pin down what defeats it.

It is not used. This pipeline has already shipped one silently-empty stack, and an
optimisation whose failure mode cannot be explained is not worth 40 seconds when the failure
is invisible.

A second attempt — having the resampler *return* the array instead of writing it — was
killed by the operating system: the previous band's 1.15 GB destination stays alive while
the next is allocated. **The disk round trip is cheaper than the memory.**

What was kept: the per-band rasters skip their overviews (they are an intermediate, nobody
looks at them) and are deleted unless `--per-band` is passed. And the cached reference patch
— 416 MB, and 208 MB now that it is `float32` rather than `float64` for 16-bit integer data
— is **freed before the resampling starts**, because holding it alongside a 1.15 GB
destination is what got the process killed in the first place.

**187 s → 144 s on a four-band run.**

The remaining large win is *band-level parallelism* — the bands are independent — but the
warp holds ~1.7 GB per band (see below), so four in flight would want 7 GB. That is a
memory problem to fix first, not a speed problem to exploit.
| Output raster | 9774 × 29403 px @ 4.0 m | 1.15 GB uncompressed, **435 MB** on disk (Deflate + floating-point predictor) |
| Elevation model | remote, cacheable | Copernicus GLO-30 over `/vsicurl`, no authentication; ~40% of the strip is open water |

Output coverage is **40%** of the raster: a 116 km strip inclined ~11° from north,
placed in a north-up projected box, is mostly NoData by area. That is inherent to a
projected single-strip product, not a defect, and the empty region compresses to almost
nothing.

### Memory: an honest regression against L1B

L1B processes a strip in **O(window)** memory — peak RAM is constant regardless of strip
height, which is the property the rest of this document is about.

**The L1C warper does not have that property.** It reads the whole source band into
memory (~507 MB as `float32`) and allocates the whole destination (~1.15 GB), so peak
usage is **~1.7 GB per band** and it grows with the strip. An earlier draft of
[`l1c_spec.md`](spec.md) claimed the warp was block-wise and bounded; that claim was
wrong, and it is recorded here rather than quietly dropped.

It is workable at this scale and it is not acceptable at the next one. The fix is
known — warp into an on-disk destination and let GDAL stream blocks through it, which
is how the geolocation-array path is meant to be driven — and it is the first thing to
address before this level meets a strip several times longer.
