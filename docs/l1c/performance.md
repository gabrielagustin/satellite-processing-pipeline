# L1C — Performance

Measured on the reference acquisition (4096 × 30948 per band, eight bands), on the same
machine as the L1B figures above.

| Stage | Cost | Notes |
|---|---|---|
| Geolocation grid | **~3.5 s** per band | 513 × 3870 lattice at step 8; includes the terrain intersection at every node |
| Lattice in memory | **32 MB** per band | `float64` longitude and latitude; `float32` would resolve longitude to only ~0.7 m |
| Warp | **~44 s** per band | Bilinear, through the geolocation array |
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
