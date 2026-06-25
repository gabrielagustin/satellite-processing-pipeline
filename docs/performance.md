# Processing Large Strips: Memory & Performance

A pushbroom acquisition is a tall, narrow strip — here **4096 × 30948 pixels per
band (~127 M pixels)**, eight bands. Loading a band in full is ~240 MB as
`uint16` and ~960 MB as `float64`; eight bands at once is several GB. This
document explains how the framework processes the full strips with **flat memory
use** (independent of raster height) in a **single pass**.

---

## The core idea: windowed streaming

Each band is processed in blocks of `window_lines` along-track (default
**2048**). For every block: read the DN window → calibrate → write the radiance
window → update QA — then release it. Memory is **O(window)**, not **O(band)**.

The loop lives in [`l1b_pipeline.py`](../src/spp/pipeline/l1b_pipeline.py)
(`_process_band`), using `rasterio.windows.Window(0, line_start, width, height)`
to read and write only the current block.

### Memory: full band vs window

| Buffer | Full band | Window (2048 lines) |
|---|---:|---:|
| DN (`uint16`) | ~240 MB | ~16 MB |
| Radiance (`float32`) | ~480 MB | ~32 MB |
| + intermediates (darkfield, etc.) | ~1 GB | tens of MB |
| **× 8 bands** | **several GB** | **flat, tens of MB** |

Peak RAM is **constant regardless of strip height** — a 300 000-line strip would
process in the same memory footprint.

---

## Supporting optimisations

| # | Technique | What it does | Where |
|---|---|---|---|
| 1 | **Lazy bands** | `Band` holds only path + raster properties, never pixels. The reader opens TIFFs only to read metadata, so assembling an 8-band acquisition costs ~no memory. | [`core/band.py`](../src/spp/core/band.py) |
| 2 | **`float32` output** | Halves memory and I/O vs `float64`; verified negligible precision loss (2.3 × 10⁻⁵). | [`l1b_calibrator.py`](../src/spp/calibration/l1b_calibrator.py) |
| 3 | **Cached per-band context** | Coefficient arrays (cast to `float32`) and the per-line temperature profile are computed **once per band** and reused across all windows, not per window. | `l1b_calibrator.py` (`_build_context`) |
| 4 | **Vectorised broadcasting** | The whole window is calibrated in one NumPy expression — no Python per-pixel loops: `darkfield = ti[None,:] + tg[None,:] * T[:,None]` → `(h, 4096)`. | `l1b_calibrator.py` (`calibrate`) |
| 5 | **Single-pass QA** | QA statistics accumulate **in the same window loop**, so quality metrics need no second read of the multi-GB output. | [`radiometric_validator.py`](../src/spp/qa/radiometric_validator.py) |
| 6 | **Block-aligned, compressed output** | Output GeoTIFF is tiled 512×512 (matching the input COG — Cloud-Optimized GeoTIFF), DEFLATE + floating-point predictor, `BIGTIFF=IF_SAFER` for >4 GB safety. Window writes align with the tile grid. | [`geotiff_writer.py`](../src/spp/products/geotiff_writer.py) |
| 7 | **Overviews for previews** | Internal overviews are built once on close; the quicklook reads **decimated** (`out_shape` + averaging) via those overviews, never full resolution. | [`quicklook.py`](../src/spp/qa/quicklook.py) |

---

## Window / tile alignment

The processing window and the GeoTIFF tiles are **two different units** that are
made compatible on purpose:

- a **tile** (512 × 512) is the *physical storage* block inside the GeoTIFF — the
  unit of compression and random access on disk;
- a **window** (`window_lines` tall, full width) is the *logical processing*
  chunk read/calibrated/written per iteration.

`window_lines = 2048` is chosen as an exact multiple of the tile size, and the
window spans the full width, so each window covers a whole number of tiles:

```text
width : 4096 / 512 = 8 tiles
height: 2048 / 512 = 4 tiles
        → one window = 8 × 4 = 32 complete tiles
```

Because every window lands on tile boundaries:

- **reads** decode only complete block-rows of the input COG (no partial-tile
  waste);
- **writes** cover complete output tiles, so each tile is compressed **once** —
  a non-aligned window (e.g. `1000`) would split a 512-line tile across two
  writes, forcing a read-modify-write that re-reads and re-compresses the
  half-filled tile.

The strip height (30948) is not a multiple of either, but `15 × 2048 = 30720 =
60 × 512`, so the first 30720 lines are perfectly aligned and the final partial
window (228 lines) coincides exactly with the final partial tile row — GDAL pads
internally. Alignment holds end to end.

## The `window_lines` tuning knob

Because `line_start` is threaded through every block, the result is **identical
regardless of window size** (verified — see
[`validation_strategy.md`](validation_strategy.md), "Window invariance"). So
`window_lines` is a free knob trading memory against per-call overhead:

- **larger** window → fewer I/O calls, more RAM per step;
- **smaller** window → minimal RAM, slightly more overhead.

Keep it a **multiple of the 512 tile size** (512, 2048, 4096 ...) so windows stay
tile-aligned (see above). Output is byte-identical either way, so it can be tuned
to the host without any risk to correctness.

---

## Reference runtime

A full 8-band strip (~4096 × 30948) calibrates in **~2 minutes** on a laptop
(≈13–18 s per band; ~123 s for the 8-band calibration, ~127 s end-to-end
including the quicklook and STAC item, on the reference scene — first, cold-cache
run), producing ~4 GB of `float32` output, with flat memory.

---

## Not yet optimised: per-band parallelism

Bands are currently processed **sequentially**, but they are **fully
independent** — there is no cross-band dependency in L1B. Processing them
concurrently (multiprocessing or Dask, one band per worker) would cut wall-clock
time roughly in proportion to the number of cores. This is the single largest
remaining performance opportunity and is noted in
[`limitations.md`](limitations.md). The current design already supports it: the
per-band work is self-contained in `_process_band`.
