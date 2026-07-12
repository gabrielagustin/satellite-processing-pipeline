"""Absolute geolocation: correcting the pointing against an independent reference.

Everything else in this level is *relative*. The bands are co-registered onto each
other, the footprint is compared against a catalogue geometry the provider derived from
the same telemetry we consume — none of it can see a bias shared by the whole product.
And this acquisition has **no satellite-navigation lock**: its positions are propagated,
not measured, so a shared bias is exactly what one expects.

On the reference acquisition that bias is **about a kilometre**, almost entirely
along-track. It was invisible to every check the framework could run, and it was found by
a user overlaying the product on a basemap.

The reference is the terrain
----------------------------
An external orthoimage is the textbook answer, and it means downloading imagery, matching
seasons, and worrying about clouds. But the pipeline **already downloads a georeferenced
reference**: the elevation model. Its coastline is ground truth, it needs no season, it
has no clouds, and it is offline once cached.

Water is near-black in the near infrared and land is not, so the image's own land/water
boundary can be matched against the terrain's. That measures the product's absolute
displacement directly.

Where the correction goes
-------------------------
Into the **boresight** — the instrument's pointing — not into a shift of the output
image. A rigid image shift corrects the average and leaves the terrain-dependent part
wrong, and it leaves the sensor model knowingly incorrect for everything computed from it
afterwards. The correction belongs where the error physically originates.

A degeneracy that must be stated, not hidden
--------------------------------------------
From a single strip, an along-track ground shift is produced *identically* by a **pitch
bias** and by a **clock offset** (`Δt = dy / v_ground`). They are not separable from this
data. The convention here is to absorb it into pitch and leave the delivered clock offset
alone — the product is correct either way, and the attribution is a convention. Saying so
is the difference between a measurement and a claim.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np

from spp.geometry.camera import Camera
from spp.geometry.sensor_model import SensorModel
from spp.geometry.terrain import DEMTerrain
from spp.refine.matcher import has_texture, phase_correlate

logger = logging.getLogger(__name__)

PERTURBATION_RAD = 1e-5
"""Step used to probe how the ground point responds to a boresight rotation."""

LAND_HEIGHT_M = 2.0
"""Orthometric height above which the terrain is called land."""


@dataclass(frozen=True)
class AbsoluteCorrection:
    """The pointing correction estimated against the terrain.

    Attributes
    ----------
    roll_rad, pitch_rad:
        Boresight rotations to apply, radians. ``pitch`` carries the along-track
        correction and is **degenerate with a clock offset** — see the module docstring.
    offset_before_m, offset_after_m:
        Ground displacement between the product and the terrain reference, before and
        after. The pair is the evidence.
    equivalent_clock_offset_s:
        The clock error that would produce the same along-track shift. Reported so the
        degeneracy is visible rather than implied.
    peak:
        Correlation peak of the match, in ``[0, 1]``.
    fitted:
        Whether a correction was accepted.
    reason:
        Why not, when ``fitted`` is ``False``.
    """

    roll_rad: float
    pitch_rad: float
    offset_before_m: float
    offset_after_m: float
    equivalent_clock_offset_s: float
    peak: float
    fitted: bool
    reason: str = ""


def _land_masks(
    model: SensorModel,
    band: str,
    image: np.ndarray,
    terrain: DEMTerrain,
    *,
    step: int,
    land_percentile: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """The image's land mask and the terrain's, sampled on the same sensor lattice.

    Sampling on the *sensor* lattice rather than the map grid means no warp is needed:
    the model already says where each sample landed, and the terrain can be read there.
    """
    n_lines, n_columns = image.shape
    lines = np.arange(0, n_lines, step, dtype=np.float64)
    columns = np.arange(0, n_columns, step, dtype=np.float64)
    mesh_l, mesh_c = np.meshgrid(lines, columns, indexing="ij")

    radiance = image[mesh_l.astype(int), mesh_c.astype(int)].astype(np.float64)
    valid = np.isfinite(radiance) & (radiance > 0)
    if valid.sum() < 100:
        raise ValueError("Too few valid samples to match against the terrain")

    threshold = np.percentile(radiance[valid], land_percentile)
    image_land = np.where(valid, radiance > threshold, 0.0).astype(np.float64)

    ground = model.locate(band, mesh_l.ravel(), mesh_c.ravel(), terrain=terrain)
    lon = ground[:, 0].reshape(mesh_l.shape)
    lat = ground[:, 1].reshape(mesh_l.shape)

    orthometric = terrain.height(lon, lat) - terrain.geoid(lon, lat)
    terrain_land = (orthometric > LAND_HEIGHT_M).astype(np.float64)

    return image_land, terrain_land, ground, float(np.nanmean(lat))


def estimate(
    model: SensorModel,
    band: str,
    image: np.ndarray,
    terrain: DEMTerrain,
    *,
    ground_speed_m_s: float,
    step: int = 16,
    land_percentile: float = 55.0,
    min_peak: float = 0.005,
) -> AbsoluteCorrection:
    """Measure the product's absolute displacement, and turn it into a pointing correction.

    Parameters
    ----------
    model:
        The sensor model, with whatever calibration is already in force.
    band:
        Band whose imagery is supplied. **Near-infrared if available** — the land/water
        contrast is what carries the signal.
    image:
        The band's raster, in sensor coordinates.
    terrain:
        Elevation model, which supplies both the surface and the reference coastline.
    ground_speed_m_s:
        Speed of the sub-satellite point, used to report the clock offset that would
        explain the same along-track shift.
    step:
        Lattice spacing for the match, in detector samples.
    min_peak:
        Reject a match weaker than this. A flat correlation surface means the scene had
        no coastline to lock onto, and an unlocked match is not a measurement.

    Returns
    -------
    AbsoluteCorrection
        With ``fitted=False`` and a reason when the scene cannot support the estimate.
    """
    image_land, terrain_land, ground, lat0 = _land_masks(
        model, band, image, terrain, step=step, land_percentile=land_percentile
    )

    land_fraction = float(terrain_land.mean())
    if not 0.05 < land_fraction < 0.95:
        return _refused(
            f"the footprint is {land_fraction * 100:.0f}% land; there is no coastline to "
            "match against"
        )

    # The image's land mask against the terrain's, in sensor-lattice cells.
    d_line, d_column, peak = phase_correlate(terrain_land, image_land)
    if peak < min_peak:
        return _refused(f"the terrain match did not lock on (peak {peak:.4f})")

    # Convert the lattice displacement into a ground displacement, using the model's own
    # local scale rather than assuming one.
    centre = np.array([image_land.shape[0] // 2, image_land.shape[1] // 2])
    offset = _lattice_to_ground(
        model, band, terrain, centre, np.array([d_line, d_column]), step=step, lat0=lat0
    )
    before = float(np.linalg.norm(offset))

    # How the ground point responds to a boresight rotation -- measured, not derived,
    # for the same reason as the line-of-sight fit: the chain of sign conventions from
    # camera to body to orbit to Earth is not worth re-deriving by hand.
    jacobian = np.empty((2, 2))
    base = _mean_ground(model, band, terrain, lat0)
    for index, axis in enumerate(("roll", "pitch")):
        perturbed = _with_boresight(model, **{axis: PERTURBATION_RAD})
        moved = _mean_ground(perturbed, band, terrain, lat0)
        jacobian[:, index] = (moved - base) / PERTURBATION_RAD

    if abs(np.linalg.det(jacobian)) < 1e-9:
        return _refused("the model is insensitive to a boresight rotation")

    roll, pitch = np.linalg.solve(jacobian, -offset)

    corrected = _with_boresight(model, roll=roll, pitch=pitch)
    after_offset = _mean_ground(corrected, band, terrain, lat0) - base - (-offset)
    after = float(np.linalg.norm(after_offset))

    along = float(offset[1])  # the strip runs north-south; the along-track component
    clock = along / ground_speed_m_s if ground_speed_m_s else float("nan")

    logger.info(
        "Absolute geolocation: the product sits %.0f m from the terrain reference "
        "(%.0f m east, %.0f m north). Correcting the boresight by roll %+.1f urad, "
        "pitch %+.1f urad. NOTE: the along-track part is equally explained by a clock "
        "offset of %+.3f s -- the two are not separable from one strip, and pitch is a "
        "convention.",
        before,
        offset[0],
        offset[1],
        roll * 1e6,
        pitch * 1e6,
        clock,
    )
    return AbsoluteCorrection(
        roll_rad=float(roll),
        pitch_rad=float(pitch),
        offset_before_m=before,
        offset_after_m=after,
        equivalent_clock_offset_s=float(clock),
        peak=float(peak),
        fitted=True,
    )


def _mean_ground(model: SensorModel, band: str, terrain, lat0: float) -> np.ndarray:
    """Where the model puts the raster's centre, in local metres."""
    timing = model.timing[band]
    line = np.array([timing.n_lines / 2.0])
    column = np.array([model.camera.intrinsics.sensor_size_px[0] / 2.0])
    ground = model.locate(band, line, column, terrain=terrain)
    return _to_metres(ground[0, :2], lat0)


def _lattice_to_ground(
    model: SensorModel,
    band: str,
    terrain,
    centre: np.ndarray,
    shift: np.ndarray,
    *,
    step: int,
    lat0: float,
) -> np.ndarray:
    """A displacement in lattice cells, expressed as a ground displacement in metres."""
    line = centre[0] * step
    column = centre[1] * step
    here = model.locate(band, np.array([line]), np.array([column]), terrain=terrain)
    there = model.locate(
        band,
        np.array([line + shift[0] * step]),
        np.array([column + shift[1] * step]),
        terrain=terrain,
    )
    return _to_metres(there[0, :2], lat0) - _to_metres(here[0, :2], lat0)


def _to_metres(lon_lat: np.ndarray, lat0: float) -> np.ndarray:
    return np.array(
        [
            lon_lat[0] * 111_320.0 * np.cos(np.radians(lat0)),
            lon_lat[1] * 110_574.0,
        ]
    )


def _with_boresight(model: SensorModel, *, roll: float = 0.0, pitch: float = 0.0) -> SensorModel:
    """A copy of the model whose boresight carries an extra rotation."""
    cos_r, sin_r = np.cos(roll), np.sin(roll)
    cos_p, sin_p = np.cos(pitch), np.sin(pitch)
    rotation = np.array(
        [[1.0, 0.0, 0.0], [0.0, cos_r, -sin_r], [0.0, sin_r, cos_r]]
    ) @ np.array([[cos_p, 0.0, sin_p], [0.0, 1.0, 0.0], [-sin_p, 0.0, cos_p]])

    camera = Camera(
        model.camera.intrinsics,
        model.camera.bands,
        boresight=rotation @ model.camera.boresight,
    )
    return SensorModel(
        camera,
        model.ephemeris,
        model.attitude,
        model.timing,
        convention=model.convention,
        rate_aware=model.rate_aware,
    )


def _refused(reason: str) -> AbsoluteCorrection:
    """No estimate, and why. Never a silent zero -- that would read as 'no error'."""
    logger.warning("Absolute geolocation: not corrected — %s", reason)
    return AbsoluteCorrection(
        roll_rad=0.0,
        pitch_rad=0.0,
        offset_before_m=float("nan"),
        offset_after_m=float("nan"),
        equivalent_clock_offset_s=float("nan"),
        peak=0.0,
        fitted=False,
        reason=reason,
    )


# -- against a reference orthoimage ------------------------------------------


def estimate_against_reference(
    model: SensorModel,
    band: str,
    image: np.ndarray,
    reference_path: str,
    terrain,
    *,
    ground_speed_m_s: float,
    step: int = 16,
    window: int = 64,
    min_matches: int = 10,
) -> AbsoluteCorrection:
    """Measure absolute displacement against an independently georeferenced image.

    Sharper than the terrain's coastline by an order of magnitude, because it matches
    *texture* rather than a binary land mask: a coastline pins you to how well a shoreline
    is defined, an orthoimage pins you to a pixel.

    The reference is sampled **at the ground points the model predicts**, on the sensor's
    own lattice. No warp is needed, and no resampling of our data is done before the
    measurement — which matters, because resampling is precisely what would smear the
    displacement being measured.

    Parameters
    ----------
    model:
        The sensor model, with whatever calibration is in force.
    band:
        Band whose imagery is supplied.
    image:
        The band's raster, in sensor coordinates.
    reference_path:
        A georeferenced raster of the same ground (see :mod:`spp.refine.reference`).
    terrain:
        Elevation model, so the ground points are the orthorectified ones.
    ground_speed_m_s:
        Speed of the sub-satellite point, to report the equivalent clock offset.
    step:
        Lattice spacing in detector samples.
    window:
        Correlation window, in lattice cells.
    min_matches:
        Refuse to estimate on fewer textured windows than this.
    """
    import rasterio
    from rasterio import Affine
    from rasterio.crs import CRS
    from rasterio.warp import transform as warp_transform
    from rasterio.windows import from_bounds

    n_lines, n_columns = image.shape
    lines = np.arange(0, n_lines, step, dtype=np.float64)
    columns = np.arange(0, n_columns, step, dtype=np.float64)
    mesh_l, mesh_c = np.meshgrid(lines, columns, indexing="ij")

    ours = image[mesh_l.astype(int), mesh_c.astype(int)].astype(np.float64)

    ground = model.locate(band, mesh_l.ravel(), mesh_c.ravel(), terrain=terrain)
    lon = ground[:, 0]
    lat = ground[:, 1]
    lat0 = float(np.nanmean(lat))

    # Sample the reference where the model says each detector sample landed.
    #
    # Only the footprint's window is read, and only as finely as the lattice can use.
    # Reading the whole mosaic at native resolution would pull hundreds of megabytes
    # across the network to answer a question posed on a grid an order of magnitude
    # coarser -- which is how this first ran, and why it never finished.
    with rasterio.open(reference_path) as src:
        xs, ys = warp_transform(CRS.from_epsg(4326), src.crs, lon.tolist(), lat.tolist())
        xs = np.asarray(xs)
        ys = np.asarray(ys)

        pad = 4.0 * abs(src.transform.a)
        read_window = from_bounds(
            float(np.nanmin(xs)) - pad,
            float(np.nanmin(ys)) - pad,
            float(np.nanmax(xs)) + pad,
            float(np.nanmax(ys)) + pad,
            transform=src.transform,
        ).round_offsets().round_lengths()

        # At NATIVE resolution. An earlier version decimated the reference to save
        # bandwidth, then sampled it with nearest-neighbour -- two quantisations stacked,
        # which put a floor of tens of metres under a measurement whose whole purpose is
        # to resolve tens of metres. The saving was imaginary anyway: only the footprint's
        # window is read, and at 10 m that is some tens of megabytes.
        patch = src.read(1, window=read_window, boundless=True, fill_value=0).astype(
            np.float64
        )
        patch_transform = src.window_transform(read_window)
        nodata = src.nodata

    if nodata is not None:
        patch[patch == nodata] = np.nan
    patch[patch == 0] = np.nan  # reference tiles pad with zero outside their footprint

    reference_flat = _bilinear_sample(patch, patch_transform, xs, ys)
    theirs = reference_flat.reshape(mesh_l.shape)
    if np.isfinite(theirs).mean() < 0.2:
        return _refused("the reference orthoimage does not cover the footprint")

    # Match window by window, so the result is a robust consensus rather than one
    # correlation over a scene that is mostly featureless water.
    # Overlapping windows, so a narrow strip still yields enough cross-track samples.
    # With non-overlapping windows the swath gives only a handful across, and the
    # cross-track estimate is then decided by four numbers.
    stride = max(1, window // 2)
    shifts = []
    for row in range(0, ours.shape[0] - window + 1, stride):
        for col in range(0, ours.shape[1] - window + 1, stride):
            a = ours[row : row + window, col : col + window]
            b = theirs[row : row + window, col : col + window]
            if not np.isfinite(b).all() or not has_texture(a) or not has_texture(b):
                continue
            d_line, d_column, peak = phase_correlate(b, a)
            if peak < 0.02:
                continue
            shifts.append((d_line, d_column, peak))

    if len(shifts) < min_matches:
        return _refused(
            f"only {len(shifts)} textured windows matched the reference (need "
            f"{min_matches}); the scene may be mostly water or cloud"
        )

    shift_array = np.array([[s[0], s[1]] for s in shifts])
    peaks = np.array([s[2] for s in shifts])

    # Reject outliers, then weight what survives by how firmly it locked on. A window
    # over a cloud edge, a ship or a moving shadow will mismatch confidently, and a plain
    # median gives it the same vote as a window that nailed a coastline.
    centre = np.median(shift_array, axis=0)
    deviation = np.abs(shift_array - centre)
    mad = np.median(deviation, axis=0) * 1.4826
    inliers = np.all(deviation <= np.maximum(3.0 * mad, 0.5), axis=1)
    if inliers.sum() >= min_matches:
        shift_array = shift_array[inliers]
        peaks = peaks[inliers]

    weights = peaks / peaks.sum()
    consensus = (shift_array * weights[:, None]).sum(axis=0)
    scatter = float(np.median(np.abs(shift_array - consensus)) * 1.4826)
    peak = float(np.median(peaks))

    centre = np.array([ours.shape[0] // 2, ours.shape[1] // 2])
    offset = _lattice_to_ground(model, band, terrain, centre, consensus, step=step, lat0=lat0)
    before = float(np.linalg.norm(offset))

    jacobian = np.empty((2, 2))
    base = _mean_ground(model, band, terrain, lat0)
    for index, axis in enumerate(("roll", "pitch")):
        perturbed = _with_boresight(model, **{axis: PERTURBATION_RAD})
        jacobian[:, index] = (
            _mean_ground(perturbed, band, terrain, lat0) - base
        ) / PERTURBATION_RAD

    if abs(np.linalg.det(jacobian)) < 1e-9:
        return _refused("the model is insensitive to a boresight rotation")

    roll, pitch = np.linalg.solve(jacobian, -offset)
    clock = float(offset[1]) / ground_speed_m_s if ground_speed_m_s else float("nan")

    logger.info(
        "Absolute geolocation vs the reference orthoimage: the product sits %.0f m away "
        "(%.0f m east, %.0f m north) over %d matched windows, scatter %.1f lattice cells. "
        "Correcting the boresight by roll %+.1f urad, pitch %+.1f urad. NOTE: the "
        "along-track part is equally explained by a clock offset of %+.3f s; the two are "
        "not separable from one strip, and pitch is a convention.",
        before,
        offset[0],
        offset[1],
        len(shifts),
        scatter,
        roll * 1e6,
        pitch * 1e6,
        clock,
    )
    return AbsoluteCorrection(
        roll_rad=float(roll),
        pitch_rad=float(pitch),
        offset_before_m=before,
        offset_after_m=float("nan"),
        equivalent_clock_offset_s=clock,
        peak=peak,
        fitted=True,
    )


def _bilinear_sample(
    grid: np.ndarray, transform, xs: np.ndarray, ys: np.ndarray
) -> np.ndarray:
    """Sample a georeferenced grid at arbitrary coordinates, bilinearly.

    Nearest-neighbour would quantise the reference to its own pixel — which puts a floor
    under a measurement whose entire purpose is to resolve displacements *smaller* than
    that pixel.
    """
    cols_f, rows_f = (~transform) * (xs, ys)

    rows = np.clip(np.floor(rows_f).astype(np.int64), 0, grid.shape[0] - 2)
    cols = np.clip(np.floor(cols_f).astype(np.int64), 0, grid.shape[1] - 2)
    dr = rows_f - rows
    dc = cols_f - cols

    top = grid[rows, cols] * (1 - dc) + grid[rows, cols + 1] * dc
    bottom = grid[rows + 1, cols] * (1 - dc) + grid[rows + 1, cols + 1] * dc
    values = top * (1 - dr) + bottom * dr

    outside = (
        (rows_f < 0)
        | (rows_f >= grid.shape[0] - 1)
        | (cols_f < 0)
        | (cols_f >= grid.shape[1] - 1)
    )
    values[outside] = np.nan
    return values


def refine_against_reference(
    model: SensorModel,
    band: str,
    image: np.ndarray,
    reference_path: str,
    terrain,
    *,
    ground_speed_m_s: float,
    steps: tuple[int, ...] = (32, 16, 8),
    tolerance_m: float = 5.0,
) -> tuple[SensorModel, list[AbsoluteCorrection]]:
    """Estimate and correct the pointing bias **coarse to fine**.

    A single pass cannot do this, and the reason is worth stating because it is the
    opposite of the intuition. Phase correlation finds a displacement by overlapping two
    windows, so it can only see a shift that is *small compared to the window*. Make the
    lattice finer and every cell covers less ground — which means a fixed displacement
    spans **more** cells, until it approaches the window size and the correlation collapses.

    On the reference acquisition the bias is ~919 m. On a 60 m lattice that is 15 cells in
    a 64-cell window: recoverable. On a 30 m lattice it is 30 cells — half the window —
    and the estimator returns a confident 25 m, which is noise. Refining the lattice made
    the measurement *worse*, and it did so silently.

    So the bias is removed in stages: a coarse lattice sees the large shift, the model is
    corrected, and only then does a finer lattice have a small shift left to resolve. Each
    stage is reported, and the sequence stops when the residual stops moving.

    Returns
    -------
    tuple
        The corrected model, and the correction from each stage.
    """
    corrections: list[AbsoluteCorrection] = []
    current = model

    for stage, step in enumerate(steps, start=1):
        correction = estimate_against_reference(
            current,
            band,
            image,
            reference_path,
            terrain,
            ground_speed_m_s=ground_speed_m_s,
            step=step,
        )
        if not correction.fitted:
            logger.warning(
                "Absolute geolocation, stage %d (lattice step %d): %s",
                stage,
                step,
                correction.reason,
            )
            break

        corrections.append(correction)
        current = _with_boresight(
            current, roll=correction.roll_rad, pitch=correction.pitch_rad
        )

        if correction.offset_before_m < tolerance_m:
            logger.info(
                "Absolute geolocation converged at stage %d: %.0f m remaining.",
                stage,
                correction.offset_before_m,
            )
            break

    return current, corrections
