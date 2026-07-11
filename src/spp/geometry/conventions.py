"""Resolving the telemetry's undocumented conventions, by measurement.

Four properties of the platform telemetry are not specified by the acquisition
package:

* which frame the attitude quaternion relates the body frame to,
* in which direction it rotates,
* whether its components are stored scalar-first or scalar-last,
* whether raster row zero is the first line acquired or the last.

There are a few dozen ways to combine them and exactly one is right. Reading them
wrong does not produce a subtly degraded product — it produces a footprint on the
wrong continent, which is the saving grace: the error is so large that a single
cheap measurement separates the candidates beyond any doubt.

The method
----------
Project the **four corners** of the raster to the ground under each candidate
convention — a few dozen rays in total, milliseconds of work — and score each by
how far the resulting footprint lands from a footprint that is already known (the
one the provider delivered in the product's catalogue entry).

The correct convention lands within the model's own error, which for an
uncalibrated model is kilometres. Every incorrect one is off by tens to thousands
of kilometres. The margin between best and second-best is the evidence: a narrow
margin means the test did not discriminate and must not be trusted, so it is
reported alongside the winner rather than hidden.

This is a *resolution* step, not a fit. It chooses among discrete readings of the
data; it does not tune any continuous parameter, and it cannot improve a model
that is otherwise wrong.
"""

from __future__ import annotations

import itertools
import logging
from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np

from spp.geometry.camera import Camera
from spp.geometry.ephemeris import Attitude, Ephemeris
from spp.geometry.sensor_model import (
    AttitudeFrame,
    AxisSign,
    ColumnAxis,
    Convention,
    EphemerisFrame,
    QuaternionDirection,
    QuaternionOrder,
    SensorModel,
)
from spp.geometry.timing import LineTiming

EPHEMERIS_FRAMES: tuple[EphemerisFrame, ...] = ("eci", "ecef")
ATTITUDE_FRAMES: tuple[AttitudeFrame, ...] = ("lvlh", "eci", "ecef")
QUATERNION_ORDERS: tuple[QuaternionOrder, ...] = ("scalar_first", "scalar_last")
QUATERNION_DIRECTIONS: tuple[QuaternionDirection, ...] = ("body_to_ref", "ref_to_body")
SCAN_DIRECTIONS: tuple[int, ...] = (1, -1)
COLUMN_AXES: tuple[ColumnAxis, ...] = ("x", "y")
AXIS_SIGNS: tuple[AxisSign, ...] = (1, -1)

logger = logging.getLogger(__name__)

MEAN_EARTH_RADIUS_M = 6371008.8
"""Mean Earth radius, for converting angular separations to ground distance."""


@dataclass(frozen=True)
class ConventionScore:
    """How well one candidate convention explains the data.

    Two independent probes, because neither alone resolves every axis.

    Attributes
    ----------
    convention:
        The candidate.
    corner_rmse_m:
        Root-mean-square distance between each modelled corner and the nearest
        corner of the reference footprint, in metres. Resolves the *absolute*
        axes — which frames the ephemeris and attitude live in — because getting
        those wrong throws the footprint hundreds or thousands of kilometres.
    band_coherence_m:
        Mean ground distance between the same raster sample located through two
        bands with widely separated detector rows. Resolves the *scan direction*,
        which the footprint cannot see: reversing the line order maps the four
        corners onto each other, leaving the footprint identical. Band coherence
        is not fooled — the per-band time stagger only cancels the per-band view
        angle when the scan runs the way the model thinks it does. Flip it and
        two bands land kilometres apart.
    centroid_offset_m:
        Distance between the modelled and reference footprint centroids, metres.
        Distinguishes a *displaced* footprint (a pointing or timing error, which
        refinement can fix) from a *malformed* one (a wrong convention, which it
        cannot).
    valid:
        Whether every ray actually intersected the Earth. A convention that
        points the sensor at the sky scores no distance at all.
    """

    convention: Convention
    corner_rmse_m: float
    band_coherence_m: float
    centroid_offset_m: float
    valid: bool

    def plausible(self, *, gross_error_m: float = 100_000.0) -> bool:
        """Whether the footprint is close enough to be worth ranking at all."""
        return self.valid and self.corner_rmse_m <= gross_error_m

    def rank_key(self, *, gross_error_m: float = 100_000.0) -> tuple[bool, bool, float]:
        """Sort key: reject the implausible, then rank the rest by coherence."""
        return (not self.valid, not self.plausible(gross_error_m=gross_error_m),
                self.band_coherence_m if self.plausible(gross_error_m=gross_error_m)
                else self.corner_rmse_m)

    def __str__(self) -> str:
        if not self.valid:
            return f"{self.convention.describe()}: rays miss the Earth"
        return (
            f"{self.convention.describe()}: "
            f"corner RMSE {self.corner_rmse_m / 1000:,.1f} km, "
            f"band coherence {self.band_coherence_m:,.0f} m"
        )


GROSS_ERROR_M = 100_000.0
"""Footprint error above which a candidate is simply wrong, in metres.

The two probes must **not** be added together, and why is the crux of the design.

The footprint probe is contaminated by an error no convention can remove: without
a GNSS (Global Navigation Satellite System) lock the platform's absolute position
is uncertain and its pointing bias is unknown, so *even the correct convention*
misses the delivered footprint by tens of kilometres. Summing the probes would let
that irreducible bias swamp the hundreds-of-metres signal that actually separates
the fine conventions — the blunt instrument would overrule the sharp one.

Band coherence has no such contamination. An absolute bias displaces both bands
identically and **cancels** in their difference, so coherence measures the model's
internal consistency and is blind to precisely the error the footprint cannot see
past. It is the sharp instrument.

So each probe does what it is good for: the footprint rejects gross failures — a
wrong reference frame throws it 700 km or more — and coherence ranks the
survivors. This threshold is the line between those regimes: far above the tens of
kilometres a *correct* model can be off by, far below the hundreds a *wrong* frame
produces. Nothing observed lands in between, which is what makes the split safe.
"""


def candidate_conventions(
    *,
    ephemeris_frames: Iterable[EphemerisFrame] = EPHEMERIS_FRAMES,
    attitude_frames: Iterable[AttitudeFrame] = ATTITUDE_FRAMES,
    quaternion_orders: Iterable[QuaternionOrder] = QUATERNION_ORDERS,
    quaternion_directions: Iterable[QuaternionDirection] = QUATERNION_DIRECTIONS,
    scan_directions: Iterable[int] = SCAN_DIRECTIONS,
    column_axes: Iterable[ColumnAxis] = COLUMN_AXES,
    column_signs: Iterable[AxisSign] = AXIS_SIGNS,
    row_signs: Iterable[AxisSign] = AXIS_SIGNS,
) -> list[Convention]:
    """Every combination of the ambiguous conventions."""
    return [
        Convention(
            ephemeris_frame=eph,
            attitude_frame=frame,
            quaternion_order=order,
            quaternion_direction=direction,
            scan_direction=scan,
            column_axis=axis,
            column_sign=col_sign,
            row_sign=row_sign,
        )
        for eph, frame, order, direction, scan, axis, col_sign, row_sign in itertools.product(
            ephemeris_frames,
            attitude_frames,
            quaternion_orders,
            quaternion_directions,
            scan_directions,
            column_axes,
            column_signs,
            row_signs,
        )
    ]


def footprint_corners(
    model: SensorModel, band: str, *, n_lines: int, n_columns: int
) -> np.ndarray:
    """Ground positions of the raster's four corners.

    Returns
    -------
    numpy.ndarray
        Shape ``(4, 2)``, ``(longitude_deg, latitude_deg)`` per corner, ordered
        first-line-first-column, first-line-last-column, last-line-last-column,
        last-line-first-column. ``NaN`` where a ray missed the Earth.
    """
    last_line = n_lines - 1
    last_col = n_columns - 1
    lines = np.array([0, 0, last_line, last_line], dtype=np.float64)
    columns = np.array([0, last_col, last_col, 0], dtype=np.float64)
    return model.locate(band, lines, columns)[:, :2]


def band_coherence(
    model: SensorModel,
    band_a: str,
    band_b: str,
    *,
    n_lines: int,
    n_columns: int,
    n_samples: int = 16,
) -> float:
    """Mean ground distance between the same raster sample seen through two bands.

    Two bands read from different detector rows look at different along-track
    angles, and the product staggers their line timestamps to compensate. The
    compensation only works in one temporal direction. So if the model runs the
    scan the right way, raster sample ``(line, column)`` lands on the *same*
    ground point in both bands, to within a few pixels; if it runs it backwards,
    the stagger is applied against itself and the two land kilometres apart.

    That makes this a probe of the scan direction which — unlike the footprint —
    cannot be fooled by the symmetry of a rectangle.

    Parameters
    ----------
    model:
        The sensor model under test.
    band_a, band_b:
        Bands to compare. Choose the two with the most widely separated detector
        rows: the further apart they look, the sharper the probe.
    n_lines, n_columns:
        Raster dimensions.
    n_samples:
        Number of samples along the strip.

    Returns
    -------
    float
        Mean distance, metres. ``inf`` if either band fails to locate.
    """
    # Stay clear of the strip ends, where one band's stagger runs off the raster.
    lines = np.linspace(0.1 * n_lines, 0.9 * n_lines, n_samples)
    columns = np.full_like(lines, 0.5 * n_columns)

    ground_a = model.locate(band_a, lines, columns)
    ground_b = model.locate(band_b, lines, columns)
    if not (np.all(np.isfinite(ground_a)) and np.all(np.isfinite(ground_b))):
        return float("inf")

    lat0 = float(np.mean(ground_a[:, 1]))
    xy_a = _local_metres(ground_a[:, :2], lat0)
    xy_b = _local_metres(ground_b[:, :2], lat0)
    return float(np.mean(np.linalg.norm(xy_a - xy_b, axis=-1)))


def score_convention(
    corners: np.ndarray, reference_corners: np.ndarray
) -> tuple[float, float, bool]:
    """Compare a modelled footprint against a reference one.

    The modelled corners are matched to the *nearest* reference corner rather
    than pairwise in order, because the two footprints need not enumerate their
    corners from the same starting vertex or in the same rotational sense. A
    correct convention lands every corner near a distinct reference corner; a
    wrong one lands them all far from every corner, so the nearest-corner
    relaxation cannot rescue it.

    Returns
    -------
    tuple
        ``(corner_rmse_m, centroid_offset_m, valid)``.
    """
    if not np.all(np.isfinite(corners)):
        return float("inf"), float("inf"), False

    lat0 = float(np.mean(reference_corners[:, 1]))
    model_xy = _local_metres(corners, lat0)
    ref_xy = _local_metres(reference_corners, lat0)

    distances = np.linalg.norm(model_xy[:, None, :] - ref_xy[None, :, :], axis=-1)
    nearest = distances.min(axis=1)
    rmse = float(np.sqrt(np.mean(nearest**2)))
    centroid = float(np.linalg.norm(model_xy.mean(axis=0) - ref_xy.mean(axis=0)))
    return rmse, centroid, True


def resolve(
    camera: Camera,
    ephemeris: Ephemeris,
    attitude: Attitude,
    timing: dict[str, LineTiming],
    *,
    band: str,
    n_lines: int,
    n_columns: int,
    reference_footprint: Sequence[Sequence[float]],
    coherence_bands: tuple[str, str] | None = None,
    candidates: Sequence[Convention] | None = None,
) -> list[ConventionScore]:
    """Rank every candidate convention against a known footprint.

    Parameters
    ----------
    camera, ephemeris, attitude, timing:
        The sensor model's components.
    band:
        Band to project for the footprint probe. Any band works; the choice
        shifts the footprint by the band's own view angle, which is kilometres at
        most — far below the margin that separates a right convention from a
        wrong one.
    n_lines, n_columns:
        Raster dimensions of that band.
    reference_footprint:
        Known ground footprint as ``(longitude, latitude)`` vertices — for
        instance the polygon in the product's catalogue entry. A closing vertex
        equal to the first is tolerated and ignored.
    coherence_bands:
        The two bands used for the band-coherence probe. Defaults to the pair
        with the most widely separated detector rows, which is the sharpest.
    candidates:
        Conventions to test. Defaults to every combination.

    Returns
    -------
    list of ConventionScore
        Sorted best-first by :attr:`ConventionScore.total_m`. Callers should
        check :func:`is_decisive` before trusting the winner — a convention that
        merely wins is not the same as one the data actually resolves.
    """
    reference = _dedupe_closing_vertex(np.asarray(reference_footprint, dtype=np.float64))
    band_a, band_b = coherence_bands or _widest_band_pair(camera)

    scores: list[ConventionScore] = []
    for convention in candidates or candidate_conventions():
        model = SensorModel(
            camera, ephemeris, attitude, timing, convention=convention
        )
        try:
            corners = footprint_corners(
                model, band, n_lines=n_lines, n_columns=n_columns
            )
            coherence = band_coherence(
                model, band_a, band_b, n_lines=n_lines, n_columns=n_columns
            )
        except (ValueError, KeyError):
            # A degenerate state or an unknown band is a failed candidate, not a
            # failed run: the point of the harness is to let bad readings lose.
            scores.append(
                ConventionScore(convention, float("inf"), float("inf"), float("inf"), False)
            )
            continue

        rmse, centroid, valid = score_convention(corners, reference)
        scores.append(ConventionScore(convention, rmse, coherence, centroid, valid))

    scores.sort(key=lambda s: s.rank_key())
    return scores


def _widest_band_pair(camera: Camera) -> tuple[str, str]:
    """The two bands whose detector rows are furthest apart.

    They see the ground from the most different angles and are staggered furthest
    apart in time, so they are the most sensitive pair for the coherence probe.
    """
    by_row = sorted(camera.bands.values(), key=lambda b: b.effective_row)
    if len(by_row) < 2:
        raise ValueError(
            "Band coherence needs at least two bands; the camera has "
            f"{len(by_row)}"
        )
    return by_row[0].name, by_row[-1].name


AXES = (
    "ephemeris_frame",
    "attitude_frame",
    "quaternion_order",
    "quaternion_direction",
    "scan_direction",
    "column_axis",
    "column_sign",
    "row_sign",
)


def resolved_axes(
    scores: Sequence[ConventionScore], *, margin: float = 10.0
) -> dict[str, bool]:
    """Which individual conventions the data actually settles.

    A single global verdict is too blunt. The winner of the ranking has a value
    for every axis, but the data may pin some of them decisively and be nearly
    blind to others — and a value the data cannot see is a guess wearing the
    winner's clothes. This asks the question separately for each axis: **holding
    the winner's other choices fixed, how much worse is the best alternative on
    this one axis?**

    An axis is resolved when flipping it costs at least ``margin`` times the
    winner's error. An unresolved axis is not a failure of the harness; it is a
    fact about the data, and one the caller must carry forward — either by
    finding a sharper probe or by recording the choice as an assumption.

    Parameters
    ----------
    scores:
        Ranked output of :func:`resolve`, over the full candidate set.
    margin:
        Required ratio between the best alternative's error and the winner's.

    Returns
    -------
    dict
        Axis name to whether the data resolves it.
    """
    if not any(s.valid for s in scores):
        return dict.fromkeys(AXES, False)

    best = scores[0]
    baseline = max(best.band_coherence_m, 1e-9)

    # Every candidate, including the ones whose rays never reach the Earth. Those
    # are not unusable — they are the most decisively refuted of all, and dropping
    # them would silently discard the strongest evidence the harness has.
    by_convention = {s.convention: s for s in scores}

    def cost(score: ConventionScore) -> float:
        """What flipping to this candidate costs, on the sharper of the two probes.

        A candidate that misses the Earth, or that the footprint rules out, is
        infinitely costly. Among the plausible ones the cost is the coherence,
        which — unlike the footprint — is not blunted by the absolute bias.
        """
        return score.band_coherence_m if score.plausible() else float("inf")

    resolved: dict[str, bool] = {}
    for axis in AXES:
        alternatives = [
            cost(score)
            for convention, score in by_convention.items()
            if getattr(convention, axis) != getattr(best.convention, axis)
            and all(
                getattr(convention, other) == getattr(best.convention, other)
                for other in AXES
                if other != axis
            )
        ]
        resolved[axis] = bool(alternatives) and min(alternatives) / baseline >= margin
    return resolved


def is_decisive(scores: Sequence[ConventionScore], *, margin: float = 10.0) -> bool:
    """Whether the data resolves **every** ambiguous convention.

    The method rests on wrong conventions being *grossly* wrong. Where that
    premise holds, the winner is beyond doubt; where it fails, picking the top
    row would be picking noise. This is true only when every axis is separately
    resolved — see :func:`resolved_axes`.
    """
    return all(resolved_axes(scores, margin=margin).values())


def _local_metres(lonlat: np.ndarray, lat0_deg: float) -> np.ndarray:
    """Project longitude/latitude onto a local equirectangular plane, in metres.

    Adequate because the comparison is over a single footprint, and because the
    distances that matter here differ by orders of magnitude, not percent.
    """
    lon = np.radians(lonlat[..., 0])
    lat = np.radians(lonlat[..., 1])
    x = MEAN_EARTH_RADIUS_M * lon * np.cos(np.radians(lat0_deg))
    y = MEAN_EARTH_RADIUS_M * lat
    return np.stack([x, y], axis=-1)


def _dedupe_closing_vertex(polygon: np.ndarray) -> np.ndarray:
    """Drop a repeated final vertex, as GeoJSON rings carry."""
    if len(polygon) > 1 and np.allclose(polygon[0], polygon[-1]):
        return polygon[:-1]
    return polygon


# -- the parities, which only image content can settle ----------------------


def resolve_parities(
    camera: Camera,
    ephemeris: Ephemeris,
    attitude: Attitude,
    timing: dict[str, LineTiming],
    *,
    band: str,
    image: np.ndarray,
    terrain,
    base: Convention,
    land_percentile: float = 55.0,
    land_height_m: float = 2.0,
    n_samples: int = 20_000,
    seed: int = 0,
) -> tuple[Convention, dict]:
    """Resolve the scan direction and the detector column sign — against the terrain.

    These two are **mirrors**: one flips the strip north–south, the other east–west.
    Both map the footprint's four corners onto each other and displace every band
    identically, so every purely geometric probe in this module returns bit-identical
    scores for all four combinations. :func:`resolved_axes` reports them as unresolved,
    and it is right to.

    Only **image content** can see them, and it needs a reference with real geolocation.
    The elevation model is one, and the pipeline has already downloaded it: its coastline
    is ground truth. Water is near-black in the near infrared and land is not, so a
    correct parity makes *bright* coincide with *high* — and a mirrored one
    anti-correlates, which is not a subtle failure.

    This is what a user checking the product against a basemap does by eye, made
    automatic and offline.

    Parameters
    ----------
    camera, ephemeris, attitude, timing:
        The sensor model's components.
    band:
        Band whose imagery is supplied. **Use a near-infrared band** if there is one:
        the water/land contrast is what carries the signal.
    image:
        The band's raster, in sensor coordinates.
    terrain:
        Elevation model, with a geoid (see :mod:`spp.geometry.terrain`).
    base:
        The conventions already resolved geometrically. Only the two parities are varied.
    land_percentile:
        Radiance percentile above which a sample is called land.
    land_height_m:
        Orthometric height above which the terrain is called land.

    Returns
    -------
    tuple
        ``(convention, report)`` — the base convention with the parities filled in, and
        the agreement and Matthews correlation of every candidate.
    """
    rng = np.random.default_rng(seed)
    n_lines, n_columns = image.shape
    lines = rng.uniform(1, n_lines - 2, n_samples)
    columns = rng.uniform(1, n_columns - 2, n_samples)

    radiance = image[lines.astype(int), columns.astype(int)].astype(np.float64)
    finite = np.isfinite(radiance)
    lines, columns, radiance = lines[finite], columns[finite], radiance[finite]
    if radiance.size < 100:
        raise ValueError("Too few valid samples to resolve the parities")

    bright = radiance > np.percentile(radiance, land_percentile)

    report: dict = {}
    best: tuple[float, Convention] | None = None
    for scan in SCAN_DIRECTIONS:
        for column_sign in AXIS_SIGNS:
            candidate = Convention(
                ephemeris_frame=base.ephemeris_frame,
                attitude_frame=base.attitude_frame,
                quaternion_order=base.quaternion_order,
                quaternion_direction=base.quaternion_direction,
                scan_direction=scan,
                column_axis=base.column_axis,
                column_sign=column_sign,
                row_sign=base.row_sign,
            )
            model = SensorModel(camera, ephemeris, attitude, timing, convention=candidate)
            ground = model.locate(band, lines, columns, terrain=terrain)

            orthometric = terrain.height(ground[:, 0], ground[:, 1]) - terrain.geoid(
                ground[:, 0], ground[:, 1]
            )
            high = orthometric > land_height_m

            agreement, mcc = _binary_agreement(bright, high)
            report[f"scan={scan:+d},column_sign={column_sign:+d}"] = {
                "agreement": round(float(agreement), 4),
                "matthews_correlation": round(float(mcc), 4),
            }
            if best is None or mcc > best[0]:
                best = (mcc, candidate)

    assert best is not None
    mcc, winner = best
    report["winner"] = winner.describe()
    report["matthews_correlation"] = round(float(mcc), 4)
    report["decisive"] = bool(mcc > 0.5)

    if mcc <= 0.5:
        logger.warning(
            "Parity resolution is weak (Matthews correlation %.2f). The scene may lack "
            "the land/water contrast this test needs. The parities remain assumptions.",
            mcc,
        )
    else:
        logger.info(
            "Parities resolved against the terrain: scan_direction=%+d, column_sign=%+d "
            "(Matthews correlation %.2f)",
            winner.scan_direction,
            winner.column_sign,
            mcc,
        )
    return winner, report


def _binary_agreement(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    """Agreement and Matthews correlation between two binary classifications.

    The Matthews correlation, not the plain agreement: it is symmetric under class
    imbalance, and a mirrored image does not merely score *lower*, it scores *negative* —
    which is a far stronger statement than "less good".
    """
    tp = float(np.sum(a & b))
    tn = float(np.sum(~a & ~b))
    fp = float(np.sum(a & ~b))
    fn = float(np.sum(~a & b))

    agreement = (tp + tn) / a.size
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = (tp * tn - fp * fn) / denominator if denominator > 0 else 0.0
    return agreement, mcc
