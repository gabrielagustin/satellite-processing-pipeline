"""The sensor model, and the harness that resolves the telemetry's conventions.

Everything here runs on a **synthetic** acquisition whose geometry is known by
construction: a circular orbit, a nadir-pointing platform, and a pushbroom
detector whose bands are staggered in time exactly as a real product staggers
them. That lets the tests assert against ground truth rather than against the
model's own output — which is the only way to catch a model that is
self-consistently wrong.
"""

from __future__ import annotations

import numpy as np
import pytest

from spp.geometry.camera import BandOptics, Camera, CameraIntrinsics
from spp.geometry.conventions import (
    band_coherence,
    candidate_conventions,
    footprint_corners,
    is_decisive,
    resolve,
    resolved_axes,
)
from spp.geometry.ephemeris import Attitude, Ephemeris
from spp.geometry.frames import (
    ecef_to_geodetic,
    lvlh_to_parent_matrix,
    quat_from_rotvec,
    quat_to_rotvec,
)
from spp.geometry.sensor_model import Convention, SensorModel, eci_to_ecef_matrix
from spp.geometry.timing import LineTiming

OMEGA_EARTH = 7.292115e-5
"""Earth rotation rate, radians per second."""

# The convention the synthetic acquisition is built with. The harness must find
# this, and nothing else.
TRUTH = Convention(
    ephemeris_frame="eci",
    attitude_frame="lvlh",
    quaternion_order="scalar_first",
    quaternion_direction="body_to_ref",
    scan_direction=1,
    column_axis="y",
    column_sign=1,
    row_sign=-1,
)

ALTITUDE_M = 400_000.0
RADIUS_M = 6_378_137.0 + ALTITUDE_M
N_LINES = 2000
N_COLUMNS = 512
FOCAL_MM = 580.0
PITCH_MM = 0.0055
PRINCIPAL = (N_COLUMNS / 2.0, 1536.0)
EPOCH = 1_774_685_700.0

# Two bands, far apart on the detector — the pair the coherence probe leans on.
BAND_ROWS = {"FWD": 700, "AFT": 2100}


def _ground_speed_ecef() -> float:
    """Speed of the sub-satellite point over the ground, metres per second.

    The *ground* speed, not the inertial one: the Earth turns beneath the orbit,
    and it is the ground track that the detector sweeps.
    """
    omega = np.sqrt(3.986004418e14 / RADIUS_M**3)
    position = np.array([[RADIUS_M, 0.0, 0.0]])
    velocity = np.array([[0.0, 0.0, RADIUS_M * omega]])

    r = eci_to_ecef_matrix(np.array([EPOCH]))
    r_ecef = np.einsum("mij,mj->mi", r, position)
    v_ecef = np.einsum("mij,mj->mi", r, velocity) - np.cross(
        np.array([0.0, 0.0, OMEGA_EARTH]), r_ecef
    )
    return float(np.linalg.norm(v_ecef) * 6_378_137.0 / RADIUS_M)


def _line_period_s() -> float:
    """The design condition: one line period advances one detector row of ground.

    A pushbroom is built so that the ground moves by exactly one row's worth of
    field between readouts. That is what lets a product stagger its bands by an
    integer number of lines and have their view angles cancel — and it is what
    the sensor model must reproduce.
    """
    ground_per_row = ALTITUDE_M * PITCH_MM / FOCAL_MM
    return ground_per_row / _ground_speed_ecef()


def _yaw_steering(t: np.ndarray, positions: np.ndarray, velocities: np.ndarray) -> np.ndarray:
    """Yaw the platform so its detector rows follow the *ground* track.

    On a rotating Earth the ground track is not the inertial track: the surface
    slides eastward beneath the orbit at up to 465 m/s. A pushbroom whose detector
    rows are not aligned with the ground track has its band stagger smeared
    sideways — the bands would fly along one line and image along another — so
    real platforms yaw by a few degrees to compensate.

    This is not synthetic-data decoration. It is why the delivered attitude is a
    small rotation rather than the identity, and without it the stagger does not
    cancel and the bands do not co-register.
    """
    r = eci_to_ecef_matrix(t)
    r_ecef = np.einsum("mij,mj->mi", r, positions)
    v_ecef = np.einsum("mij,mj->mi", r, velocities) - np.cross(
        np.array([0.0, 0.0, OMEGA_EARTH]), r_ecef
    )
    # The ground-track velocity, resolved back onto the inertial frame's axes.
    v_ground = np.einsum("mji,mj->mi", r, v_ecef)

    lvlh = lvlh_to_parent_matrix(positions, velocities)  # columns: x, y, z in ECI
    along = np.einsum("mi,mi->m", v_ground, lvlh[:, :, 0])
    across = np.einsum("mi,mi->m", v_ground, lvlh[:, :, 1])
    yaw = np.arctan2(across, along)

    return quat_from_rotvec(yaw[:, None] * np.array([0.0, 0.0, 1.0]))


def _synthetic_acquisition() -> tuple[Camera, Ephemeris, Attitude, dict[str, LineTiming]]:
    """A circular, yaw-steered pass with two bands staggered as a real product staggers them."""
    period_s = _line_period_s()
    duration = N_LINES * period_s

    # Telemetry brackets the imaging window, as real telemetry does.
    t = np.arange(EPOCH - 5.0, EPOCH + duration + 5.0, 0.25)
    omega = np.sqrt(3.986004418e14 / RADIUS_M**3)
    angle = omega * (t - EPOCH)
    positions = RADIUS_M * np.stack([np.cos(angle), np.zeros_like(angle), np.sin(angle)], -1)
    velocities = (
        RADIUS_M * omega * np.stack([-np.sin(angle), np.zeros_like(angle), np.cos(angle)], -1)
    )
    ephemeris = Ephemeris(times=t, positions=positions, velocities=velocities)

    quaternions = _yaw_steering(t, positions, velocities)
    rates = np.gradient(quat_to_rotvec(quaternions), t, axis=0)
    attitude = Attitude(times=t, quaternions=quaternions, rates=rates)

    intrinsics = CameraIntrinsics(
        focal_length_mm=FOCAL_MM,
        pixel_size_mm=PITCH_MM,
        principal_point_px=PRINCIPAL,
        sensor_size_px=(N_COLUMNS, 3072),
        detector_to_body=np.diag([-1.0, -1.0, 1.0]),
    )
    camera = Camera(
        intrinsics,
        {name: BandOptics(name=name, start_row=row) for name, row in BAND_ROWS.items()},
    )

    # The stagger: a band's line times are offset by the rows between it and the
    # reference, so its forward/backward view angle is cancelled by the platform
    # having flown that far.
    reference_row = max(BAND_ROWS.values())
    timing = {
        name: LineTiming(
            times=EPOCH
            + (np.arange(N_LINES) + (reference_row - row)) * period_s,
            nominal_period_s=period_s,
        )
        for name, row in BAND_ROWS.items()
    }
    return camera, ephemeris, attitude, timing


@pytest.fixture(scope="module")
def acquisition():
    return _synthetic_acquisition()


@pytest.fixture(scope="module")
def truth_model(acquisition):
    camera, ephemeris, attitude, timing = acquisition
    return SensorModel(camera, ephemeris, attitude, timing, convention=TRUTH)


@pytest.fixture(scope="module")
def truth_footprint(truth_model):
    return footprint_corners(truth_model, "AFT", n_lines=N_LINES, n_columns=N_COLUMNS)


# -- the sensor model -------------------------------------------------------


def test_nadir_column_lands_beneath_the_platform(truth_model):
    """The centre column must land at the sub-satellite point, not near it."""
    line = np.array([N_LINES // 2])
    column = np.array([PRINCIPAL[0]])

    origins, _ = truth_model.rays_ecef("AFT", line, column)
    ground = truth_model.locate("AFT", line, column)

    subsatellite = ecef_to_geodetic(origins[0])
    assert abs(ground[0, 0] - subsatellite[0]) < 0.02
    assert abs(ground[0, 1] - subsatellite[1]) < 0.02


def test_swath_width_matches_the_optics(truth_model):
    """Ground swath must equal 2 * altitude * tan(half-field), from first principles."""
    lines = np.array([N_LINES // 2, N_LINES // 2])
    columns = np.array([0.0, N_COLUMNS - 1.0])
    ground = truth_model.locate("AFT", lines, columns)

    half_field = np.arctan(0.5 * N_COLUMNS * PITCH_MM / FOCAL_MM)
    expected = 2.0 * ALTITUDE_M * np.tan(half_field)

    # Great-circle distance between the two swath edges.
    lat = np.radians(ground[:, 1])
    dx = np.radians(ground[1, 0] - ground[0, 0]) * 6_371_008.8 * np.cos(lat.mean())
    dy = np.radians(ground[1, 1] - ground[0, 1]) * 6_371_008.8
    assert abs(np.hypot(dx, dy) - expected) / expected < 0.02


def test_the_stagger_co_registers_the_bands(truth_model):
    """The whole design, in one assertion.

    Two bands look at different along-track angles and are exposed at different
    times. Under the true convention those two effects cancel, so the same raster
    sample lands on the same ground point in both bands. This is what the L1C
    level is *for*, and if it does not hold, nothing downstream can be trusted.
    """
    coherence = band_coherence(
        truth_model, "FWD", "AFT", n_lines=N_LINES, n_columns=N_COLUMNS
    )
    ground_pixel_m = ALTITUDE_M * PITCH_MM / FOCAL_MM
    assert coherence < 2.0 * ground_pixel_m


def test_flipping_the_row_sign_makes_the_bands_diverge(truth_model, acquisition):
    """With the view angle reversed, the stagger reinforces instead of cancelling.

    This is the failure the harness exists to catch, and it is worth pinning: the
    bands do not merely drift, they land *twice* as far apart as if the view angle
    had been ignored entirely.
    """
    camera, ephemeris, attitude, timing = acquisition
    flipped = SensorModel(
        camera,
        ephemeris,
        attitude,
        timing,
        convention=Convention(**{**TRUTH.__dict__, "row_sign": 1}),
    )
    coherence = band_coherence(
        flipped, "FWD", "AFT", n_lines=N_LINES, n_columns=N_COLUMNS
    )

    row_gap = BAND_ROWS["AFT"] - BAND_ROWS["FWD"]
    naive_separation = ALTITUDE_M * np.tan(np.arctan(row_gap * PITCH_MM / FOCAL_MM))
    assert coherence > 1.5 * naive_separation


def test_each_band_needs_its_own_timing(acquisition):
    """Reusing one band's line times for another is a bug, not a shortcut."""
    camera, ephemeris, attitude, timing = acquisition
    model = SensorModel(camera, ephemeris, attitude, {"AFT": timing["AFT"]})
    with pytest.raises(KeyError, match="No line timing for band"):
        model.locate("FWD", np.array([0]), np.array([0]))


# -- the convention harness -------------------------------------------------


def test_harness_recovers_the_true_convention(acquisition, truth_footprint):
    """The central claim of Phase 0: the conventions are measurable, not guesswork."""
    camera, ephemeris, attitude, timing = acquisition
    scores = resolve(
        camera,
        ephemeris,
        attitude,
        timing,
        band="AFT",
        n_lines=N_LINES,
        n_columns=N_COLUMNS,
        reference_footprint=truth_footprint,
        coherence_bands=("FWD", "AFT"),
    )
    best = scores[0].convention

    # The axes the geometry can see. The parities it cannot (see below) are
    # deliberately excluded: demanding them here would be demanding the harness
    # invent information the data does not contain.
    assert best.ephemeris_frame == TRUTH.ephemeris_frame
    assert best.attitude_frame == TRUTH.attitude_frame
    assert best.quaternion_order == TRUTH.quaternion_order
    assert best.column_axis == TRUTH.column_axis
    assert best.row_sign == TRUTH.row_sign


def test_harness_resolves_the_frames_decisively(acquisition, truth_footprint):
    """A win is not enough; the margin must be large or the answer is noise."""
    camera, ephemeris, attitude, timing = acquisition
    scores = resolve(
        camera, ephemeris, attitude, timing,
        band="AFT", n_lines=N_LINES, n_columns=N_COLUMNS,
        reference_footprint=truth_footprint, coherence_bands=("FWD", "AFT"),
    )
    axes = resolved_axes(scores)
    assert axes["ephemeris_frame"]
    assert axes["attitude_frame"]
    assert axes["column_axis"]
    assert axes["row_sign"]


def test_harness_admits_what_it_cannot_see(acquisition, truth_footprint):
    """The parities are invisible to geometry, and the harness must say so.

    Reversing the scan mirrors the strip north-south; flipping the column sign
    mirrors it east-west. Both map the footprint's four corners onto each other,
    and both displace every band identically so band coherence cancels them too.
    No amount of geometry recovers them — only image content can. A harness that
    claimed otherwise would be fabricating a result, so this test pins the
    *admission* rather than the answer.
    """
    camera, ephemeris, attitude, timing = acquisition
    scores = resolve(
        camera, ephemeris, attitude, timing,
        band="AFT", n_lines=N_LINES, n_columns=N_COLUMNS,
        reference_footprint=truth_footprint, coherence_bands=("FWD", "AFT"),
    )
    axes = resolved_axes(scores)
    assert not axes["scan_direction"]
    assert not axes["column_sign"]
    assert not is_decisive(scores)


def test_wrong_ephemeris_frame_is_grossly_wrong(acquisition, truth_footprint):
    """Reading an inertial ephemeris as Earth-fixed must fail loudly, not subtly.

    This is the error that cost the real run 5,000 km. It must never be a close
    call, because a close call is what lets it slip through.
    """
    camera, ephemeris, attitude, timing = acquisition
    scores = resolve(
        camera, ephemeris, attitude, timing,
        band="AFT", n_lines=N_LINES, n_columns=N_COLUMNS,
        reference_footprint=truth_footprint, coherence_bands=("FWD", "AFT"),
        candidates=candidate_conventions(ephemeris_frames=("ecef",)),
    )
    assert all(not s.plausible() for s in scores)


def test_harness_scores_every_candidate(acquisition, truth_footprint):
    """No candidate may be silently dropped — a skipped one is an unexamined one."""
    camera, ephemeris, attitude, timing = acquisition
    scores = resolve(
        camera, ephemeris, attitude, timing,
        band="AFT", n_lines=N_LINES, n_columns=N_COLUMNS,
        reference_footprint=truth_footprint, coherence_bands=("FWD", "AFT"),
    )
    assert len(scores) == len(candidate_conventions())
