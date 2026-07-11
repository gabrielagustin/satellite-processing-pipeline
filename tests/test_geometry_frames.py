"""Numerical correctness of the geometric primitives.

These are the load-bearing computations of the L1C level, and every one of them
is checked against something *independent* — a closed form, an analytic inverse,
or a property that must hold regardless of implementation. A rotation library
that is subtly wrong produces a plausible product on the wrong piece of ground,
so "it looks right" is not evidence here.
"""

from __future__ import annotations

import numpy as np
import pytest

from spp.geometry.ephemeris import Attitude, Ephemeris
from spp.geometry.frames import (
    WGS84_A,
    WGS84_B,
    ecef_to_geodetic,
    geodetic_to_ecef,
    lvlh_to_parent_matrix,
    quat_conjugate,
    quat_from_rotvec,
    quat_multiply,
    quat_normalize,
    quat_to_matrix,
    quat_to_rotvec,
    ray_ellipsoid_intersection,
)


# -- quaternion algebra -----------------------------------------------------


def test_quaternion_to_matrix_is_a_rotation():
    """Any unit quaternion must map to an orthonormal, right-handed matrix."""
    rng = np.random.default_rng(0)
    q = quat_normalize(rng.normal(size=(20, 4)))
    r = quat_to_matrix(q)

    identity = np.einsum("nij,nkj->nik", r, r)
    np.testing.assert_allclose(identity, np.broadcast_to(np.eye(3), (20, 3, 3)), atol=1e-12)
    np.testing.assert_allclose(np.linalg.det(r), 1.0, atol=1e-12)


def test_quaternion_rotation_matches_known_case():
    """A 90 degree rotation about z takes x to y."""
    q = quat_from_rotvec(np.array([0.0, 0.0, np.pi / 2]))
    rotated = quat_to_matrix(q) @ np.array([1.0, 0.0, 0.0])
    np.testing.assert_allclose(rotated, [0.0, 1.0, 0.0], atol=1e-12)


def test_quaternion_composition_matches_matrix_composition():
    """`quat_multiply(a, b)` must mean "apply b, then a" — as matrices do.

    Getting this backwards silently transposes the whole body-to-ground chain.
    """
    rng = np.random.default_rng(1)
    a = quat_normalize(rng.normal(size=4))
    b = quat_normalize(rng.normal(size=4))
    np.testing.assert_allclose(
        quat_to_matrix(quat_multiply(a, b)),
        quat_to_matrix(a) @ quat_to_matrix(b),
        atol=1e-12,
    )


def test_conjugate_inverts():
    """`q` composed with its conjugate is the identity rotation."""
    rng = np.random.default_rng(2)
    q = quat_normalize(rng.normal(size=(10, 4)))
    np.testing.assert_allclose(
        quat_to_matrix(quat_multiply(q, quat_conjugate(q))),
        np.broadcast_to(np.eye(3), (10, 3, 3)),
        atol=1e-12,
    )


def test_rotvec_round_trip():
    """exp and log are mutual inverses, over the representable range."""
    rng = np.random.default_rng(3)
    axis = rng.normal(size=(20, 3))
    axis /= np.linalg.norm(axis, axis=-1, keepdims=True)
    angle = rng.uniform(0.0, np.pi - 1e-6, size=(20, 1))
    rotvec = axis * angle
    np.testing.assert_allclose(quat_to_rotvec(quat_from_rotvec(rotvec)), rotvec, atol=1e-10)


def test_rotvec_small_angle_is_stable():
    """The small-angle branch must not divide by zero, and must stay accurate.

    Attitude interpolation lives here: between telemetry samples the relative
    rotation is milliradians, so this is the branch that runs in production, not
    an edge case.
    """
    tiny = np.array([[1e-12, 0.0, 0.0], [0.0, 1e-9, 0.0], [0.0, 0.0, 0.0]])
    q = quat_from_rotvec(tiny)
    assert np.all(np.isfinite(q))
    np.testing.assert_allclose(np.linalg.norm(q, axis=-1), 1.0, atol=1e-12)
    np.testing.assert_allclose(quat_to_rotvec(q), tiny, atol=1e-15)


# -- orbital frame ----------------------------------------------------------


def test_lvlh_triad_is_orthonormal_and_nadir_pointing():
    """z must point at the Earth's centre; the triad must be right-handed."""
    position = np.array([[7000e3, 0.0, 0.0]])
    velocity = np.array([[0.0, 7.5e3, 0.0]])
    m = lvlh_to_parent_matrix(position, velocity)[0]

    np.testing.assert_allclose(m.T @ m, np.eye(3), atol=1e-12)
    np.testing.assert_allclose(np.linalg.det(m), 1.0, atol=1e-12)

    nadir = -position[0] / np.linalg.norm(position[0])
    np.testing.assert_allclose(m[:, 2], nadir, atol=1e-12)  # z column
    # x completes the triad along the ground track, i.e. with the velocity.
    assert float(m[:, 0] @ velocity[0]) > 0.0


def test_lvlh_rejects_a_degenerate_state():
    """A velocity parallel to the radius leaves the orbit normal undefined."""
    with pytest.raises(ValueError, match="Degenerate"):
        lvlh_to_parent_matrix(np.array([[7000e3, 0.0, 0.0]]), np.array([[1.0, 0.0, 0.0]]))


# -- Earth geometry ---------------------------------------------------------


def test_geodetic_round_trip():
    """ecef_to_geodetic must invert geodetic_to_ecef to millimetres."""
    rng = np.random.default_rng(4)
    g = np.stack(
        [
            rng.uniform(-180.0, 180.0, 50),
            rng.uniform(-85.0, 85.0, 50),
            rng.uniform(-500.0, 9000.0, 50),
        ],
        axis=-1,
    )
    np.testing.assert_allclose(ecef_to_geodetic(geodetic_to_ecef(g)), g, atol=1e-6)


def test_geodetic_poles_and_equator():
    """The ellipsoid's own axes are the cases a wrong formula gets wrong."""
    equator = ecef_to_geodetic(np.array([WGS84_A, 0.0, 0.0]))
    np.testing.assert_allclose(equator, [0.0, 0.0, 0.0], atol=1e-6)

    pole = ecef_to_geodetic(np.array([0.0, 0.0, WGS84_B]))
    np.testing.assert_allclose(pole[1:], [90.0, 0.0], atol=1e-6)


def test_ray_hits_the_ellipsoid_beneath_the_platform():
    """A nadir ray lands at the platform's own longitude and latitude."""
    platform = geodetic_to_ecef(np.array([30.0, 45.0, 500e3]))
    direction = -platform / np.linalg.norm(platform)  # straight down-ish

    hit = ecef_to_geodetic(ray_ellipsoid_intersection(platform, direction))
    # A geocentric nadir ray is not exactly geodetic nadir, so latitude differs
    # slightly; the height is what pins the intersection to the surface.
    np.testing.assert_allclose(hit[2], 0.0, atol=1e-6)
    np.testing.assert_allclose(hit[0], 30.0, atol=1e-9)


def test_ray_intersection_returns_the_near_side():
    """The ray must stop at the visible surface, not exit through the far side."""
    platform = np.array([10_000e3, 0.0, 0.0])
    hit = ray_ellipsoid_intersection(platform, np.array([-1.0, 0.0, 0.0]))
    np.testing.assert_allclose(hit, [WGS84_A, 0.0, 0.0], atol=1e-6)


def test_ray_at_height_lands_at_that_height():
    """Intersecting an inflated ellipsoid must return that height."""
    platform = geodetic_to_ecef(np.array([0.0, 20.0, 500e3]))
    direction = -platform / np.linalg.norm(platform)
    hit = ecef_to_geodetic(ray_ellipsoid_intersection(platform, direction, height=1000.0))
    assert abs(hit[2] - 1000.0) < 1.0


def test_ray_that_misses_the_earth_is_nan():
    """A ray pointed at the sky must not silently produce a ground point."""
    platform = np.array([10_000e3, 0.0, 0.0])
    hit = ray_ellipsoid_intersection(platform, np.array([1.0, 0.0, 0.0]))
    assert np.all(np.isnan(hit))


# -- state interpolation ----------------------------------------------------


def _circular_orbit(times: np.ndarray, radius: float = 6_871_000.0):
    """An exact circular orbit — a trajectory whose truth we know everywhere."""
    omega = np.sqrt(3.986004418e14 / radius**3)
    angle = omega * times
    positions = radius * np.stack([np.cos(angle), np.sin(angle), np.zeros_like(angle)], -1)
    velocities = (
        radius * omega * np.stack([-np.sin(angle), np.cos(angle), np.zeros_like(angle)], -1)
    )
    return positions, velocities


def test_ephemeris_reproduces_its_samples_exactly():
    """Interpolation must be exact at the nodes, or it is not interpolation."""
    t = np.arange(0.0, 10.0, 0.25)
    p, v = _circular_orbit(t)
    eph = Ephemeris(times=t, positions=p, velocities=v)

    pi, vi = eph.interpolate(t)
    np.testing.assert_allclose(pi, p, atol=1e-6)
    np.testing.assert_allclose(vi, v, atol=1e-6)


def test_hermite_beats_linear_on_a_known_curve():
    """Using the delivered velocities must actually buy accuracy.

    The claim in the spec is that rate-aware interpolation is worth it. On an
    exact circular orbit, sampled as coarsely as the real telemetry, Hermite must
    land orders of magnitude closer to the truth than a straight line between
    samples — otherwise the extra machinery earns nothing.
    """
    t_nodes = np.arange(0.0, 10.0, 0.25)  # ~4 Hz, as delivered
    p, v = _circular_orbit(t_nodes)
    eph = Ephemeris(times=t_nodes, positions=p, velocities=v)

    t_mid = t_nodes[:-1] + 0.125  # worst case: halfway between samples
    truth, _ = _circular_orbit(t_mid)

    hermite, _ = eph.interpolate(t_mid)
    linear = np.stack([np.interp(t_mid, t_nodes, p[:, i]) for i in range(3)], axis=-1)

    hermite_err = np.abs(hermite - truth).max()
    linear_err = np.abs(linear - truth).max()
    assert hermite_err < linear_err / 100.0
    assert hermite_err < 1e-3  # sub-millimetre


def test_attitude_reproduces_its_samples_exactly():
    """Both interpolation paths must pass through the sampled quaternions."""
    t = np.arange(0.0, 5.0, 0.25)
    rate = np.array([0.0, 0.01, 0.0])
    q = quat_from_rotvec(rate * t[:, None])
    rates = np.broadcast_to(rate, (t.size, 3))
    att = Attitude(times=t, quaternions=q, rates=rates)

    for rate_aware in (True, False):
        got = att.interpolate(t, rate_aware=rate_aware)
        # q and -q are the same rotation, so compare rotations, not components.
        np.testing.assert_allclose(quat_to_matrix(got), quat_to_matrix(q), atol=1e-9)


def test_attitude_interpolates_a_constant_rotation_exactly():
    """Under a constant body rate, the attitude at any time is known in closed form."""
    t = np.arange(0.0, 5.0, 0.25)
    rate = np.array([0.0, 0.0, 0.02])
    q = quat_from_rotvec(rate * t[:, None])
    att = Attitude(times=t, quaternions=q, rates=np.broadcast_to(rate, (t.size, 3)))

    t_mid = t[:-1] + 0.125
    expected = quat_to_matrix(quat_from_rotvec(rate * t_mid[:, None]))
    np.testing.assert_allclose(quat_to_matrix(att.interpolate(t_mid)), expected, atol=1e-9)


def test_state_histories_reject_bad_input():
    """Unordered or ragged telemetry is a bug upstream; fail loudly, not quietly."""
    t = np.array([0.0, 2.0, 1.0])
    p = np.zeros((3, 3))
    with pytest.raises(ValueError, match="strictly increasing"):
        Ephemeris(times=t, positions=p, velocities=p)

    with pytest.raises(ValueError, match="parallel"):
        Ephemeris(
            times=np.array([0.0, 1.0, 2.0]),
            positions=np.zeros((2, 3)),
            velocities=np.zeros((3, 3)),
        )


# -- angular rates that are not the attitude's derivative -------------------


def test_consistent_rates_are_recognised():
    """Rates that really are the quaternion's slope must be accepted and used."""
    t = np.arange(0.0, 5.0, 0.25)
    rate = np.array([0.0, 0.0, 0.02])
    att = Attitude(
        times=t,
        quaternions=quat_from_rotvec(rate * t[:, None]),
        rates=np.broadcast_to(rate, (t.size, 3)),
    )
    assert att.rates_are_consistent


def test_inconsistent_rates_are_rejected_and_slerp_takes_over():
    """The trap the reference acquisition set, pinned.

    The delivered angular rates were ~55x larger than the rotation the quaternion
    sequence actually undergoes: they describe motion in a different frame. Used
    as Hermite endpoint slopes they overshoot wildly, and the model gets *worse*.

    The rates must therefore be checked against the quaternions and ignored when
    they disagree — and asking for `rate_aware=True` must not be able to override
    that, because the caller cannot know what the data contains.
    """
    t = np.arange(0.0, 5.0, 0.25)
    true_rate = np.array([0.0, 0.0, 0.002])
    quaternions = quat_from_rotvec(true_rate * t[:, None])
    bogus = np.broadcast_to(true_rate * 55.0, (t.size, 3))

    att = Attitude(times=t, quaternions=quaternions, rates=bogus)
    assert not att.rates_are_consistent

    # Asking for rate-awareness must yield the SLERP answer, not the overshoot.
    slerp = Attitude(times=t, quaternions=quaternions).interpolate(t[:-1] + 0.125)
    guarded = att.interpolate(t[:-1] + 0.125, rate_aware=True)
    np.testing.assert_allclose(quat_to_matrix(guarded), quat_to_matrix(slerp), atol=1e-12)


def test_the_overshoot_the_guard_prevents_is_real():
    """The guard is not defensive boilerplate — without it the error is large.

    Pinning the magnitude keeps the rationale honest: if bogus rates were harmless,
    the check would be dead weight and should be deleted.

    The mechanism is worth stating, because it is not the obvious one. A wrong rate
    that merely has the wrong *magnitude* is largely harmless: the two Hermite
    endpoint terms carry weights `h10` and `h11` that sum to `s(2s-1)(s-1)`, which
    vanishes at the midpoint, so a constant too-steep slope mostly cancels itself.
    The damage comes from a rate pointing along a **different axis** — which is
    exactly what the reference acquisition delivers, its rates describing the body's
    motion in some other frame rather than the drift of the attitude offset. That
    injects rotation the quaternions never contained, and it does not cancel.
    """
    t = np.arange(0.0, 5.0, 0.25)
    true_rate = np.array([0.0, 0.0, 0.002])
    quaternions = quat_from_rotvec(true_rate * t[:, None])

    # As delivered: far larger, and about a different axis.
    bogus_rates = np.broadcast_to(np.array([0.0, 0.11, 0.0]), (t.size, 3))

    # Sample away from the midpoint, where the endpoint weights do not cancel.
    dt = 0.25
    t_eval = t[:-1] + 0.21 * dt
    truth = quat_from_rotvec(true_rate * t_eval[:, None])

    def peak_error(got: np.ndarray) -> float:
        """Largest rotation angle between `got` and the truth, radians."""
        delta = quat_multiply(quat_conjugate(truth), got)
        return float(np.linalg.norm(quat_to_rotvec(delta), axis=-1).max())

    # What the guard produces: SLERP, exact on a constant-rate rotation.
    guarded = Attitude(times=t, quaternions=quaternions, rates=bogus_rates)
    assert not guarded.rates_are_consistent
    assert peak_error(guarded.interpolate(t_eval, rate_aware=True)) < 1e-9

    # Rates that *are* the derivative are accepted, and Hermite is exact too.
    honest = Attitude(
        times=t, quaternions=quaternions, rates=np.broadcast_to(true_rate, (t.size, 3))
    )
    assert honest.rates_are_consistent
    assert peak_error(honest.interpolate(t_eval, rate_aware=True)) < 1e-9

    # And what believing the delivered rates would have cost: an error larger than
    # the entire rotation being interpolated.
    segment_rotation = float(np.linalg.norm(true_rate) * dt)
    forced = _hermite_with_rates(quaternions, t, bogus_rates, t_eval)
    assert peak_error(forced) > 2.0 * segment_rotation


def _hermite_with_rates(
    quaternions: np.ndarray, t: np.ndarray, rates: np.ndarray, t_eval: np.ndarray
) -> np.ndarray:
    """Rate-aware interpolation with the consistency guard deliberately bypassed.

    Only a test has any business doing this — it exists to show what the guard is
    for, by reproducing the failure it prevents.
    """
    idx = np.clip(np.searchsorted(t, t_eval, side="right") - 1, 0, t.size - 2)
    dt = (t[idx + 1] - t[idx])[:, None]
    s = ((t_eval - t[idx]) / dt[:, 0])[:, None]

    q0, q1 = quaternions[idx], quaternions[idx + 1]
    delta = quat_to_rotvec(quat_multiply(quat_conjugate(q0), q1))

    h10 = s**3 - 2 * s**2 + s
    h01 = -2 * s**3 + 3 * s**2
    h11 = s**3 - s**2
    rotvec = h10 * dt * rates[idx] + h01 * delta + h11 * dt * rates[idx + 1]
    return quat_normalize(quat_multiply(q0, quat_from_rotvec(rotvec)))
