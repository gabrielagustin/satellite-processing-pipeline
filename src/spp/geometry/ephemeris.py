"""Platform state: ephemeris and attitude, sampled and interpolated.

The platform telemetry arrives as a sparse, time-tagged history — typically a
few samples per second — while the sensor model needs the state at the exact
instant every raster line was exposed. This module holds the sampled states and
interpolates them.

The telemetry delivers derivatives alongside the states — velocity with position,
angular rates with attitude — which makes cubic Hermite interpolation available
for free, and at a few hertz that accuracy is worth having.

**But a derivative is only useful if it is the derivative of the thing you are
interpolating**, and that must be *checked*, not assumed. The position and its
velocity agree, so position is interpolated with Hermite. The attitude and its
angular rates do **not** agree in the reference acquisition — the rates are ~55x
larger than the quaternion sequence's own rotation — so attitude falls back to
SLERP, automatically and with a warning (see :attr:`Attitude.rates_are_consistent`).

Feeding Hermite a slope that is 55x too steep does not degrade the result a
little; it overshoots between every pair of samples. On the reference acquisition
it made the band-to-band error **seven times worse** (18 m to 123 m) — a
"refinement" that quietly wrecked the model. Hence the check.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np

from spp.geometry.frames import (
    quat_conjugate,
    quat_from_rotvec,
    quat_multiply,
    quat_normalize,
    quat_to_rotvec,
)

logger = logging.getLogger(__name__)


def _hermite_basis(s: np.ndarray) -> tuple[np.ndarray, ...]:
    """Cubic Hermite basis functions evaluated at ``s`` in ``[0, 1]``."""
    s2 = s * s
    s3 = s2 * s
    return (
        2 * s3 - 3 * s2 + 1,  # h00 -> value at the left node
        s3 - 2 * s2 + s,      # h10 -> derivative at the left node
        -2 * s3 + 3 * s2,     # h01 -> value at the right node
        s3 - s2,              # h11 -> derivative at the right node
    )


def _segments(sample_times: np.ndarray, t: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Locate ``t`` within ``sample_times``.

    Returns the left node index, the segment duration and the normalised
    position ``s`` within the segment. Times outside the sampled span are
    clamped to the end segments, which extrapolates rather than failing; callers
    are expected to have checked coverage (:meth:`StateHistory.covers`).
    """
    n = sample_times.size
    idx = np.searchsorted(sample_times, t, side="right") - 1
    idx = np.clip(idx, 0, n - 2)
    t0 = sample_times[idx]
    t1 = sample_times[idx + 1]
    dt = t1 - t0
    s = (t - t0) / dt
    return idx, dt, s


@dataclass(frozen=True)
class Ephemeris:
    """Sampled platform position and velocity in an Earth-fixed frame.

    Attributes
    ----------
    times:
        Sample times, shape ``(n,)``, seconds (Unix epoch). Strictly increasing.
    positions:
        ECEF positions, shape ``(n, 3)``, metres.
    velocities:
        ECEF velocities, shape ``(n, 3)``, metres per second.
    """

    times: np.ndarray = field(repr=False)
    positions: np.ndarray = field(repr=False)
    velocities: np.ndarray = field(repr=False)

    def __post_init__(self) -> None:
        _validate_history(self.times, (self.positions, self.velocities), "Ephemeris")

    def covers(self, t: np.ndarray | float) -> bool:
        """Whether ``t`` lies within the sampled span (no extrapolation needed)."""
        t = np.asarray(t, dtype=np.float64)
        return bool(t.min() >= self.times[0] and t.max() <= self.times[-1])

    def interpolate(self, t: np.ndarray | float) -> tuple[np.ndarray, np.ndarray]:
        """Interpolate position and velocity at times ``t``.

        Cubic Hermite on position, using the sampled velocities as the
        derivatives; the velocity returned is the analytic derivative of that
        same cubic, so position and velocity stay mutually consistent.

        Parameters
        ----------
        t:
            Times, seconds (Unix epoch). Scalar or shape ``(m,)``.

        Returns
        -------
        tuple of numpy.ndarray
            ``(positions, velocities)``, each of shape ``(m, 3)``.
        """
        t = np.atleast_1d(np.asarray(t, dtype=np.float64))
        idx, dt, s = _segments(self.times, t)
        h00, h10, h01, h11 = _hermite_basis(s)

        p0, p1 = self.positions[idx], self.positions[idx + 1]
        v0, v1 = self.velocities[idx], self.velocities[idx + 1]
        dt_c = dt[:, None]

        position = (
            h00[:, None] * p0
            + h10[:, None] * dt_c * v0
            + h01[:, None] * p1
            + h11[:, None] * dt_c * v1
        )

        # Derivative of the Hermite cubic, in units of time (not of s).
        s2 = s * s
        d00 = 6 * s2 - 6 * s
        d10 = 3 * s2 - 4 * s + 1
        d01 = -6 * s2 + 6 * s
        d11 = 3 * s2 - 2 * s
        velocity = (
            d00[:, None] * p0 / dt_c
            + d10[:, None] * v0
            + d01[:, None] * p1 / dt_c
            + d11[:, None] * v1
        )
        return position, velocity


@dataclass(frozen=True)
class Attitude:
    """Sampled platform attitude and angular rates.

    The quaternions are stored scalar-first and normalised, in the convention of
    :mod:`spp.geometry.frames`. What frame they rotate *between*, and in which
    direction, is a property of the payload and is resolved separately — see
    :mod:`spp.geometry.conventions`. This class only interpolates them.

    Attributes
    ----------
    times:
        Sample times, shape ``(n,)``, seconds (Unix epoch). Strictly increasing.
    quaternions:
        Unit quaternions, shape ``(n, 4)``, scalar-first.
    rates:
        Angular rates, shape ``(n, 3)``, radians per second, expressed in the
        rotating (body) frame. ``None`` if not delivered, in which case
        interpolation falls back to SLERP.
    """

    times: np.ndarray = field(repr=False)
    quaternions: np.ndarray = field(repr=False)
    rates: np.ndarray | None = field(default=None, repr=False)
    rates_are_consistent: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        arrays = (self.quaternions,) if self.rates is None else (self.quaternions, self.rates)
        _validate_history(self.times, arrays, "Attitude")
        object.__setattr__(self, "quaternions", quat_normalize(self.quaternions))

        # Decided once, here, because it is a property of the *data*, not of any one
        # call. Evaluating it per-interpolation would recompute the whole quaternion
        # sequence on every geolocation node -- and would emit the same warning dozens
        # of times, which is how a warning stops being read.
        object.__setattr__(self, "rates_are_consistent", self._check_rates())
        if self.rates is not None and not self.rates_are_consistent:
            logger.warning(
                "Attitude: the delivered angular rates are not the derivative of the "
                "delivered quaternions, so they cannot be used as interpolation slopes; "
                "falling back to SLERP. Using them would overshoot between samples and "
                "make the model worse, not better."
            )

    def _check_rates(self) -> bool:
        """Whether the delivered angular rates are the derivative of the quaternions.

        They need not be, and in the reference acquisition they are **not**: the
        delivered rates are ~55x larger than the rotation the quaternion sequence
        actually undergoes between samples. They evidently describe the body's
        motion in some other frame (an inertial one, most likely, which for an
        Earth-pointing platform is dominated by the orbital rate) rather than the
        drift of the small attitude offset the quaternions encode.

        This is not a defect in the data, it is an unstated convention — but it is
        a trap. Rate-aware interpolation uses the rates as the endpoint slopes of
        the quaternion curve. Feed it slopes that are 55x too steep and it does not
        degrade gracefully; it overshoots violently between every pair of samples.
        On the reference acquisition that inflated the band-to-band error from 18 m
        to 123 m — a model made worse by a "refinement".

        So the rates are **checked, not trusted**. When they disagree with the
        quaternion sequence, interpolation falls back to SLERP, which needs no
        derivative and cannot be misled by a wrong one.

        Evaluated once, at construction; the answer is stored in
        :attr:`rates_are_consistent`.
        """
        if self.rates is None:
            return False
        implied = quat_to_rotvec(
            quat_multiply(quat_conjugate(self.quaternions[:-1]), self.quaternions[1:])
        ) / np.diff(self.times)[:, None]
        delivered = 0.5 * (self.rates[:-1] + self.rates[1:])

        implied_speed = np.linalg.norm(implied, axis=-1).mean()
        delivered_speed = np.linalg.norm(delivered, axis=-1).mean()
        if implied_speed < 1e-12:
            return delivered_speed < 1e-9  # a static attitude needs a zero rate
        return bool(0.5 <= delivered_speed / implied_speed <= 2.0)

    def covers(self, t: np.ndarray | float) -> bool:
        """Whether ``t`` lies within the sampled span (no extrapolation needed)."""
        t = np.asarray(t, dtype=np.float64)
        return bool(t.min() >= self.times[0] and t.max() <= self.times[-1])

    def interpolate(self, t: np.ndarray | float, *, rate_aware: bool = True) -> np.ndarray:
        """Interpolate attitude at times ``t``.

        Parameters
        ----------
        t:
            Times, seconds (Unix epoch). Scalar or shape ``(m,)``.
        rate_aware:
            Use the sampled angular rates as endpoint derivatives (cubic Hermite
            in rotation-vector space). Falls back to SLERP (Spherical Linear
            intERPolation) when the rates are absent, when set to ``False``, or —
            crucially — when the delivered rates turn out **not** to be the
            derivative of the delivered quaternions (see
            :attr:`rates_are_consistent`). Asking for rate-awareness is a request,
            not an instruction: the data gets a veto.

        Returns
        -------
        numpy.ndarray
            Unit quaternions, shape ``(m, 4)``, scalar-first.

        Notes
        -----
        Both paths reproduce the sampled quaternions **exactly** at the sample
        times. The rate-aware path additionally matches the sampled angular rate
        at each node — which is an improvement only if that rate really is the
        curve's slope, and a serious regression if it is not.
        """
        t = np.atleast_1d(np.asarray(t, dtype=np.float64))
        idx, dt, s = _segments(self.times, t)

        q0 = self.quaternions[idx]
        q1 = self.quaternions[idx + 1]

        # Relative rotation from q0 to q1, as a rotation vector in q0's frame.
        # The shortest-arc convention in quat_to_rotvec resolves the double
        # cover (q and -q are the same rotation), so no sign fix is needed here.
        delta = quat_to_rotvec(quat_multiply(quat_conjugate(q0), q1))

        use_rates = rate_aware and self.rates is not None and self.rates_are_consistent
        if not use_rates:
            rotvec = s[:, None] * delta  # SLERP, expressed in the same algebra
        else:
            w0 = self.rates[idx]
            w1 = self.rates[idx + 1]
            h00, h10, h01, h11 = _hermite_basis(s)
            dt_c = dt[:, None]
            # Hermite on the rotation vector: zero at the left node, `delta` at
            # the right node, with the angular rates as the endpoint slopes.
            rotvec = (
                h00[:, None] * 0.0
                + h10[:, None] * dt_c * w0
                + h01[:, None] * delta
                + h11[:, None] * dt_c * w1
            )

        return quat_normalize(quat_multiply(q0, quat_from_rotvec(rotvec)))


def _validate_history(times: np.ndarray, arrays: tuple[np.ndarray, ...], what: str) -> None:
    """Guard the shared invariants of a sampled state history."""
    if times.ndim != 1 or times.size < 2:
        raise ValueError(f"{what} needs at least two time samples, got {times.size}")
    if not np.all(np.diff(times) > 0):
        raise ValueError(f"{what} sample times must be strictly increasing")
    for a in arrays:
        if a.shape[0] != times.size:
            raise ValueError(
                f"{what} sample arrays must be parallel to the times "
                f"({times.size}), got {a.shape[0]}"
            )
