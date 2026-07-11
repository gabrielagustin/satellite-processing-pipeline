"""Reference frames, rotations and Earth geometry.

The geometric level composes a chain of rotations to turn a detector sample into
a viewing ray in an Earth-fixed frame::

    detector --> camera --> body --> LVLH --> ECEF

where LVLH is the Local Vertical, Local Horizontal orbital frame and ECEF is the
Earth-Centred, Earth-Fixed frame. This module provides the primitives for that
chain (quaternion algebra, frame construction, ellipsoid intersection and
geodetic conversion) and nothing else: it is pure NumPy, has no I/O, and knows
nothing about acquisitions.

Quaternion convention
---------------------
Quaternions are held **scalar-first**, ``[w, x, y, z]``, normalised, and are
interpreted as *active* rotations: :func:`quat_to_matrix` returns the matrix
``R`` such that ``R @ v`` rotates the vector ``v``. Composition follows
``quat_multiply(a, b)`` = "apply ``b`` first, then ``a``", matching matrix
composition ``R(a) @ R(b)``.

The *input* convention of a given payload (component order, rotation direction,
reference frame) is a separate question, and a deliberately open one — see
:mod:`spp.geometry.conventions`.
"""

from __future__ import annotations

import numpy as np

# WGS84 (World Geodetic System 1984) ellipsoid.
WGS84_A = 6378137.0
"""Semi-major axis (equatorial radius), metres."""

WGS84_F = 1.0 / 298.257223563
"""Flattening."""

WGS84_B = WGS84_A * (1.0 - WGS84_F)
"""Semi-minor axis (polar radius), metres."""

WGS84_E2 = WGS84_F * (2.0 - WGS84_F)
"""First eccentricity squared."""


# -- quaternion algebra -----------------------------------------------------


def quat_normalize(q: np.ndarray) -> np.ndarray:
    """Normalise quaternions to unit norm.

    Parameters
    ----------
    q:
        Array of shape ``(..., 4)``, scalar-first.

    Returns
    -------
    numpy.ndarray
        Unit quaternions of the same shape.

    Raises
    ------
    ValueError
        If any quaternion has (near-)zero norm and cannot be normalised.
    """
    q = np.asarray(q, dtype=np.float64)
    norm = np.linalg.norm(q, axis=-1, keepdims=True)
    if np.any(norm < 1e-12):
        raise ValueError("Cannot normalise a zero-norm quaternion")
    return q / norm


def quat_conjugate(q: np.ndarray) -> np.ndarray:
    """Return the conjugate (inverse, for unit quaternions) of ``q``."""
    q = np.asarray(q, dtype=np.float64)
    out = q.copy()
    out[..., 1:] *= -1.0
    return out


def quat_multiply(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Hamilton product ``a ⊗ b`` (apply ``b`` first, then ``a``).

    Parameters
    ----------
    a, b:
        Arrays of shape ``(..., 4)``, scalar-first, broadcastable against each
        other.
    """
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    aw, ax, ay, az = a[..., 0], a[..., 1], a[..., 2], a[..., 3]
    bw, bx, by, bz = b[..., 0], b[..., 1], b[..., 2], b[..., 3]
    return np.stack(
        [
            aw * bw - ax * bx - ay * by - az * bz,
            aw * bx + ax * bw + ay * bz - az * by,
            aw * by - ax * bz + ay * bw + az * bx,
            aw * bz + ax * by - ay * bx + az * bw,
        ],
        axis=-1,
    )


def quat_to_matrix(q: np.ndarray) -> np.ndarray:
    """Convert unit quaternions to rotation matrices.

    Parameters
    ----------
    q:
        Array of shape ``(..., 4)``, scalar-first. Normalised internally.

    Returns
    -------
    numpy.ndarray
        Array of shape ``(..., 3, 3)``. ``R @ v`` applies the rotation to ``v``.
    """
    q = quat_normalize(q)
    w, x, y, z = q[..., 0], q[..., 1], q[..., 2], q[..., 3]
    return np.stack(
        [
            np.stack(
                [
                    1 - 2 * (y * y + z * z),
                    2 * (x * y - w * z),
                    2 * (x * z + w * y),
                ],
                axis=-1,
            ),
            np.stack(
                [
                    2 * (x * y + w * z),
                    1 - 2 * (x * x + z * z),
                    2 * (y * z - w * x),
                ],
                axis=-1,
            ),
            np.stack(
                [
                    2 * (x * z - w * y),
                    2 * (y * z + w * x),
                    1 - 2 * (x * x + y * y),
                ],
                axis=-1,
            ),
        ],
        axis=-2,
    )


def quat_from_rotvec(rotvec: np.ndarray) -> np.ndarray:
    """Exponential map: rotation vector (axis · angle, radians) to quaternion.

    Parameters
    ----------
    rotvec:
        Array of shape ``(..., 3)``.

    Returns
    -------
    numpy.ndarray
        Unit quaternions of shape ``(..., 4)``, scalar-first.
    """
    rotvec = np.asarray(rotvec, dtype=np.float64)
    angle = np.linalg.norm(rotvec, axis=-1, keepdims=True)
    half = 0.5 * angle
    # sin(half)/angle, guarding the small-angle limit (-> 1/2).
    small = angle < 1e-9
    scale = np.where(small, 0.5 - angle**2 / 48.0, np.sin(half) / np.where(small, 1.0, angle))
    return np.concatenate([np.cos(half), scale * rotvec], axis=-1)


def quat_to_rotvec(q: np.ndarray) -> np.ndarray:
    """Logarithmic map: quaternion to rotation vector (axis · angle, radians).

    The shortest-arc representative is always returned (the sign of ``q`` is
    normalised so that the scalar part is non-negative), so the result has an
    angle in ``[0, π]``.
    """
    q = quat_normalize(q)
    q = np.where(q[..., :1] < 0.0, -q, q)  # shortest arc
    vec = q[..., 1:]
    sin_half = np.linalg.norm(vec, axis=-1, keepdims=True)
    cos_half = q[..., :1]
    angle = 2.0 * np.arctan2(sin_half, cos_half)
    small = sin_half < 1e-9
    scale = np.where(small, 2.0, angle / np.where(small, 1.0, sin_half))
    return scale * vec


# -- orbital frame ----------------------------------------------------------


def lvlh_to_parent_matrix(position: np.ndarray, velocity: np.ndarray) -> np.ndarray:
    """Build the rotation from the LVLH orbital frame to the state's own frame.

    The LVLH (Local Vertical, Local Horizontal) triad used here is::

        z  =  nadir                    = -normalize(r)
        y  = -orbit normal             = -normalize(r x v)
        x  =  y x z                    (completes the right-handed triad,
                                        pointing along the ground track)

    which is the common "z-down, x-forward" convention for an Earth-pointing
    platform. The returned matrix has these axes as its **columns**, so it maps
    a vector expressed in LVLH into the parent frame.

    The **parent frame is whichever frame the state is expressed in** — inertial
    or Earth-fixed. That is not a detail to gloss over: an orbital frame built
    from an inertial velocity and one built from an Earth-fixed velocity differ
    by the Earth-rotation term, which is hundreds of metres per second and tilts
    the frame appreciably. Which one the telemetry delivers is not documented, so
    it is resolved by measurement (:mod:`spp.geometry.conventions`) and carried
    explicitly in :class:`~spp.geometry.sensor_model.Convention`.

    Parameters
    ----------
    position, velocity:
        Platform state, arrays of shape ``(..., 3)``, in metres and metres per
        second, both in the same (unspecified) frame.

    Returns
    -------
    numpy.ndarray
        Array of shape ``(..., 3, 3)``.

    Raises
    ------
    ValueError
        If the state is degenerate (zero radius, or velocity parallel to the
        radius vector, which leaves the orbit normal undefined).
    """
    r = np.asarray(position, dtype=np.float64)
    v = np.asarray(velocity, dtype=np.float64)

    r_norm = np.linalg.norm(r, axis=-1, keepdims=True)
    h = np.cross(r, v)
    h_norm = np.linalg.norm(h, axis=-1, keepdims=True)
    if np.any(r_norm < 1e-6) or np.any(h_norm < 1e-6):
        raise ValueError(
            "Degenerate platform state: cannot build an LVLH frame from a zero "
            "radius or a velocity parallel to the radius vector"
        )

    z = -r / r_norm
    y = -h / h_norm
    x = np.cross(y, z)
    return np.stack([x, y, z], axis=-1)


# -- Earth geometry ---------------------------------------------------------


def ray_ellipsoid_intersection(
    origin: np.ndarray,
    direction: np.ndarray,
    height: np.ndarray | float = 0.0,
) -> np.ndarray:
    """Intersect rays with an ellipsoid of constant height above WGS84.

    The surface is the WGS84 ellipsoid with both semi-axes inflated by
    ``height``. That is not rigorously a surface of constant geodetic height,
    but the difference is millimetric at the heights of interest, and the
    terrain intersection iterates anyway (see :mod:`spp.geometry.terrain`).

    Parameters
    ----------
    origin:
        Ray origins in ECEF, shape ``(..., 3)``.
    direction:
        Ray directions in ECEF, shape ``(..., 3)``. Need not be unit length.
    height:
        Height of the target surface above the ellipsoid, in metres. Scalar or
        broadcastable to the ray batch shape.

    Returns
    -------
    numpy.ndarray
        ECEF intersection points, shape ``(..., 3)``. Rays that miss the
        surface yield ``NaN``.

    Notes
    -----
    The **near** intersection (smallest positive root) is returned: the ray from
    a platform enters the Earth at the visible surface.
    """
    o = np.asarray(origin, dtype=np.float64)
    d = np.asarray(direction, dtype=np.float64)
    d = d / np.linalg.norm(d, axis=-1, keepdims=True)

    a = WGS84_A + np.asarray(height, dtype=np.float64)
    b = WGS84_B + np.asarray(height, dtype=np.float64)

    # Scale into a space where the ellipsoid is the unit sphere.
    scale = np.stack([1.0 / a * np.ones_like(o[..., 0]),
                      1.0 / a * np.ones_like(o[..., 1]),
                      1.0 / b * np.ones_like(o[..., 2])], axis=-1)
    os_ = o * scale
    ds = d * scale

    qa = np.sum(ds * ds, axis=-1)
    qb = 2.0 * np.sum(os_ * ds, axis=-1)
    qc = np.sum(os_ * os_, axis=-1) - 1.0

    disc = qb * qb - 4.0 * qa * qc
    real = disc >= 0.0
    sqrt_disc = np.sqrt(np.where(real, disc, 0.0))
    t = (-qb - sqrt_disc) / (2.0 * qa)  # near root

    # A real root is not enough: the *line* through the ray always meets the
    # ellipsoid when the ray points away from it, but the intersection lies
    # behind the platform. Requiring t > 0 is what distinguishes "looking at the
    # Earth" from "looking at the sky".
    hit = real & (t > 0.0)

    point = o + t[..., None] * d
    return np.where(hit[..., None], point, np.nan)


def _rot_x(angle: np.ndarray) -> np.ndarray:
    """Rotation about the x axis, shape ``(..., 3, 3)``."""
    c, s = np.cos(angle), np.sin(angle)
    zero, one = np.zeros_like(c), np.ones_like(c)
    return np.stack(
        [
            np.stack([one, zero, zero], -1),
            np.stack([zero, c, s], -1),
            np.stack([zero, -s, c], -1),
        ],
        axis=-2,
    )


def _rot_y(angle: np.ndarray) -> np.ndarray:
    """Rotation about the y axis, shape ``(..., 3, 3)``."""
    c, s = np.cos(angle), np.sin(angle)
    zero, one = np.zeros_like(c), np.ones_like(c)
    return np.stack(
        [
            np.stack([c, zero, -s], -1),
            np.stack([zero, one, zero], -1),
            np.stack([s, zero, c], -1),
        ],
        axis=-2,
    )


def _rot_z(angle: np.ndarray) -> np.ndarray:
    """Rotation about the z axis, shape ``(..., 3, 3)``."""
    c, s = np.cos(angle), np.sin(angle)
    zero, one = np.zeros_like(c), np.ones_like(c)
    return np.stack(
        [
            np.stack([c, s, zero], -1),
            np.stack([-s, c, zero], -1),
            np.stack([zero, zero, one], -1),
        ],
        axis=-2,
    )


ARCSEC = np.pi / (180.0 * 3600.0)
"""One arcsecond in radians."""

J2000_UNIX = 946728000.0
"""Unix epoch seconds of J2000.0 (2000-01-01T12:00:00 TT, near enough)."""


def _julian_centuries(t: np.ndarray) -> np.ndarray:
    """Julian centuries since J2000.0."""
    return (t - J2000_UNIX) / (86400.0 * 36525.0)


def _precession_matrix(tc: np.ndarray) -> np.ndarray:
    """IAU-76 precession: J2000 mean equator/equinox to mean-of-date.

    **This is not a small correction.** Precession accumulates at about 50
    arcseconds per year, so a quarter-century after J2000 it is of order a third
    of a degree — which at orbital radius is **tens of kilometres** on the ground.
    Neglecting it, as a first pass at this model did, displaces the footprint by
    roughly 40 km: an error large enough to matter and small enough to look like a
    plausible platform bias, which is the worst combination.
    """
    zeta = (2306.2181 * tc + 0.30188 * tc**2 + 0.017998 * tc**3) * ARCSEC
    theta = (2004.3109 * tc - 0.42665 * tc**2 - 0.041833 * tc**3) * ARCSEC
    z = (2306.2181 * tc + 1.09468 * tc**2 + 0.018203 * tc**3) * ARCSEC
    return _rot_z(-z) @ _rot_y(theta) @ _rot_z(-zeta)


def _mean_obliquity(tc: np.ndarray) -> np.ndarray:
    """Mean obliquity of the ecliptic, radians (IAU-76)."""
    return (
        84381.448 - 46.8150 * tc - 0.00059 * tc**2 + 0.001813 * tc**3
    ) * ARCSEC


def _nutation_angles(tc: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Nutation in longitude and obliquity, radians (IAU-80, leading terms).

    Nutation is a ~17 arcsecond wobble — about 500 m at orbital radius. That is
    well below the geolocation error this level carries anyway (§3.2 of the
    spec), so only the dominant terms are kept: the full 106-term series would
    buy metres against an error budget measured in kilometres.
    """
    omega = np.radians(125.04452 - 1934.136261 * tc)  # Moon's ascending node
    lsun = np.radians(280.4665 + 36000.7698 * tc)  # mean longitude of the Sun
    lmoon = np.radians(218.3165 + 481267.8813 * tc)  # mean longitude of the Moon

    d_psi = (
        -17.20 * np.sin(omega)
        - 1.32 * np.sin(2 * lsun)
        - 0.23 * np.sin(2 * lmoon)
        + 0.21 * np.sin(2 * omega)
    ) * ARCSEC
    d_eps = (
        9.20 * np.cos(omega)
        + 0.57 * np.cos(2 * lsun)
        + 0.10 * np.cos(2 * lmoon)
        - 0.09 * np.cos(2 * omega)
    ) * ARCSEC
    return d_psi, d_eps


def _gmst(t: np.ndarray, tc: np.ndarray) -> np.ndarray:
    """Greenwich Mean Sidereal Time, radians."""
    gmst_s = (
        67310.54841
        + (876600.0 * 3600.0 + 8640184.812866) * tc
        + 0.093104 * tc**2
        - 6.2e-6 * tc**3
    )
    return np.radians((gmst_s % 86400.0) / 240.0)  # 240 s of time = 1 degree


def eci_to_ecef_matrix(t: np.ndarray) -> np.ndarray:
    """Rotation from the J2000 inertial frame to ECEF.

    The full chain — **precession, nutation, then Earth rotation**::

        ECEF  <--  R_z(GAST)  <--  nutation  <--  precession  <--  J2000

    Earth rotation alone is not enough, and the gap is not academic. Between J2000
    and the present, precession has turned the celestial reference frame by a
    third of a degree; skipping it puts the platform tens of kilometres from where
    it was. Polar motion (a few metres) is neglected, and the nutation series is
    truncated to its leading terms (a few hundred metres) — both are far below the
    kilometre-scale geolocation error this level carries without a GNSS lock, so
    refining them would be false precision.

    Parameters
    ----------
    t:
        Times, shape ``(m,)``, seconds (Unix epoch). Treated as UT1 ≈ UTC; the
        difference is under a second, i.e. metres of Earth rotation.

    Returns
    -------
    numpy.ndarray
        Array of shape ``(m, 3, 3)``.
    """
    t = np.atleast_1d(np.asarray(t, dtype=np.float64))
    tc = _julian_centuries(t)

    precession = _precession_matrix(tc)

    eps0 = _mean_obliquity(tc)
    d_psi, d_eps = _nutation_angles(tc)
    nutation = _rot_x(-(eps0 + d_eps)) @ _rot_z(-d_psi) @ _rot_x(eps0)

    # Apparent sidereal time: mean, plus the equation of the equinoxes.
    gast = _gmst(t, tc) + d_psi * np.cos(eps0)

    return _rot_z(gast) @ nutation @ precession


def ecef_to_geodetic(points: np.ndarray) -> np.ndarray:
    """Convert ECEF coordinates to geodetic longitude, latitude and height.

    Uses Bowring's method, which converges to sub-millimetre accuracy for
    terrestrial heights without iteration.

    Parameters
    ----------
    points:
        ECEF points, shape ``(..., 3)``, in metres.

    Returns
    -------
    numpy.ndarray
        Array of shape ``(..., 3)`` holding ``(longitude_deg, latitude_deg,
        height_m)``. ``NaN`` inputs propagate.
    """
    p = np.asarray(points, dtype=np.float64)
    x, y, z = p[..., 0], p[..., 1], p[..., 2]

    lon = np.arctan2(y, x)
    rho = np.hypot(x, y)

    # Bowring's auxiliary (parametric) latitude.
    ep2 = (WGS84_A**2 - WGS84_B**2) / WGS84_B**2
    beta = np.arctan2(WGS84_A * z, WGS84_B * rho)
    lat = np.arctan2(
        z + ep2 * WGS84_B * np.sin(beta) ** 3,
        rho - WGS84_E2 * WGS84_A * np.cos(beta) ** 3,
    )

    n = WGS84_A / np.sqrt(1.0 - WGS84_E2 * np.sin(lat) ** 2)
    height = rho / np.cos(lat) - n

    # Near the poles rho -> 0 and the height expression above is ill-conditioned.
    polar = rho < 1.0
    if np.any(polar):
        height = np.where(polar, np.abs(z) - WGS84_B, height)

    return np.stack([np.degrees(lon), np.degrees(lat), height], axis=-1)


def geodetic_to_ecef(lon_lat_h: np.ndarray) -> np.ndarray:
    """Convert geodetic ``(longitude_deg, latitude_deg, height_m)`` to ECEF.

    The exact inverse of :func:`ecef_to_geodetic`.
    """
    g = np.asarray(lon_lat_h, dtype=np.float64)
    lon = np.radians(g[..., 0])
    lat = np.radians(g[..., 1])
    h = g[..., 2]

    n = WGS84_A / np.sqrt(1.0 - WGS84_E2 * np.sin(lat) ** 2)
    return np.stack(
        [
            (n + h) * np.cos(lat) * np.cos(lon),
            (n + h) * np.cos(lat) * np.sin(lon),
            (n * (1.0 - WGS84_E2) + h) * np.sin(lat),
        ],
        axis=-1,
    )
