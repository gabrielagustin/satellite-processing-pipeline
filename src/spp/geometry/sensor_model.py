"""The forward sensor model: detector sample to ground point.

This is the heart of the geometric level. It composes the pieces —
:mod:`~spp.geometry.timing` (when was this line exposed),
:mod:`~spp.geometry.ephemeris` (where was the platform, how was it pointing),
:mod:`~spp.geometry.camera` (where was this detector column looking) — into the
one question the level exists to answer::

    (band, line, column)  -->  (longitude, latitude, height)

Undocumented conventions are parameters, not assumptions
--------------------------------------------------------
Several properties of the telemetry are not specified by the delivered package:
which frame the attitude quaternion refers to, in which direction it rotates,
how its components are ordered, and whether raster row zero is the first or the
last line acquired. Guessing wrong does not produce a subtly wrong product — it
produces a footprint on the wrong continent.

So they are not guessed. They are collected in :class:`Convention`, passed in
explicitly, and resolved once by measurement against a known footprint (see
:mod:`spp.geometry.conventions`). Any code reading this model can see exactly
which reading of the telemetry it depends on.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from spp.geometry.camera import Camera
from spp.geometry.ephemeris import Attitude, Ephemeris
from spp.geometry.frames import (
    ecef_to_geodetic,
    eci_to_ecef_matrix,
    lvlh_to_parent_matrix,
    quat_conjugate,
    quat_to_matrix,
    ray_ellipsoid_intersection,
)
from spp.geometry.timing import LineTiming

ColumnAxis = Literal["x", "y"]
AxisSign = Literal[1, -1]
EphemerisFrame = Literal["eci", "ecef"]
AttitudeFrame = Literal["lvlh", "eci", "ecef"]
QuaternionOrder = Literal["scalar_first", "scalar_last"]
QuaternionDirection = Literal["body_to_ref", "ref_to_body"]


@dataclass(frozen=True)
class Convention:
    """How to read the platform telemetry.

    Attributes
    ----------
    ephemeris_frame:
        The frame the position and velocity are expressed in: ``"eci"`` for an
        inertial frame, ``"ecef"`` for the Earth-fixed frame. This is the single
        largest lever in the model — reading an inertial ephemeris as Earth-fixed
        misplaces the platform by the angle the Earth has turned since the epoch,
        which is thousands of kilometres.
    attitude_frame:
        The frame the attitude quaternion relates the body frame to. ``"lvlh"``
        for the local orbital frame (the quaternion is then a small offset from
        nadir pointing), ``"eci"`` for an inertial frame, ``"ecef"`` for the
        Earth-fixed frame.
    quaternion_order:
        Component order as stored: scalar part first or last.
    quaternion_direction:
        Whether the stored quaternion rotates body vectors into the reference
        frame, or the reverse.
    scan_direction:
        ``+1`` if raster row zero is the first line acquired, ``-1`` if it is
        the last (the readout ran against the raster's row order).
    column_axis, column_sign, row_sign:
        How the detector's two axes map into the camera frame — which axis each
        takes, and with what sign (see :meth:`spp.geometry.camera.Camera.rays`).
        Together they decide whether a band's along-track view angle opposes its
        time stagger, as the instrument's design intends, or reinforces it.
    """

    ephemeris_frame: EphemerisFrame = "eci"
    attitude_frame: AttitudeFrame = "lvlh"
    quaternion_order: QuaternionOrder = "scalar_first"
    quaternion_direction: QuaternionDirection = "body_to_ref"
    scan_direction: int = 1
    column_axis: ColumnAxis = "y"
    column_sign: AxisSign = 1
    row_sign: AxisSign = 1

    def describe(self) -> str:
        """One-line human-readable summary, for logs and the decision record."""
        return (
            f"ephemeris_frame={self.ephemeris_frame}, "
            f"attitude_frame={self.attitude_frame}, "
            f"quaternion_order={self.quaternion_order}, "
            f"quaternion_direction={self.quaternion_direction}, "
            f"scan_direction={self.scan_direction:+d}, "
            f"column_axis={self.column_axis}{'+' if self.column_sign > 0 else '-'}, "
            f"row_sign={'+' if self.row_sign > 0 else '-'}"
        )


class SensorModel:
    """Locate detector samples on the ground.

    Parameters
    ----------
    camera:
        Interior orientation (detector column to body-frame ray).
    ephemeris:
        Sampled platform position and velocity, ECEF.
    attitude:
        Sampled platform attitude and angular rates.
    timing:
        Exposure time of every raster line, **per band**. Each band has its own
        timestamp series and they are not interchangeable: the product assembler
        staggers the bands in time to compensate their differing view angles, so
        a band's line times are what tie its raster rows to the platform state.
        Using one band's timing for another silently reintroduces a
        several-hundred-line offset.
    convention:
        How to read the telemetry (see :class:`Convention`).
    rate_aware:
        Use the sampled angular rates when interpolating attitude.
    """

    def __init__(
        self,
        camera: Camera,
        ephemeris: Ephemeris,
        attitude: Attitude,
        timing: dict[str, LineTiming],
        *,
        convention: Convention | None = None,
        rate_aware: bool = True,
    ) -> None:
        self.camera = camera
        self.ephemeris = ephemeris
        self.attitude = attitude
        self.timing = timing
        self.convention = convention or Convention()
        self.rate_aware = rate_aware

    # -- public API ---------------------------------------------------------

    def locate(
        self,
        band: str,
        lines: np.ndarray,
        columns: np.ndarray,
        *,
        height: np.ndarray | float = 0.0,
    ) -> np.ndarray:
        """Ground point observed by each detector sample.

        Parameters
        ----------
        band:
            Band name.
        lines, columns:
            Parallel arrays of raster line and detector column indices, shape
            ``(m,)``. Zero-based; may be fractional.
        height:
            Height above the WGS84 ellipsoid of the surface to intersect, in
            metres. Scalar, or shape ``(m,)``. Terrain-aware location iterates
            this (see ``docs/l1c_spec.md``); with the default of zero this is
            georeferencing on the ellipsoid.

        Returns
        -------
        numpy.ndarray
            Array of shape ``(m, 3)`` holding ``(longitude_deg, latitude_deg,
            height_m)``. Samples whose ray misses the surface yield ``NaN``.
        """
        origins, directions = self.rays_ecef(band, lines, columns)
        ground = ray_ellipsoid_intersection(origins, directions, height=height)
        return ecef_to_geodetic(ground)

    def rays_ecef(
        self, band: str, lines: np.ndarray, columns: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Viewing rays in ECEF: platform positions and unit directions.

        Exposed separately from :meth:`locate` because terrain intersection
        needs to re-intersect the *same* rays at successive height estimates,
        and rebuilding them each iteration would be wasteful.

        Parameters
        ----------
        band:
            Band name.
        lines, columns:
            Parallel arrays, shape ``(m,)``.

        Returns
        -------
        tuple of numpy.ndarray
            ``(origins, directions)``, each of shape ``(m, 3)``, in ECEF.
        """
        lines = np.atleast_1d(np.asarray(lines, dtype=np.float64))
        columns = np.atleast_1d(np.asarray(columns, dtype=np.float64))
        if lines.shape != columns.shape:
            raise ValueError(
                f"lines and columns must be parallel, got {lines.shape} and {columns.shape}"
            )

        timing = self._timing_for(band)
        t = timing.at(self._raster_to_readout(lines, timing))

        # The state arrives in whichever frame the telemetry uses; everything
        # downstream is Earth-fixed, so rotate both the origin and the ray.
        position, velocity = self.ephemeris.interpolate(t)
        eph_to_ecef = self._ephemeris_to_ecef(t)
        origins = np.einsum("mij,mj->mi", eph_to_ecef, position)

        body_to_ecef = self._body_to_ecef(t, position, velocity, eph_to_ecef)
        d_body = self.camera.rays(
            band,
            columns,
            column_axis=self.convention.column_axis,
            column_sign=self.convention.column_sign,
            row_sign=self.convention.row_sign,
        )
        directions = np.einsum("mij,mj->mi", body_to_ecef, d_body)
        return origins, directions

    # -- internals ----------------------------------------------------------

    def _timing_for(self, band: str) -> LineTiming:
        try:
            return self.timing[band]
        except KeyError:
            raise KeyError(
                f"No line timing for band {band!r}; known bands: {sorted(self.timing)}. "
                "Each band carries its own exposure timestamps and they are not "
                "interchangeable."
            ) from None

    def _raster_to_readout(self, lines: np.ndarray, timing: LineTiming) -> np.ndarray:
        """Map raster row indices onto readout order.

        When the scan ran against the raster's row order, row zero holds the
        *last* line acquired, and time must be read from the other end.
        """
        if self.convention.scan_direction >= 0:
            return lines
        return (timing.n_lines - 1) - lines

    def _ephemeris_to_ecef(self, t: np.ndarray) -> np.ndarray:
        """Rotation from the ephemeris frame to ECEF, per sample time."""
        if self.convention.ephemeris_frame == "ecef":
            return np.broadcast_to(np.eye(3), (t.size, 3, 3))
        return eci_to_ecef_matrix(t)

    def _body_to_ecef(
        self,
        t: np.ndarray,
        position: np.ndarray,
        velocity: np.ndarray,
        eph_to_ecef: np.ndarray,
    ) -> np.ndarray:
        """Rotation from the platform body frame to ECEF, per sample time.

        Parameters
        ----------
        t:
            Sample times.
        position, velocity:
            Platform state **in the ephemeris frame** (not yet rotated to ECEF).
            The orbital frame must be built from the state in its own frame: an
            LVLH triad derived from an inertial velocity and one derived from an
            Earth-fixed velocity are genuinely different frames.
        eph_to_ecef:
            Rotation from the ephemeris frame to ECEF, as returned by
            :meth:`_ephemeris_to_ecef`.
        """
        q = self.attitude.interpolate(t, rate_aware=self.rate_aware)
        q = _reorder(q, self.convention.quaternion_order)
        if self.convention.quaternion_direction == "ref_to_body":
            q = quat_conjugate(q)

        body_to_ref = quat_to_matrix(q)

        frame = self.convention.attitude_frame
        if frame == "lvlh":
            # The orbital frame is expressed in the ephemeris frame, so it needs
            # the same rotation into ECEF that the position does.
            lvlh_to_eph = lvlh_to_parent_matrix(position, velocity)
            ref_to_ecef = np.einsum("mij,mjk->mik", eph_to_ecef, lvlh_to_eph)
        elif frame == "eci":
            ref_to_ecef = eci_to_ecef_matrix(t)
        elif frame == "ecef":
            return body_to_ref
        else:  # pragma: no cover - guarded by the Literal type
            raise ValueError(f"Unknown attitude frame: {frame!r}")

        return np.einsum("mij,mjk->mik", ref_to_ecef, body_to_ref)


def _reorder(q: np.ndarray, order: QuaternionOrder) -> np.ndarray:
    """Normalise the stored component order to scalar-first."""
    if order == "scalar_first":
        return q
    return np.roll(q, 1, axis=-1)  # (x, y, z, w) -> (w, x, y, z)
