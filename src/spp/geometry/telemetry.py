"""Adapters: raw package documents to typed geometric entities.

The acquisition package delivers its platform telemetry, imager configuration and
optics as loosely-structured JSON. This module is the single place that knows the
shape of those documents; everything downstream works with the typed entities
(:class:`~spp.geometry.ephemeris.Ephemeris`,
:class:`~spp.geometry.ephemeris.Attitude`,
:class:`~spp.geometry.camera.Camera`,
:class:`~spp.geometry.timing.LineTiming`) and never re-reads a raw dictionary.

The functions here are pure: they take the parsed documents, not file paths, so
they can be tested against synthetic packages without touching a filesystem.
"""

from __future__ import annotations

import logging

import numpy as np

from spp.geometry.camera import BandOptics, Camera, CameraIntrinsics
from spp.geometry.ephemeris import Attitude, Ephemeris
from spp.geometry.timing import LineTiming, line_timing_from_exposures

logger = logging.getLogger(__name__)

POSITION_VELOCITY_FRAME = 3
"""Frame tag the telemetry uses for Earth-fixed position and velocity."""

ATTITUDE_FRAME = 1
"""Frame tag the telemetry uses for attitude and angular rates.

What this frame *is* — inertial, Earth-fixed, or local-orbital — is not stated by
the package. It is resolved by measurement in :mod:`spp.geometry.conventions`.
"""


def _sample_time(entry: dict) -> float:
    """Unix epoch seconds of a telemetry history entry."""
    return float(entry["time_s"]) + float(entry.get("time_ns", 0)) * 1e-9


def _sorted_unique(times: list[float], rows: list[list[float]]) -> tuple[np.ndarray, np.ndarray]:
    """Sort samples by time and drop duplicate timestamps.

    Telemetry histories are not guaranteed to arrive ordered, and a repeated
    timestamp would make interpolation ill-posed (a zero-length segment).
    """
    t = np.asarray(times, dtype=np.float64)
    v = np.asarray(rows, dtype=np.float64)
    order = np.argsort(t, kind="stable")
    t, v = t[order], v[order]
    keep = np.concatenate([[True], np.diff(t) > 0])
    if not keep.all():
        logger.warning(
            "Telemetry: dropped %d sample(s) with duplicate timestamps", int((~keep).sum())
        )
    return t[keep], v[keep]


def ephemeris_from_ancillary(ancillary: dict) -> Ephemeris:
    """Build an :class:`Ephemeris` from the platform telemetry history.

    Position and velocity are delivered as separate entries on their own time
    grids. They are parsed independently and the velocity is resampled onto the
    position times, so the two are consistent by construction.

    Raises
    ------
    ValueError
        If the history carries no position or no velocity samples.
    """
    hist = ancillary.get("extrinsics", {}).get("hist", [])

    p_times, p_vals = [], []
    v_times, v_vals = [], []
    for entry in hist:
        if "position" in entry:
            d = entry["position"]
            p_times.append(_sample_time(entry))
            p_vals.append([d["x_m"], d["y_m"], d["z_m"]])
        elif "velocity" in entry:
            d = entry["velocity"]
            v_times.append(_sample_time(entry))
            v_vals.append([d["x_m"], d["y_m"], d["z_m"]])

    if not p_times or not v_times:
        raise ValueError(
            "Platform telemetry carries no position and/or velocity history; "
            "the geometric level cannot proceed without an ephemeris"
        )

    pt, positions = _sorted_unique(p_times, p_vals)
    vt, velocities = _sorted_unique(v_times, v_vals)

    if not np.array_equal(pt, vt):
        velocities = np.stack(
            [np.interp(pt, vt, velocities[:, i]) for i in range(3)], axis=-1
        )
    return Ephemeris(times=pt, positions=positions, velocities=velocities)


def attitude_from_ancillary(ancillary: dict) -> Attitude:
    """Build an :class:`Attitude` from the platform telemetry history.

    The quaternion is stored with keys ``w, q0, q1, q2``. It is read as
    **scalar-first** ``[w, q0, q1, q2]`` — which the data supports (the four
    components are unit-norm together) — but whether that reading, and the
    rotation's direction and reference frame, are correct is settled empirically
    by :mod:`spp.geometry.conventions`, not asserted here.

    Angular rates are resampled onto the attitude times when the two histories
    are sampled separately.
    """
    hist = ancillary.get("extrinsics", {}).get("hist", [])

    q_times, q_vals = [], []
    r_times, r_vals = [], []
    for entry in hist:
        if "attitude" in entry:
            d = entry["attitude"]
            q_times.append(_sample_time(entry))
            q_vals.append([d["w"], d["q0"], d["q1"], d["q2"]])
        elif "angular_rates" in entry:
            d = entry["angular_rates"]
            r_times.append(_sample_time(entry))
            r_vals.append([d["x_rads"], d["y_rads"], d["z_rads"]])

    if not q_times:
        raise ValueError(
            "Platform telemetry carries no attitude history; the geometric level "
            "cannot proceed without it"
        )

    qt, quaternions = _sorted_unique(q_times, q_vals)

    rates = None
    if r_times:
        rt, r = _sorted_unique(r_times, r_vals)
        rates = (
            r
            if np.array_equal(qt, rt)
            else np.stack([np.interp(qt, rt, r[:, i]) for i in range(3)], axis=-1)
        )
    else:
        logger.warning(
            "Platform telemetry carries no angular rates; attitude interpolation "
            "falls back to SLERP"
        )

    return Attitude(times=qt, quaternions=quaternions, rates=rates)


def gnss_lock(ancillary: dict) -> bool:
    """Whether the platform had a GNSS fix during the acquisition.

    ``False`` means the positions are propagated rather than measured, and the
    product's absolute geolocation cannot be trusted without external refinement.
    """
    return bool(ancillary.get("extrinsics", {}).get("gnss_lock", False))


def camera_from_package(
    ancillary: dict,
    geometric_calibration: dict | None,
    band_start_rows: dict[str, int],
    band_tdi: dict[str, int],
) -> Camera:
    """Build the :class:`Camera` from the intrinsics and geometric calibration.

    Parameters
    ----------
    ancillary:
        Platform ancillary document (supplies ``intrinsics``).
    geometric_calibration:
        The ``geometric`` block of the Calibration Parameter File, or ``None``.
        Supplies the boresight matrix and the per-band line-of-sight correction
        coefficients — both unpopulated in the reference package, which is why
        the interior orientation is reconstructed from the intrinsics and the
        detector rows instead (see :mod:`spp.geometry.camera`).
    band_start_rows, band_tdi:
        Per-band detector start row and Time-Delay Integration count, from the
        imager configuration.
    """
    intr = ancillary["intrinsics"]
    intrinsics = CameraIntrinsics(
        focal_length_mm=float(intr["focal_length_mm"]),
        pixel_size_mm=float(intr["pixel_size_mm"]),
        principal_point_px=tuple(float(v) for v in intr["principal_point_px"]),
        sensor_size_px=tuple(int(v) for v in intr["sensor_size_px"]),
        detector_to_body=np.asarray(intr["reference_frame"], dtype=np.float64),
    )

    geometric_calibration = geometric_calibration or {}
    los = {
        entry["band"]: entry
        for entry in geometric_calibration.get("lineOfSight", [])
        if "band" in entry
    }

    bands = {
        name: BandOptics(
            name=name,
            start_row=int(start_row),
            tdi=int(band_tdi.get(name, 1)),
            los_along=tuple(los.get(name, {}).get("along", ()) or ()),
            los_across=tuple(los.get(name, {}).get("across", ()) or ()),
        )
        for name, start_row in band_start_rows.items()
    }

    boresight = geometric_calibration.get("boreSightAlignment")
    return Camera(
        intrinsics,
        bands,
        boresight=np.asarray(boresight, dtype=np.float64) if boresight else None,
    )


def line_timing_from_session(
    session: dict,
    band_name_to_id: dict[str, int],
    *,
    time_sync_offset_s: float = 0.0,
    scene_index: str = "0",
) -> dict[str, LineTiming]:
    """Build per-band :class:`LineTiming` from the session metadata.

    **Every band carries its own exposure timestamps**, offset from its
    neighbours' by the time the platform takes to fly the gap between their
    detector rows. That stagger is how the product compensates the bands'
    differing view angles, so each band's timing must be read separately: reusing
    one band's series for another reintroduces a several-hundred-line offset.

    Parameters
    ----------
    session:
        The session document for this acquisition (the value under
        ``Sessions/<id>``).
    band_name_to_id:
        Band name to detector band id, as the raw band streams are keyed.
    time_sync_offset_s:
        Calibrated offset between the imager and platform clocks, seconds.
    scene_index:
        Which scene of the session to read.

    Returns
    -------
    dict
        Band name to :class:`LineTiming`. Bands absent from the session metadata
        are skipped with a warning rather than failing the run.

    Raises
    ------
    ValueError
        If the clock synchronisation anchor is missing, without which the imager
        clock cannot be placed on the platform timeline.
    """
    anchor = _time_sync_anchor(session)
    if anchor is None:
        raise ValueError(
            "Session metadata carries no usable TimeSync record (an ImagerTime "
            "paired with a PlatformTime); line times cannot be placed on the "
            "platform timeline"
        )
    anchor_ticks, anchor_epoch_s = anchor

    line_period_s = float(session["ImagerConfiguration"]["LinePeriod"]) * 1e-6
    raw_bands = session["Scenes"][scene_index]["RawBands"]

    timings: dict[str, LineTiming] = {}
    for name, band_id in band_name_to_id.items():
        stream = raw_bands.get(str(band_id))
        if stream is None or "Lines" not in stream:
            logger.warning("Band %s (id %d): no line timing in session metadata", name, band_id)
            continue
        lines = stream["Lines"]
        ticks = np.array(
            [lines[k]["ExposureTimestamp"] for k in sorted(lines, key=int)],
            dtype=np.float64,
        )
        timings[name] = line_timing_from_exposures(
            ticks,
            anchor_ticks=anchor_ticks,
            anchor_epoch_s=anchor_epoch_s,
            nominal_period_s=line_period_s,
            time_sync_offset_s=time_sync_offset_s,
        )
    return timings


def _time_sync_anchor(session: dict) -> tuple[float, float] | None:
    """First synchronisation record pairing the imager and platform clocks.

    Most ``TimeSync`` entries are bare pulse-per-second ticks with no platform
    time; only the records that carry both clocks can anchor one to the other.
    """
    for entry in session.get("TimeSync", []):
        if "ImagerTime" in entry and "PlatformTime" in entry:
            return float(entry["ImagerTime"]), float(entry["PlatformTime"]) / 1000.0
    return None
