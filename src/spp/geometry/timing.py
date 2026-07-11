"""Per-line acquisition timing.

A pushbroom raster has no spatial along-track axis: its rows are *instants*. The
sensor model therefore needs, for every raster line, the exact time at which it
was exposed. This module builds that mapping.

The imager stamps each line with its own free-running clock; a synchronisation
record pairs one imager tick with the platform clock, and a calibration offset
corrects the residual skew between them. Line time is then::

    t(line) = platform_epoch_of(anchor)
            + (exposure_ticks[line] - anchor_ticks) * tick_seconds
            + time_sync_offset

Per-line timestamps are **authoritative**. The nominal line period is used only
to *detect* gaps, never to generate times: a dropped line silently compresses
the along-track scale, and a model that regenerates times from a constant period
cannot see that happen.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

TICK_SECONDS = 1e-6
"""Duration of one imager clock tick, seconds (the clock counts microseconds)."""


@dataclass(frozen=True)
class LineTiming:
    """Exposure time of every raster line of an acquisition.

    Attributes
    ----------
    times:
        Exposure time per raster line, shape ``(n_lines,)``, seconds (Unix
        epoch), already corrected by the clock offset.
    nominal_period_s:
        Configured line readout period, seconds. Retained for gap detection and
        for reporting; not used to generate :attr:`times`.
    """

    times: np.ndarray = field(repr=False)
    nominal_period_s: float

    def __post_init__(self) -> None:
        if self.times.ndim != 1 or self.times.size < 2:
            raise ValueError(f"LineTiming needs at least two lines, got {self.times.size}")
        if not np.all(np.diff(self.times) > 0):
            raise ValueError(
                "Line exposure times must be strictly increasing; a non-monotonic "
                "sequence means the timestamps or the scan direction are misread"
            )

    @property
    def n_lines(self) -> int:
        """Number of raster lines."""
        return int(self.times.size)

    @property
    def duration_s(self) -> float:
        """Elapsed time between the first and last line exposure, seconds."""
        return float(self.times[-1] - self.times[0])

    def at(self, lines: np.ndarray | float) -> np.ndarray:
        """Exposure time at (possibly fractional) line indices.

        Fractional indices are linearly interpolated, which is exact wherever
        the readout cadence is uniform and is the sensible reading of "halfway
        between two lines" where it is not.

        Parameters
        ----------
        lines:
            Line indices, scalar or shape ``(m,)``. Zero-based.

        Returns
        -------
        numpy.ndarray
            Times, shape ``(m,)``, seconds (Unix epoch).
        """
        lines = np.atleast_1d(np.asarray(lines, dtype=np.float64))
        return np.interp(lines, np.arange(self.n_lines, dtype=np.float64), self.times)

    def gaps(self, *, tolerance: float = 0.5) -> np.ndarray:
        """Indices of lines preceded by an anomalous readout interval.

        A gap is a step between consecutive exposures that departs from the
        modal step by more than ``tolerance`` of a nominal line period — the
        signature of dropped or duplicated lines.

        Parameters
        ----------
        tolerance:
            Allowed departure, as a fraction of :attr:`nominal_period_s`.

        Returns
        -------
        numpy.ndarray
            Zero-based indices ``i`` such that the interval between line
            ``i - 1`` and line ``i`` is anomalous. Empty when the cadence is
            clean.
        """
        steps = np.diff(self.times)
        modal = float(np.median(steps))
        deviation = np.abs(steps - modal)
        return np.flatnonzero(deviation > tolerance * self.nominal_period_s) + 1


def line_timing_from_exposures(
    exposure_ticks: np.ndarray,
    *,
    anchor_ticks: float,
    anchor_epoch_s: float,
    nominal_period_s: float,
    time_sync_offset_s: float = 0.0,
    tick_seconds: float = TICK_SECONDS,
) -> LineTiming:
    """Build :class:`LineTiming` from raw per-line imager timestamps.

    Parameters
    ----------
    exposure_ticks:
        Imager-clock timestamp of each raster line, shape ``(n_lines,)``, in
        ticks. Must be ordered as the raster rows are stored.
    anchor_ticks:
        Imager-clock tick of the synchronisation record.
    anchor_epoch_s:
        Platform-clock time of that same record, seconds (Unix epoch).
    nominal_period_s:
        Configured line readout period, seconds.
    time_sync_offset_s:
        Calibrated residual offset between the imager and platform clocks,
        seconds. Small (milliseconds) but not negligible: at orbital ground
        speed a millisecond is metres of along-track shift.
    tick_seconds:
        Duration of one imager tick, seconds.

    Returns
    -------
    LineTiming
    """
    ticks = np.asarray(exposure_ticks, dtype=np.float64)
    times = anchor_epoch_s + (ticks - anchor_ticks) * tick_seconds + time_sync_offset_s
    return LineTiming(times=times, nominal_period_s=float(nominal_period_s))
