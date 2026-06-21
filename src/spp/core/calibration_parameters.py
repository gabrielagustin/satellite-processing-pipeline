"""Domain entities describing the radiometric calibration model.

These entities mirror the contents of the per-payload Calibration Parameter
File (``CPF_payload_0.json``). They are pure data holders: the physics (how the
coefficients are combined to turn DN into radiance) lives in the calibration
stage, not here. This keeps the core layer free of processing assumptions.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass(frozen=True)
class RadiometricCalibration:
    """Radiometric calibration coefficients for one ``(band, start_row, tdi)``.

    The three per-pixel coefficient arrays each hold one value per detector
    column. Their physical meaning, as documented in the calibration package,
    is::

        darkfield(T) = thermal_intercept + thermal_gradient * T
        relative_corrected = (DN - darkfield(T)) * non_uniformity
        radiance = absolute * relative_corrected + absolute_offset

    where ``T`` is the per-line detector temperature.

    Attributes
    ----------
    band:
        Band name this entry applies to (e.g. ``"R"``).
    start_row:
        Detector start row this entry applies to. Together with ``band`` and
        ``tdi`` it uniquely identifies the entry for a given acquisition.
    tdi:
        Time-Delay Integration count this entry applies to.
    absolute:
        Absolute calibration scalar (radiance per DN) converting relatively
        corrected counts into at-aperture spectral radiance.
    absolute_offset:
        Absolute calibration additive offset (default ``0.0``).
    non_uniformity:
        Per-column multiplicative correction (detector sensitivity variation).
        Shape ``(width,)``.
    thermal_intercept:
        Per-column, temperature-independent dark-signal term. Shape ``(width,)``.
    thermal_gradient:
        Per-column rate of change of dark signal with temperature.
        Shape ``(width,)``.
    uuid:
        Identifier of the calibration entry, for provenance/traceability.
    """

    band: str
    start_row: int
    tdi: int
    absolute: float
    absolute_offset: float
    non_uniformity: np.ndarray = field(repr=False)
    thermal_intercept: np.ndarray = field(repr=False)
    thermal_gradient: np.ndarray = field(repr=False)
    uuid: str | None = None

    @property
    def width(self) -> int:
        """Number of detector columns the coefficient arrays span."""
        return int(self.non_uniformity.shape[0])


@dataclass(frozen=True)
class CalibrationParameters:
    """Collection of calibration assets loaded from the CPF.

    For the current L1B (radiance) scope only the radiometric entries are
    consumed. The raw ``geometric`` block (boresight alignment + per-band
    line-of-sight arrays) is retained untyped so that future geolocation levels
    (L1C) can build on it without re-reading the file.

    Attributes
    ----------
    radiometric:
        One :class:`RadiometricCalibration` per band present in the CPF.
    geometric:
        Raw geometric calibration block from the CPF, or ``None``. Kept as a
        plain mapping to avoid modelling structure not yet needed.
    """

    radiometric: list[RadiometricCalibration]
    geometric: dict | None = None

    def lookup(self, band: str, start_row: int, tdi: int) -> RadiometricCalibration:
        """Return the radiometric entry matching ``(band, start_row, tdi)``.

        This is the authoritative lookup: the imager configuration of the scene
        provides ``start_row`` and ``tdi`` per band, and the matching entry is
        the one whose ``band``, ``row`` and ``tdi`` all agree.

        Raises
        ------
        KeyError
            If no entry matches the requested tuple.
        """
        for r in self.radiometric:
            if r.band == band and r.start_row == start_row and r.tdi == tdi:
                return r
        raise KeyError(
            f"No radiometric calibration for band={band!r} "
            f"start_row={start_row} tdi={tdi}"
        )
