"""Radiometric calibration to L1B at-aperture spectral radiance.

The calibrator turns raw detector counts (DN) into top-of-atmosphere spectral
radiance, applying — in order — a temperature-dependent dark-signal subtraction,
a per-pixel non-uniformity (relative) correction, and the absolute calibration
scaling::

    darkfield(T)       = thermal_intercept[col] + thermal_gradient[col] * T
    relative_corrected = (DN - darkfield(T)) * non_uniformity[col]
    radiance           = absolute * relative_corrected + absolute_offset

The detector temperature ``T`` varies along track. It is provided by a small
number of telemetry samples, which are interpolated to a per-line temperature
profile (see :meth:`per_line_temperature`).
"""

from __future__ import annotations

import numpy as np

from spp.calibration.base import Calibrator
from spp.core.acquisition import Acquisition


class L1BCalibrator(Calibrator):
    """Convert raw DN blocks into L1B at-aperture spectral radiance.

    Parameters
    ----------
    acquisition:
        The acquisition providing the calibration model, imager configuration
        and detector temperature telemetry.
    nodata:
        Value written where the input is NoData. Defaults to NaN.
    dtype:
        Output dtype. Defaults to ``float32`` (radiance fits comfortably and it
        halves memory versus ``float64``).

    Notes
    -----
    Radiance is **not** clipped: pixels near the dark level may yield slightly
    negative values after dark subtraction. These are preserved so that quality
    assessment can flag them rather than silently masking calibration bias.
    """

    #: Physical units of the output radiance (TOA spectral radiance, per micrometre).
    UNITS = "W / (m2 sr um)"

    def __init__(
        self,
        acquisition: Acquisition,
        *,
        nodata: float = float("nan"),
        dtype: np.dtype | str = np.float32,
    ) -> None:
        self._acq = acquisition
        self._nodata = np.float32(nodata)
        self._dtype = np.dtype(dtype)
        self._context: dict[str, _BandContext] = {}

    @property
    def units(self) -> str:
        """Physical units of the calibrated radiance."""
        return self.UNITS

    @property
    def nodata(self) -> float:
        """NoData value written for masked pixels."""
        return float(self._nodata)

    def calibrate(
        self,
        band_name: str,
        dn: np.ndarray,
        line_start: int = 0,
    ) -> np.ndarray:
        """Calibrate one block of DN for ``band_name`` into radiance."""
        ctx = self._band_context(band_name)
        if dn.ndim != 2:
            raise ValueError(f"Expected a 2-D DN block, got shape {dn.shape}")
        h, w = dn.shape
        if w != ctx.width:
            raise ValueError(
                f"Band {band_name!r}: block width {w} != calibration width {ctx.width}"
            )
        line_stop = line_start + h
        if line_stop > ctx.temperature.size:
            raise ValueError(
                f"Block lines [{line_start}:{line_stop}] exceed band height "
                f"{ctx.temperature.size}"
            )

        temperature = ctx.temperature[line_start:line_stop]  # (h,)
        dn_f = dn.astype(np.float32, copy=False)

        # darkfield(line, col) = intercept[col] + gradient[col] * T[line]
        darkfield = ctx.thermal_intercept[None, :] + (
            ctx.thermal_gradient[None, :] * temperature[:, None]
        )
        relative = (dn_f - darkfield) * ctx.non_uniformity[None, :]
        radiance = ctx.absolute * relative + ctx.absolute_offset

        nodata = self._acq.band(band_name).nodata
        if nodata is not None:
            radiance[dn == nodata] = self._nodata

        return radiance.astype(self._dtype, copy=False)

    def per_line_temperature(self, band_name: str) -> np.ndarray:
        """Per-line detector temperature profile for ``band_name``."""
        return self._band_context(band_name).temperature

    # -- internals ----------------------------------------------------------

    def _band_context(self, band_name: str) -> "_BandContext":
        ctx = self._context.get(band_name)
        if ctx is None:
            ctx = self._build_context(band_name)
            self._context[band_name] = ctx
        return ctx

    def _build_context(self, band_name: str) -> "_BandContext":
        band = self._acq.band(band_name)
        cfg = self._acq.imager_config
        rc = self._acq.calibration.lookup(
            band_name,
            cfg.start_row_for(band.band_id),
            cfg.tdi_for(band.band_id),
        )
        temperature = self._interpolate_temperature(band.height)
        return _BandContext(
            width=rc.width,
            non_uniformity=rc.non_uniformity.astype(np.float32),
            thermal_intercept=rc.thermal_intercept.astype(np.float32),
            thermal_gradient=rc.thermal_gradient.astype(np.float32),
            absolute=np.float32(rc.absolute),
            absolute_offset=np.float32(rc.absolute_offset),
            temperature=temperature.astype(np.float32),
        )

    def _interpolate_temperature(self, n_lines: int) -> np.ndarray:
        """Interpolate the telemetry temperatures onto ``n_lines`` lines.

        When the reader has resolved the along-track line position of each
        telemetry sample from its timestamp (``temperature_sample_lines``), the
        samples are interpolated **at their true line positions** — exact even if
        the telemetry is irregularly spaced or brackets the acquisition window.

        When that timing is unavailable, the samples are spread **uniformly**
        across the lines. This is a first-order approximation (it assumes the
        samples are equally spaced in time and exactly span the acquisition) and
        is documented as a known limitation.
        """
        temps = np.asarray(self._acq.sensor_temperatures, dtype=np.float64)
        if temps.size == 0:
            raise ValueError("No detector temperature samples available")
        if temps.size == 1:
            return np.full(n_lines, temps[0], dtype=np.float64)

        sample_positions = self._acq.temperature_sample_lines
        if sample_positions is None:
            # Fallback: assume samples are uniformly spaced over the lines.
            sample_positions = np.linspace(0.0, n_lines - 1, temps.size)
        else:
            sample_positions = np.asarray(sample_positions, dtype=np.float64)

        # np.interp needs ascending sample positions; sort jointly to be safe.
        order = np.argsort(sample_positions)
        line_positions = np.arange(n_lines, dtype=np.float64)
        return np.interp(
            line_positions, sample_positions[order], temps[order]
        )


class _BandContext:
    """Precomputed, float32 calibration context for a single band."""

    __slots__ = (
        "width",
        "non_uniformity",
        "thermal_intercept",
        "thermal_gradient",
        "absolute",
        "absolute_offset",
        "temperature",
    )

    def __init__(
        self,
        width: int,
        non_uniformity: np.ndarray,
        thermal_intercept: np.ndarray,
        thermal_gradient: np.ndarray,
        absolute: np.float32,
        absolute_offset: np.float32,
        temperature: np.ndarray,
    ) -> None:
        self.width = width
        self.non_uniformity = non_uniformity
        self.thermal_intercept = thermal_intercept
        self.thermal_gradient = thermal_gradient
        self.absolute = absolute
        self.absolute_offset = absolute_offset
        self.temperature = temperature
