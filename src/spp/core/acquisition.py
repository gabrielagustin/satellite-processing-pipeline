"""Domain entities describing a full image acquisition package."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from spp.core.band import Band
from spp.core.calibration_parameters import CalibrationParameters


@dataclass(frozen=True)
class ImagerConfiguration:
    """Per-acquisition imager configuration (from ``metadata.json``).

    The per-band arrays are indexed by detector band id (0..N-1). They are the
    bridge between a band and its calibration entry: the ``(start_row, tdi)``
    tuple for a band selects the matching radiometric calibration.

    Attributes
    ----------
    line_period_us:
        Detector line readout period in microseconds.
    spectral_bands:
        Number of spectral bands.
    band_setup:
        TDI count per band id.
    band_start_row:
        Detector start row per band id.
    band_cwl:
        Central wavelength (nm) per band id.
    scan_direction:
        Along-track scan direction (``1`` = forward, ``-1`` = reverse).
    binning_factor:
        On-detector binning factor (``0`` = none).
    """

    line_period_us: int
    spectral_bands: int
    band_setup: list[int]
    band_start_row: list[int]
    band_cwl: list[float]
    scan_direction: int
    binning_factor: int

    def tdi_for(self, band_id: int) -> int:
        """TDI count configured for ``band_id``."""
        return self.band_setup[band_id]

    def start_row_for(self, band_id: int) -> int:
        """Detector start row configured for ``band_id``."""
        return self.band_start_row[band_id]


@dataclass(frozen=True)
class Acquisition:
    """An image acquisition package ready for processing.

    This is the common data structure produced by a reader and consumed by the
    downstream stages. It bundles the per-band rasters (lazily referenced), the
    calibration model, the imager configuration and the ancillary telemetry
    needed by calibration — without holding any pixel data in memory.

    Attributes
    ----------
    scene_id:
        Product/scene identifier (the L0 STAC item id).
    bands:
        Mapping of band name to :class:`Band`.
    calibration:
        Calibration model for this acquisition.
    imager_config:
        Imager configuration for this acquisition.
    sensor_temperatures:
        Detector temperature telemetry samples (degrees, as stored). Used to
        evaluate the temperature-dependent darkfield. May be coarsely sampled
        (one value per telemetry tick, not per line).
    temperature_sample_lines:
        Along-track line index (fractional) of each detector-temperature sample,
        derived from the sample timestamps and the per-line timing. Parallel to
        :attr:`sensor_temperatures`. ``None`` when timing information is
        unavailable, in which case downstream stages fall back to spreading the
        samples uniformly across the lines.
    acquired_at:
        Acquisition datetime in ISO-8601 UTC, if known.
    metadata:
        Raw session metadata mapping, retained for provenance and for fields
        not promoted to typed attributes.
    ancillary:
        Raw platform ancillary mapping (position/velocity/attitude history),
        retained for future geolocation levels.
    """

    scene_id: str
    bands: dict[str, Band]
    calibration: CalibrationParameters
    imager_config: ImagerConfiguration
    sensor_temperatures: np.ndarray = field(repr=False)
    temperature_sample_lines: np.ndarray | None = field(default=None, repr=False)
    acquired_at: str | None = None
    metadata: dict = field(default_factory=dict, repr=False)
    ancillary: dict = field(default_factory=dict, repr=False)

    @property
    def band_names(self) -> list[str]:
        """Band names available in this acquisition."""
        return list(self.bands.keys())

    def band(self, name: str) -> Band:
        """Return the :class:`Band` named ``name``."""
        return self.bands[name]
