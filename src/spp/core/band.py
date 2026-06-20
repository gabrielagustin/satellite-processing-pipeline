"""Domain entity describing a single spectral band of an acquisition."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Band:
    """A single spectral band within an acquisition package.

    The entity is deliberately *lazy*: it carries the location and the raster
    properties of the band image, but never the pixel data itself. Pushbroom
    band rasters can be very large (thousands of columns by tens of thousands of
    lines, hundreds of MB each), so loading every band into memory at once is
    wasteful. Downstream stages read the pixels in windows directly from
    :attr:`path` using their I/O library of choice, which keeps this core entity
    free of any raster-library dependency.

    Attributes
    ----------
    name:
        Short band name as used in the package (e.g. ``"R"``, ``"NIR"``).
    band_id:
        Zero-based detector band id. Used to index the per-band arrays in the
        imager configuration (``BandSetup``, ``BandStartRow`` ...).
    cwl_nm:
        Central wavelength of the band in nanometres.
    path:
        Filesystem path to the band raster (a GeoTIFF in the L0 package).
    width:
        Raster width in pixels (cross-track / detector column axis).
    height:
        Raster height in pixels (along-track / acquisition line axis).
    dtype:
        NumPy dtype name of the stored samples (e.g. ``"uint16"``).
    nodata:
        NoData value of the raster, or ``None`` if undefined. A common
        convention is ``0`` to flag detector pixels that read out as no-signal.
    """

    name: str
    band_id: int
    cwl_nm: float
    path: Path
    width: int
    height: int
    dtype: str
    nodata: float | None = None

    @property
    def shape(self) -> tuple[int, int]:
        """Raster shape as ``(height, width)`` — i.e. ``(lines, columns)``."""
        return (self.height, self.width)
