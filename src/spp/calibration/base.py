"""Abstract calibrator interface."""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np


class Calibrator(ABC):
    """Abstract interface for radiometric calibrators.

    A calibrator converts a block of raw detector values (digital numbers) for a
    single band into calibrated values. It operates on *blocks* rather than whole
    rasters so that large acquisitions can be processed window-by-window without
    loading a full band into memory.

    Calibrators perform pure numerical work: no file I/O and no product
    assembly. Reading the windows, writing the outputs and assembling the
    :class:`~spp.core.product.Product` are the responsibility of the pipeline.
    """

    @abstractmethod
    def calibrate(
        self,
        band_name: str,
        dn: np.ndarray,
        line_start: int = 0,
    ) -> np.ndarray:
        """Calibrate one block of digital numbers.

        Parameters
        ----------
        band_name:
            Name of the band the block belongs to (selects its calibration).
        dn:
            2-D array ``(lines, columns)`` of raw digital numbers for the block.
        line_start:
            Index of the first line of the block within the full band raster.
            Needed because some corrections (e.g. temperature) vary along track.

        Returns
        -------
        numpy.ndarray
            Calibrated block, same shape as ``dn``.
        """
        raise NotImplementedError
