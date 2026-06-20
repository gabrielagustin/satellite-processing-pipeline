"""Abstract product writer interface."""

from __future__ import annotations

from abc import ABC, abstractmethod
from contextlib import AbstractContextManager

import numpy as np


class ProductWriter(ABC):
    """Abstract interface for raster product writers.

    Writers persist a calibrated raster product to storage. To support the
    windowed (streaming) processing of large rasters, a writer exposes a band as
    a context-managed handle that accepts blocks of lines, rather than a single
    whole-array call::

        with writer.open_band("R", width=W, height=H, dtype="float32",
                              nodata=float("nan")) as band:
            for block, line_start in blocks:
                band.write_block(block, line_start)

    The handle yielded by :meth:`open_band` must provide:

    - ``write_block(array, line_start)`` — write a ``(lines, width)`` block whose
      first row maps to ``line_start`` in the full raster;
    - ``path`` — the filesystem path the band was written to.

    Writers contain no processing logic. Catalogue/metadata persistence (e.g.
    STAC) is a separate concern handled by the metadata layer.
    """

    @abstractmethod
    def open_band(
        self,
        name: str,
        *,
        width: int,
        height: int,
        dtype: np.dtype | str,
        nodata: float | None,
    ) -> AbstractContextManager:
        """Open a band for windowed writing; see class docstring."""
        raise NotImplementedError
