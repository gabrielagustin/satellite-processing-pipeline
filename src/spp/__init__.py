"""Satellite Processing Pipeline (spp): modular EO processing framework."""

from __future__ import annotations

from importlib import metadata

try:
    __version__ = metadata.version("satellite-processing-pipeline")
except metadata.PackageNotFoundError:  # running from source without an install
    __version__ = "0.0.0+unknown"

__all__ = ["__version__"]
