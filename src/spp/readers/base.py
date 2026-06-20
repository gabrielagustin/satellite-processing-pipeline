"""Abstract reader interface."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from spp.core.acquisition import Acquisition


class Reader(ABC):
    """Abstract interface for acquisition readers.

    A reader is responsible for turning the contents of an acquisition package
    (on the filesystem, in object storage, ...) into a common
    :class:`~spp.core.acquisition.Acquisition` domain object.

    Readers perform **no** calibration, processing or quality assessment. Their
    sole responsibility is discovery and parsing.
    """

    @abstractmethod
    def read(self) -> "Acquisition":
        """Read the package and return an :class:`Acquisition`."""
        raise NotImplementedError
