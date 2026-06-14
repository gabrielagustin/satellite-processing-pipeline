from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass
class Acquisition:
    """
    Represents a raw satellite acquisition.

    Attributes
    ----------
    image : np.ndarray
        Raw detector measurements.

    metadata : dict[str, Any]
        Acquisition metadata.
    """

    image: np.ndarray
    metadata: dict[str, Any]
