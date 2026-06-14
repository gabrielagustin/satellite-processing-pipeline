from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass
class Product:
    """
    Represents a generated processing product.
    """

    level: str
    image: np.ndarray
    metadata: dict[str, Any]
