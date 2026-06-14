from dataclasses import dataclass


@dataclass
class CalibrationModel:
    """
    Radiometric calibration parameters.
    """

    gain: float
    offset: float
