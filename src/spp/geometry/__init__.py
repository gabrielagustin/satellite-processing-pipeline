"""Geometric processing: the sensor model behind the L1C level.

Turns a detector sample into a ground point, by composing the platform's
trajectory and pointing with the instrument's interior orientation. See
``docs/l1c/spec.md`` for the design and the error budget.
"""

from spp.geometry.camera import BandOptics, Camera, CameraIntrinsics
from spp.geometry.conventions import ConventionScore, is_decisive, resolve
from spp.geometry.ephemeris import Attitude, Ephemeris
from spp.geometry.sensor_model import Convention, SensorModel
from spp.geometry.timing import LineTiming

__all__ = [
    "Attitude",
    "BandOptics",
    "Camera",
    "CameraIntrinsics",
    "Convention",
    "ConventionScore",
    "Ephemeris",
    "LineTiming",
    "SensorModel",
    "is_decisive",
    "resolve",
]
