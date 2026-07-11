"""Camera model: detector sample to viewing ray.

This is the *interior orientation* of the instrument — the part of the geometry
that depends on the optics and the detector, not on where the platform was or
how it was pointing.

Built from first principles
---------------------------
The delivered geometric calibration is a **placeholder**: an identity boresight
matrix, and line-of-sight coefficients that are identical for all bands. That
cannot be physically true — the bands sit on different detector rows and
therefore demonstrably look in different directions. So the interior orientation
is reconstructed from the pinhole model implied by the intrinsics (focal length,
pixel pitch, principal point) and each band's detector row, and the calibration's
line-of-sight term enters as an **additive angular correction** layered on top —
zero as delivered, and estimable from the imagery itself (see
``docs/l1c_spec.md``, band co-registration).

That per-band row offset is not a detail. It is what makes each band look
slightly forward or backward of the others, so it is simultaneously the cause of
the band-to-band ground offset and, once modelled, the cure.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np

logger = logging.getLogger(__name__)

UNPOPULATED_LOS_SENTINEL = (1.0, 0.0, 0.0, 0.0, 0.0, 0.0)
"""Line-of-sight coefficients that mean "not calibrated".

Read as an additive angular polynomial, these coefficients would be a constant
offset of **one radian** — 57 degrees, an absurd correction for a boresighted
instrument, and identical for every band besides. They are an unpopulated
default, and are treated as a zero correction (with a warning), never applied.
"""


@dataclass(frozen=True)
class CameraIntrinsics:
    """Interior geometry of the imager.

    Attributes
    ----------
    focal_length_mm:
        Effective focal length, millimetres.
    pixel_size_mm:
        Detector pixel pitch, millimetres.
    principal_point_px:
        Principal point ``(column, row)`` in detector pixels.
    sensor_size_px:
        Detector size ``(columns, rows)`` in pixels.
    detector_to_body:
        3x3 rotation from the detector/camera frame to the platform body frame.
    """

    focal_length_mm: float
    pixel_size_mm: float
    principal_point_px: tuple[float, float]
    sensor_size_px: tuple[int, int]
    detector_to_body: np.ndarray = field(repr=False)

    def __post_init__(self) -> None:
        if self.focal_length_mm <= 0 or self.pixel_size_mm <= 0:
            raise ValueError("Focal length and pixel pitch must be positive")
        m = np.asarray(self.detector_to_body, dtype=np.float64)
        if m.shape != (3, 3):
            raise ValueError(f"detector_to_body must be 3x3, got {m.shape}")
        object.__setattr__(self, "detector_to_body", m)

    @property
    def half_field_rad(self) -> float:
        """Cross-track half-field of view, radians.

        A narrow field is what makes terrain displacement small in absolute
        terms — and what makes the *band-to-band* view-angle spread, small as it
        looks, a comparable fraction of the field.
        """
        half_width_mm = 0.5 * self.sensor_size_px[0] * self.pixel_size_mm
        return float(np.arctan(half_width_mm / self.focal_length_mm))


@dataclass(frozen=True)
class BandOptics:
    """Where one spectral band sits on the detector, and how it is corrected.

    Attributes
    ----------
    name:
        Band name (e.g. ``"NIR"``).
    start_row:
        Detector row at which the band's readout window begins.
    tdi:
        Time-Delay Integration count: the number of detector rows accumulated
        into each output line. The band's effective row is the centre of that
        window, not its first row.
    los_along, los_across:
        Additive angular correction coefficients, a polynomial in normalised
        detector column, in radians. Empty (or the unpopulated sentinel) means
        no correction.
    """

    name: str
    start_row: int
    tdi: int = 1
    los_along: tuple[float, ...] = ()
    los_across: tuple[float, ...] = ()

    @property
    def effective_row(self) -> float:
        """Detector row the band effectively images along.

        With Time-Delay Integration the output line is the sum of ``tdi``
        consecutive detector rows, so the effective row is the centre of that
        window.
        """
        return self.start_row + 0.5 * (self.tdi - 1)


class Camera:
    """Turns detector columns into unit viewing rays in the body frame.

    Parameters
    ----------
    intrinsics:
        Interior geometry of the imager.
    bands:
        Per-band detector placement and line-of-sight correction, keyed by band
        name.
    boresight:
        3x3 boresight alignment matrix (body frame). Identity if uncalibrated.
    """

    def __init__(
        self,
        intrinsics: CameraIntrinsics,
        bands: dict[str, BandOptics],
        *,
        boresight: np.ndarray | None = None,
    ) -> None:
        self.intrinsics = intrinsics
        self.bands = bands
        self.boresight = (
            np.eye(3) if boresight is None else np.asarray(boresight, dtype=np.float64)
        )
        if self.boresight.shape != (3, 3):
            raise ValueError(f"boresight must be 3x3, got {self.boresight.shape}")

    # -- public API ---------------------------------------------------------

    def view_angles(self, band: str, columns: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Across-track and along-track view angles of detector columns.

        Parameters
        ----------
        band:
            Band name.
        columns:
            Detector column indices, shape ``(m,)``. Zero-based, may be
            fractional.

        Returns
        -------
        tuple of numpy.ndarray
            ``(across_rad, along_rad)``, each of shape ``(m,)``. The across-track
            angle sweeps with the column; the along-track angle is **constant per
            band** (it is the band's detector row offset) and is what separates
            the bands on the ground.
        """
        optics = self._optics(band)
        cols = np.atleast_1d(np.asarray(columns, dtype=np.float64))

        c0, r0 = self.intrinsics.principal_point_px
        p = self.intrinsics.pixel_size_mm
        f = self.intrinsics.focal_length_mm

        across = np.arctan((cols - c0) * p / f)
        along = np.full_like(across, np.arctan((optics.effective_row - r0) * p / f))

        u = self._normalised_column(cols)
        across = across + _polyval(_usable_coefficients(optics.los_across, band, "across"), u)
        along = along + _polyval(_usable_coefficients(optics.los_along, band, "along"), u)
        return across, along

    def rays(
        self,
        band: str,
        columns: np.ndarray,
        *,
        column_axis: str = "y",
        column_sign: int = 1,
        row_sign: int = 1,
    ) -> np.ndarray:
        """Unit viewing rays for detector columns, in the platform body frame.

        Parameters
        ----------
        band:
            Band name.
        columns:
            Detector column indices, shape ``(m,)``. Zero-based, may be
            fractional.
        column_axis:
            Which camera axis the detector's **column** direction lies along:
            ``"x"`` or ``"y"``. The detector's row direction takes the other.

            This is not a free choice, it is an undocumented one. The package
            gives a detector-to-body matrix but never says which way round the
            detector's two axes feed into it, and the two readings are not
            interchangeable: one puts the band's along-track view angle on the
            platform's along-track axis, where the per-band time stagger cancels
            it, and the other puts it across-track, where nothing cancels it and
            the bands land kilometres apart. It is resolved by measurement in
            :mod:`spp.geometry.conventions`.

        Returns
        -------
        numpy.ndarray
            Unit vectors, shape ``(m, 3)``, in the body frame.
        """
        across, along = self.view_angles(band, columns)
        ones = np.ones_like(across)

        col_component = column_sign * np.tan(across)
        row_component = row_sign * np.tan(along)

        # Camera frame: +z along the optical axis; the two detector axes take x
        # and y, in the order `column_axis` selects.
        if column_axis == "x":
            d_cam = np.stack([col_component, row_component, ones], axis=-1)
        elif column_axis == "y":
            d_cam = np.stack([row_component, col_component, ones], axis=-1)
        else:
            raise ValueError(f"column_axis must be 'x' or 'y', got {column_axis!r}")
        d_cam /= np.linalg.norm(d_cam, axis=-1, keepdims=True)

        d_body = d_cam @ self.intrinsics.detector_to_body.T
        d_body = d_body @ self.boresight.T
        return d_body / np.linalg.norm(d_body, axis=-1, keepdims=True)

    def effective_los(self, band: str) -> tuple[tuple[float, ...], tuple[float, ...]]:
        """The line-of-sight coefficients actually in force for a band.

        The *stored* coefficients may be the unpopulated sentinel
        (:data:`UNPOPULATED_LOS_SENTINEL`), which is not a correction and is never
        applied. Anything building on the current calibration — the self-calibration,
        above all — must start from what is **in force**, not from what is stored.

        Reading the raw coefficients instead would treat the sentinel's leading ``1.0``
        as a one-radian offset (57 degrees) and inherit it as the baseline for a
        correction, producing a confidently-computed, catastrophically wrong result.
        """
        optics = self._optics(band)
        return (
            _usable_coefficients(optics.los_along, band, "along"),
            _usable_coefficients(optics.los_across, band, "across"),
        )

    def with_boresight(self, boresight: np.ndarray) -> Camera:
        """Return a copy of this camera with a different boresight matrix.

        Used by the absolute refinement, which corrects pointing by adjusting
        the boresight rather than by shifting the output image.
        """
        return Camera(self.intrinsics, self.bands, boresight=boresight)

    def with_los(self, corrections: dict[str, tuple[tuple[float, ...], tuple[float, ...]]]) -> Camera:
        """Return a copy with per-band line-of-sight corrections replaced.

        Parameters
        ----------
        corrections:
            Mapping of band name to ``(los_along, los_across)`` coefficient
            tuples, in radians. This is the seam the band co-registration writes
            back into: the calibration the instrument did not ship, estimated
            from the imagery.
        """
        bands = dict(self.bands)
        for name, (along, across) in corrections.items():
            optics = self.bands[name]
            bands[name] = BandOptics(
                name=optics.name,
                start_row=optics.start_row,
                tdi=optics.tdi,
                los_along=tuple(along),
                los_across=tuple(across),
            )
        return Camera(self.intrinsics, bands, boresight=self.boresight)

    # -- internals ----------------------------------------------------------

    def _optics(self, band: str) -> BandOptics:
        try:
            return self.bands[band]
        except KeyError:
            raise KeyError(
                f"No optics for band {band!r}; known bands: {sorted(self.bands)}"
            ) from None

    def _normalised_column(self, columns: np.ndarray) -> np.ndarray:
        """Map detector columns onto ``[-1, 1]`` about the principal point."""
        c0 = self.intrinsics.principal_point_px[0]
        half_width = 0.5 * self.intrinsics.sensor_size_px[0]
        return (columns - c0) / half_width


def _polyval(coefficients: tuple[float, ...], u: np.ndarray) -> np.ndarray:
    """Evaluate ``sum(c_k * u**k)``. Empty coefficients mean a zero correction."""
    if not coefficients:
        return np.zeros_like(u)
    return np.polyval(tuple(reversed(coefficients)), u)


def _usable_coefficients(
    coefficients: tuple[float, ...], band: str, axis: str
) -> tuple[float, ...]:
    """Reject the unpopulated-calibration sentinel, keep anything genuine."""
    if tuple(coefficients) == UNPOPULATED_LOS_SENTINEL:
        logger.warning(
            "Band %s: %s line-of-sight calibration is unpopulated (a constant "
            "1 rad offset is not a physical correction); using zero. The "
            "correction can be estimated from band co-registration.",
            band,
            axis,
        )
        return ()
    return tuple(coefficients)
