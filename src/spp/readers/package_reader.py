"""Reader for filesystem acquisition packages.

The :class:`PackageReader` discovers the assets of an acquisition package and
assembles them into an :class:`~spp.core.acquisition.Acquisition`. It avoids
hard-coded file names by preferring the STAC item's asset table for band
discovery and glob patterns for the calibration assets, falling back to
conventional names when a STAC item is not present.

Only raster *metadata* (size, dtype, nodata) is read here — the pixel data is
left on disk and read in windows by downstream stages.
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import rasterio
from rasterio.errors import NotGeoreferencedWarning

from spp.core.acquisition import Acquisition, ImagerConfiguration
from spp.core.band import Band
from spp.core.calibration_parameters import (
    CalibrationParameters,
    RadiometricCalibration,
)
from spp.readers.base import Reader


class PackageReader(Reader):
    """Read an acquisition package directory into an :class:`Acquisition`.

    Parameters
    ----------
    package_dir:
        Path to the package root (the directory containing the band rasters,
        ``metadata.json``, ``ancillary.json`` and the ``calibration/`` subtree).
    calibration_subdir:
        Name of the calibration subdirectory inside the package.

    Notes
    -----
    The reader is tolerant about exact file names: band rasters are discovered
    from the STAC item assets (or, failing that, from ``*.tif*`` files), and the
    calibration documents are located by glob so that payload-specific suffixes
    do not need to be hard-coded.
    """

    def __init__(
        self,
        package_dir: Path | str,
        *,
        calibration_subdir: str = "calibration",
    ) -> None:
        self.package_dir = Path(package_dir)
        self.calibration_dir = self.package_dir / calibration_subdir
        if not self.package_dir.is_dir():
            raise FileNotFoundError(f"Package directory not found: {self.package_dir}")

    # -- public API ---------------------------------------------------------

    def read(self) -> Acquisition:
        """Discover and parse the package, returning an :class:`Acquisition`."""
        stac = self._read_stac_item()
        session = self._read_session_metadata()
        imager_config = self._build_imager_config(session)
        calibration = self._read_calibration()
        ancillary = self._read_json_if_present(self.package_dir / "ancillary.json")

        name_to_id = self._band_name_to_id()
        bands = self._build_bands(stac, name_to_id, imager_config)
        temps, sample_times = self._temperature_samples(session)

        scene_id, acquired_at = self._scene_identity(stac)
        sample_lines = self._temperature_sample_lines(
            session,
            imager_config,
            acquired_at,
            sample_times,
            self._time_sync_offset(stac),
        )

        return Acquisition(
            scene_id=scene_id,
            bands=bands,
            calibration=calibration,
            imager_config=imager_config,
            sensor_temperatures=temps,
            temperature_sample_lines=sample_lines,
            acquired_at=acquired_at,
            metadata=self._provenance_metadata(session),
            ancillary=ancillary or {},
        )

    # -- discovery helpers --------------------------------------------------

    def _read_stac_item(self) -> dict | None:
        """Return the STAC item document at the package root, if any."""
        for path in sorted(self.package_dir.glob("*.json")):
            doc = self._read_json_if_present(path)
            if isinstance(doc, dict) and "stac_version" in doc:
                return doc
        return None

    def _band_raster_paths(self, stac: dict | None) -> dict[str, Path]:
        """Map band name to raster path, preferring STAC assets."""
        if stac is not None:
            assets = stac.get("assets", {})
            paths = {
                name: self.package_dir / asset["href"]
                for name, asset in assets.items()
                if "tiff" in str(asset.get("type", "")).lower()
            }
            if paths:
                return paths
        # Fallback: every GeoTIFF at the package root, keyed by file stem.
        return {
            p.stem: p
            for p in sorted(self.package_dir.glob("*.tif*"))
        }

    def _calibration_file(self, pattern: str) -> Path:
        """Locate a single calibration document by glob ``pattern``."""
        matches = sorted(self.calibration_dir.glob(pattern))
        if not matches:
            raise FileNotFoundError(
                f"No calibration file matching {pattern!r} in {self.calibration_dir}"
            )
        return matches[0]

    # -- parsing helpers ----------------------------------------------------

    def _read_session_metadata(self) -> dict:
        """Return the (single) session block from ``metadata.json``."""
        doc = self._read_json_if_present(self.package_dir / "metadata.json")
        if not doc or "Sessions" not in doc:
            raise ValueError("metadata.json missing or has no 'Sessions' block")
        sessions = doc["Sessions"]
        session_id = next(iter(sessions))
        return sessions[session_id]

    def _build_imager_config(self, session: dict) -> ImagerConfiguration:
        cfg = session["ImagerConfiguration"]
        return ImagerConfiguration(
            line_period_us=cfg["LinePeriod"],
            spectral_bands=cfg["SpectralBands"],
            band_setup=list(cfg["BandSetup"]),
            band_start_row=list(cfg["BandStartRow"]),
            band_cwl=list(cfg["BandCWL"]),
            scan_direction=cfg["ScanDirection"],
            binning_factor=cfg["BinningFactor"],
        )

    def _read_calibration(self) -> CalibrationParameters:
        cpf = self._read_json(self._calibration_file("CPF*.json"))
        radiometric = [
            RadiometricCalibration(
                band=entry["band"],
                start_row=entry["row"],
                tdi=entry["tdi"],
                absolute=float(entry["absolute"]),
                # The absolute offset is optional: the reference CPF omits the
                # key entirely (offset = 0), but other CPFs may carry it as an
                # explicit JSON ``null``. ``dict.get`` only substitutes the
                # default for a *missing* key, so a present-but-null value would
                # otherwise reach ``float(None)`` and raise. ``_optional_float``
                # treats both missing and null as "no offset".
                absolute_offset=self._optional_float(entry.get("absolute_offset")),
                non_uniformity=np.asarray(entry["coefficients"][0], dtype=np.float64),
                thermal_intercept=np.asarray(entry["coefficients"][1], dtype=np.float64),
                thermal_gradient=np.asarray(entry["coefficients"][2], dtype=np.float64),
                uuid=entry.get("UUID"),
            )
            for entry in cpf["radiometric"]
        ]
        return CalibrationParameters(
            radiometric=radiometric,
            geometric=cpf.get("geometric"),
        )

    def _band_name_to_id(self) -> dict[str, int]:
        """Authoritative band name -> detector id map, from the filter file."""
        filters = self._read_json(self._calibration_file("filters*.json"))
        return {entry["name"]: int(entry["id"]) for entry in filters}

    def _build_bands(
        self,
        stac: dict | None,
        name_to_id: dict[str, int],
        imager_config: ImagerConfiguration,
    ) -> dict[str, Band]:
        raster_paths = self._band_raster_paths(stac)
        bands: dict[str, Band] = {}
        for name, path in raster_paths.items():
            if not path.exists():
                continue  # e.g. the 'raw' asset references a file not shipped
            band_id = name_to_id.get(name)
            cwl = (
                float(imager_config.band_cwl[band_id])
                if band_id is not None
                and band_id < len(imager_config.band_cwl)
                else float("nan")
            )
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", NotGeoreferencedWarning)
                with rasterio.open(path) as ds:
                    width, height = ds.width, ds.height
                    dtype = ds.dtypes[0]
                    nodata = ds.nodata
            bands[name] = Band(
                name=name,
                band_id=band_id if band_id is not None else -1,
                cwl_nm=cwl,
                path=path,
                width=width,
                height=height,
                dtype=dtype,
                nodata=nodata,
            )
        # Order bands by detector id for stable, predictable iteration.
        return dict(sorted(bands.items(), key=lambda kv: kv[1].band_id))

    @staticmethod
    def _provenance_metadata(session: dict) -> dict:
        """Lightweight copy of the session metadata for provenance.

        The per-scene timing block (``Scenes``) holds per-line records and can be
        tens of megabytes; it is not needed for radiometric processing and would
        otherwise be carried in memory on every acquisition. It is replaced by a
        small summary; the full record remains on disk for levels that need it
        (e.g. geolocation).
        """
        trimmed = {k: v for k, v in session.items() if k != "Scenes"}
        scenes = session.get("Scenes")
        if isinstance(scenes, dict):
            trimmed["Scenes"] = {"_omitted": True, "count": len(scenes)}
        return trimmed

    @staticmethod
    def _temperature_samples(session: dict) -> tuple[np.ndarray, np.ndarray]:
        """Detector temperature samples and their imager timestamps.

        Returns two parallel arrays ``(temperatures, imager_times)``. Only
        samples carrying both a temperature and a timestamp are kept, so the two
        arrays stay aligned; ``imager_times`` is empty when no timestamps are
        present.
        """
        telemetry = session.get("ImagerTelemetry", [])
        temps: list[float] = []
        times: list[float] = []
        for t in telemetry:
            if not isinstance(t, dict) or "SensorTemperature" not in t:
                continue
            temps.append(t["SensorTemperature"])
            times.append(t["ImagerTime"] if "ImagerTime" in t else np.nan)
        temps_arr = np.asarray(temps, dtype=np.float64)
        times_arr = np.asarray(times, dtype=np.float64)
        if times_arr.size and np.isnan(times_arr).any():
            times_arr = np.asarray([], dtype=np.float64)  # incomplete timing
        return temps_arr, times_arr

    def _temperature_sample_lines(
        self,
        session: dict,
        imager_config: ImagerConfiguration,
        acquired_at: str | None,
        sample_times: np.ndarray,
        time_sync_offset_s: float,
    ) -> np.ndarray | None:
        """Map each temperature sample to its along-track line index.

        Uses the detector timestamps (``ImagerTime``) and the per-line timing to
        place each sample at its true line, which is exact even when the
        telemetry is irregularly spaced or extends beyond the imaging window.
        Returns ``None`` when any required timing input is missing, so the
        calibrator falls back to uniform spreading.

        The line clock is anchored to the imager clock via the ``TimeSync``
        block (which ties an ``ImagerTime`` to a platform epoch time) and the
        acquisition start time; line ``l`` is then at
        ``t_line0 + l * line_period``.
        """
        if sample_times.size == 0 or acquired_at is None:
            return None
        anchor = self._time_sync_anchor(session)
        line_period = imager_config.line_period_us
        start_ms = self._iso_to_epoch_ms(acquired_at)
        if anchor is None or not line_period or start_ms is None:
            return None
        imager_ref_us, platform_ref_ms = anchor
        # ImagerTime of the first image line: convert the acquisition start
        # (epoch ms) into the imager clock (µs) through the sync anchor.
        t_line0_us = (
            imager_ref_us
            + (start_ms - platform_ref_ms) * 1000.0
            + time_sync_offset_s * 1e6
        )
        return (sample_times - t_line0_us) / float(line_period)

    @staticmethod
    def _time_sync_anchor(session: dict) -> tuple[float, float] | None:
        """Return ``(imager_time_us, platform_time_ms)`` from the TimeSync block.

        The anchor is the entry tying the imager clock to the platform epoch
        clock (``TimeFormat == 1``); the bare PPS entries are ignored.
        """
        for entry in session.get("TimeSync", []):
            if (
                isinstance(entry, dict)
                and entry.get("TimeFormat") == 1
                and "ImagerTime" in entry
                and "PlatformTime" in entry
            ):
                return float(entry["ImagerTime"]), float(entry["PlatformTime"])
        return None

    @staticmethod
    def _iso_to_epoch_ms(iso: str) -> float | None:
        """Convert an ISO-8601 UTC timestamp to epoch milliseconds."""
        from datetime import datetime, timezone

        try:
            dt = datetime.fromisoformat(iso.replace("Z", "+00:00"))
        except (ValueError, AttributeError):
            return None
        if dt.tzinfo is None:
            dt = dt.replace(tzinfo=timezone.utc)
        return dt.timestamp() * 1000.0

    @staticmethod
    def _time_sync_offset(stac: dict | None) -> float:
        """Calibration clock offset (seconds) from the STAC properties, if any."""
        if not stac:
            return 0.0
        for key, value in stac.get("properties", {}).items():
            if key.endswith("time_sync_offset"):
                try:
                    return float(value)
                except (TypeError, ValueError):
                    return 0.0
        return 0.0

    @staticmethod
    def _optional_float(value: object, default: float = 0.0) -> float:
        """Coerce an optional CPF numeric field to ``float``.

        A field such as the absolute offset may be **absent** from a CPF entry or
        present as an explicit JSON ``null``; both mean "not specified" and map to
        ``default``. Any other value is converted with ``float`` so a malformed
        (non-numeric) value still fails loudly rather than silently.
        """
        if value is None:
            return default
        return float(value)

    @staticmethod
    def _scene_identity(stac: dict | None) -> tuple[str, str | None]:
        if stac is None:
            return ("unknown", None)
        props = stac.get("properties", {})
        acquired_at = props.get("start_datetime") or props.get("datetime")
        return (stac.get("id", "unknown"), acquired_at)

    # -- io -----------------------------------------------------------------

    @staticmethod
    def _read_json(path: Path) -> object:
        with path.open("r", encoding="utf-8") as fh:
            return json.load(fh)

    @classmethod
    def _read_json_if_present(cls, path: Path) -> object | None:
        if not path.is_file():
            return None
        try:
            return cls._read_json(path)
        except (json.JSONDecodeError, UnicodeDecodeError):
            return None
