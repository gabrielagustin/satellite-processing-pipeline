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
        temps = self._sensor_temperatures(session)

        scene_id, acquired_at = self._scene_identity(stac)

        return Acquisition(
            scene_id=scene_id,
            bands=bands,
            calibration=calibration,
            imager_config=imager_config,
            sensor_temperatures=temps,
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
                absolute_offset=float(entry.get("absolute_offset", 0.0)),
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
    def _sensor_temperatures(session: dict) -> np.ndarray:
        telemetry = session.get("ImagerTelemetry", [])
        temps = [
            t["SensorTemperature"]
            for t in telemetry
            if isinstance(t, dict) and "SensorTemperature" in t
        ]
        return np.asarray(temps, dtype=np.float64)

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
