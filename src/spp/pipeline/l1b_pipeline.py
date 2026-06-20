"""End-to-end L0/L1A -> L1B radiance pipeline.

The pipeline wires together the reader, calibrator, quality validator and
writer. It owns only orchestration and the windowing loop; every algorithm lives
in its own module. For each band it streams windows of DN through the calibrator
into the writer while a QA accumulator measures the radiance in the same pass.
"""

from __future__ import annotations

import logging
import time
import warnings
from pathlib import Path
from typing import Callable

import rasterio
from rasterio.errors import NotGeoreferencedWarning
from rasterio.windows import Window

from spp.calibration.l1b_calibrator import L1BCalibrator
from spp.core.acquisition import Acquisition
from spp.core.product import Product
from spp.pipeline.base import ProcessingPipeline
from spp.products.base import ProductWriter
from spp.qa.radiometric_validator import RadiometricValidator
from spp.readers.base import Reader

logger = logging.getLogger(__name__)

#: Factory mapping an acquisition to the calibrator that processes it.
CalibratorFactory = Callable[[Acquisition], L1BCalibrator]


class L1BPipeline(ProcessingPipeline):
    """Produce an L1B at-aperture radiance product from an acquisition package.

    Parameters
    ----------
    reader:
        Reader that yields the :class:`Acquisition` to process.
    writer:
        Raster writer for the calibrated bands.
    validator:
        Radiometric QA validator (provides per-band accumulators).
    calibrator_factory:
        Callable building a calibrator from the acquisition. Defaults to
        :class:`~spp.calibration.l1b_calibrator.L1BCalibrator`.
    bands:
        Subset of band names to process. ``None`` processes all bands.
    window_lines:
        Number of along-track lines processed per window.
    """

    LEVEL = "L1B"

    def __init__(
        self,
        reader: Reader,
        writer: ProductWriter,
        validator: RadiometricValidator,
        *,
        calibrator_factory: CalibratorFactory = L1BCalibrator,
        bands: list[str] | None = None,
        window_lines: int = 2048,
    ) -> None:
        self.reader = reader
        self.writer = writer
        self.validator = validator
        self.calibrator_factory = calibrator_factory
        self.bands = bands
        self.window_lines = window_lines

    def run(self) -> Product:
        """Execute the pipeline and return the generated :class:`Product`."""
        t0 = time.perf_counter()
        acquisition = self.reader.read()
        calibrator = self.calibrator_factory(acquisition)
        band_names = self.bands or acquisition.band_names
        logger.info(
            "Read acquisition %s (%d bands to process)",
            acquisition.scene_id,
            len(band_names),
        )

        product = Product(
            level=self.LEVEL,
            scene_id=acquisition.scene_id,
            units=calibrator.units,
        )
        band_reports: dict[str, dict] = {}

        for index, name in enumerate(band_names, start=1):
            report, path = self._process_band(
                acquisition, calibrator, name, index, len(band_names)
            )
            band_reports[name] = report
            product.bands[name] = path

        product.metadata.update(self._build_metadata(acquisition, calibrator))
        product.metadata["qa"] = self.validator.summarise(
            band_reports, scene_id=acquisition.scene_id, level=self.LEVEL
        )
        product.provenance.update(self._build_provenance(acquisition, calibrator))

        elapsed = time.perf_counter() - t0
        product.metadata["runtime_seconds"] = round(elapsed, 2)
        logger.info(
            "L1B product complete: %d bands in %.1f s (QA passed=%s)",
            len(band_names),
            elapsed,
            product.metadata["qa"]["passed"],
        )
        return product

    # -- internals ----------------------------------------------------------

    def _process_band(
        self,
        acquisition: Acquisition,
        calibrator: L1BCalibrator,
        name: str,
        index: int,
        total: int,
    ) -> tuple[dict, "Path"]:
        band = acquisition.band(name)
        accumulator = self.validator.accumulator(name)
        t0 = time.perf_counter()
        logger.info("Calibrating band %s (%d/%d) %dx%d", name, index, total, band.width, band.height)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", NotGeoreferencedWarning)
            with rasterio.open(band.path) as src, self.writer.open_band(
                name,
                width=band.width,
                height=band.height,
                dtype="float32",
                nodata=calibrator.nodata,
            ) as band_writer:
                for line_start in range(0, band.height, self.window_lines):
                    height = min(self.window_lines, band.height - line_start)
                    window = Window(0, line_start, band.width, height)
                    dn = src.read(1, window=window)
                    radiance = calibrator.calibrate(name, dn, line_start=line_start)
                    band_writer.write_block(radiance, line_start)
                    accumulator.update(radiance, dn=dn)
                band_path = band_writer.path

        report = accumulator.result()
        logger.info(
            "Band %s done in %.1f s: mean=%.3f flags=%s",
            name,
            time.perf_counter() - t0,
            report["mean"],
            report["flags"],
        )
        return report, band_path

    def _build_metadata(self, acquisition: Acquisition, calibrator: L1BCalibrator) -> dict:
        cfg = acquisition.imager_config
        return {
            "scene_id": acquisition.scene_id,
            "acquired_at": acquisition.acquired_at,
            "level": self.LEVEL,
            "units": calibrator.units,
            "imager_configuration": {
                "line_period_us": cfg.line_period_us,
                "band_setup": cfg.band_setup,
                "band_start_row": cfg.band_start_row,
                "scan_direction": cfg.scan_direction,
            },
        }

    def _build_provenance(self, acquisition: Acquisition, calibrator: L1BCalibrator) -> dict:
        cal = acquisition.calibration
        return {
            "source_scene_id": acquisition.scene_id,
            "calibration_uuids": {
                r.band: r.uuid for r in cal.radiometric
            },
            "processing": {
                "temperature_model": "per_line_linear_interpolation",
                "window_lines": self.window_lines,
                "nodata": calibrator.nodata,
                "clipping": "none",
            },
        }
