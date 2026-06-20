"""Radiometric quality assessment for calibrated radiance products.

The validator computes per-band statistics and raises quality flags. It is
**non-destructive**: it only measures, never modifies. To fit the windowed
pipeline it accumulates statistics incrementally — blocks are fed in as they are
produced, then finalised — so a band never needs to be held in memory in full.
The same validator can also assess an already-written product by re-reading it.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from spp.core.product import Product
from spp.qa.base import QAValidator


@dataclass(frozen=True)
class QAThresholds:
    """Thresholds controlling which conditions raise a quality flag.

    Attributes
    ----------
    saturation_dn:
        DN value at/above which a pixel is considered saturated. ``None``
        disables the check (e.g. when DN is not available post-hoc).
    max_negative_fraction:
        Fraction of valid pixels with negative radiance above which the
        ``negative_radiance`` flag is raised.
    max_nodata_fraction:
        Fraction of NoData pixels above which the ``high_nodata`` flag is raised.
    """

    saturation_dn: int | None = None
    max_negative_fraction: float = 0.05
    max_nodata_fraction: float = 0.50


@dataclass
class _BandAccumulator:
    """Incremental statistics accumulator for a single band."""

    band: str
    thresholds: QAThresholds
    n_total: int = 0
    n_valid: int = 0
    n_nodata: int = 0
    n_nonfinite: int = 0
    n_negative: int = 0
    n_saturated: int = 0
    _sum: float = 0.0
    _sumsq: float = 0.0
    _min: float = math.inf
    _max: float = -math.inf
    _saturation_checked: bool = field(default=False)

    def update(self, radiance: np.ndarray, dn: np.ndarray | None = None) -> None:
        """Accumulate one block of radiance (and optionally its source DN)."""
        a = np.asarray(radiance, dtype=np.float64)
        self.n_total += a.size

        nan_mask = np.isnan(a)
        inf_mask = np.isinf(a)
        self.n_nodata += int(nan_mask.sum())
        self.n_nonfinite += int(inf_mask.sum())

        valid = a[~(nan_mask | inf_mask)]
        if valid.size:
            self.n_valid += valid.size
            self._sum += float(valid.sum())
            self._sumsq += float(np.square(valid).sum())
            self._min = min(self._min, float(valid.min()))
            self._max = max(self._max, float(valid.max()))
            self.n_negative += int((valid < 0).sum())

        if dn is not None and self.thresholds.saturation_dn is not None:
            self._saturation_checked = True
            self.n_saturated += int((np.asarray(dn) >= self.thresholds.saturation_dn).sum())

    def result(self) -> dict:
        """Finalise into a metrics + flags dictionary."""
        flags: list[str] = []
        nodata_fraction = self.n_nodata / self.n_total if self.n_total else 0.0
        negative_fraction = self.n_negative / self.n_valid if self.n_valid else 0.0

        if self.n_valid == 0:
            flags.append("all_invalid")
            mean = std = vmin = vmax = float("nan")
        else:
            mean = self._sum / self.n_valid
            var = max(self._sumsq / self.n_valid - mean * mean, 0.0)
            std = math.sqrt(var)
            vmin, vmax = self._min, self._max
            if vmax == vmin:
                flags.append("zero_dynamic_range")

        if nodata_fraction > self.thresholds.max_nodata_fraction:
            flags.append("high_nodata")
        if negative_fraction > self.thresholds.max_negative_fraction:
            flags.append("negative_radiance")
        if self.n_nonfinite > 0:
            flags.append("non_finite_values")
        if self._saturation_checked and self.n_saturated > 0:
            flags.append("saturation")

        metrics = {
            "band": self.band,
            "pixels": self.n_total,
            "valid": self.n_valid,
            "nodata": self.n_nodata,
            "nodata_fraction": nodata_fraction,
            "min": vmin,
            "max": vmax,
            "mean": mean,
            "std": std,
            "negative": self.n_negative,
            "negative_fraction": negative_fraction,
            "non_finite": self.n_nonfinite,
            "flags": flags,
            "passed": not flags,
        }
        if self._saturation_checked:
            metrics["saturated"] = self.n_saturated
        return metrics


class RadiometricValidator(QAValidator):
    """Assess the radiometric quality of a calibrated radiance product.

    Parameters
    ----------
    thresholds:
        Flagging thresholds (see :class:`QAThresholds`).
    window_lines:
        Block height used when re-reading a product in :meth:`validate`.
    """

    def __init__(
        self,
        thresholds: QAThresholds | None = None,
        *,
        window_lines: int = 2048,
    ) -> None:
        self.thresholds = thresholds or QAThresholds()
        self.window_lines = window_lines

    def accumulator(self, band_name: str) -> _BandAccumulator:
        """Create a fresh streaming accumulator for ``band_name``."""
        return _BandAccumulator(band=band_name, thresholds=self.thresholds)

    @staticmethod
    def summarise(band_results: dict[str, dict], *, scene_id: str, level: str) -> dict:
        """Assemble per-band results into a product-level QA report."""
        return {
            "scene_id": scene_id,
            "level": level,
            "passed": all(r["passed"] for r in band_results.values()),
            "bands": band_results,
        }

    def validate(self, product: Product) -> dict:
        """Validate an already-written product by re-reading its rasters.

        Note: DN is not available here, so the saturation check is skipped.
        """
        import rasterio  # local import keeps the entity/QA layer import-light

        band_results: dict[str, dict] = {}
        for name, path in product.bands.items():
            acc = self.accumulator(name)
            with rasterio.open(Path(path)) as ds:
                for _, window in ds.block_windows(1):
                    acc.update(ds.read(1, window=window))
            band_results[name] = acc.result()
        return self.summarise(
            band_results, scene_id=product.scene_id, level=product.level
        )
