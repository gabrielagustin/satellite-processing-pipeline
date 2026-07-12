"""The L1C pipeline: L1B radiance in sensor coordinates to a projected map grid.

Orchestration only. Every decision this stage makes lives in the modules it calls —
the sensor model, the terrain, the self-calibration — and every one of them reports
what it could *not* establish as loudly as what it could.

The stage runs whether or not each piece succeeds, and records which did:

* Without a reachable elevation model it georeferences on the ellipsoid, and says so.
* Without enough textured imagery it skips the self-calibration, and says so.
* It never resolves the two telemetry parities that geometry cannot see, and says so.

A product that is missing a correction and admits it is worth more than one that is
missing the same correction and looks finished.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import rasterio
from rasterio.enums import Resampling

from spp.geometry import conventions, geoloc_grid
from spp.geometry.camera import Camera
from spp.geometry.ephemeris import Attitude, Ephemeris
from spp.geometry.frames import ecef_to_geodetic, eci_to_ecef_matrix
from spp.geometry.sensor_model import Convention, SensorModel
from spp.geometry.terrain import DEMTerrain, EllipsoidTerrain, TerrainModel
from spp.geometry.timing import LineTiming
from spp.refine import absolute, matcher, reference, relative, scene
from spp.resample.grid import TargetGrid, native_gsd
from spp.resample.warper import GeolocWarper

logger = logging.getLogger(__name__)


class _Stages:
    """Announces each stage as it starts, and says plainly when one is skipped."""

    def __init__(self, plan: list[str]) -> None:
        self._plan = plan
        self._number = 0

    def begin(self, title: str) -> None:
        if title not in self._plan:
            return
        self._number = self._plan.index(title) + 1
        logger.info("[%d/%d] %s", self._number, len(self._plan), title)

    def skip(self, reason: str) -> None:
        logger.info("      skipped — %s", reason)


RESOLVED_CONVENTION = Convention(
    ephemeris_frame="eci",
    attitude_frame="lvlh",
    quaternion_order="scalar_first",
    quaternion_direction="body_to_ref",
    scan_direction=1,
    column_axis="y",
    column_sign=-1,
    row_sign=-1,
)
"""The telemetry conventions.

Six were resolved by geometry (``spp.geometry.conventions.resolve``). The remaining two —
the scan direction and the detector column sign — are **mirrors**, invisible to every
geometric probe, and they are resolved against the **terrain**: water is near-black in the
near infrared, so a correct parity makes bright coincide with high ground and a mirrored
one anti-correlates (``resolve_parities``).

``column_sign`` was **wrong here until a user overlaid the product on a basemap and saw
that it did not match.** Every geometric check had passed. That is precisely what the
harness meant when it reported the parities as unresolved rather than guessing them, and
it is why the pipeline now resolves them from image content instead of carrying them as
assumptions.
"""


@dataclass
class L1CResult:
    """What the run produced, and what it could not."""

    products: dict[str, Path] = field(default_factory=dict)
    qa: dict = field(default_factory=dict)
    grid: TargetGrid | None = None
    stack: Path | None = None


class L1CPipeline:
    """Georeference, orthorectify and co-register an acquisition.

    Parameters
    ----------
    camera, ephemeris, attitude, timing:
        The sensor model's components (see :mod:`spp.geometry.telemetry`).
    terrain:
        Elevation model. ``None`` georeferences on the ellipsoid — correct over water,
        and displacing every metre of relief elsewhere.
    convention:
        How to read the telemetry. Defaults to the resolved set above.
    reference_band:
        The band the others are co-registered onto.
    refine:
        Estimate the missing per-band line-of-sight calibration from the imagery.
    scene_correct:
        Also correct the part of the band misregistration that a constant offset cannot
        explain — the along-track-varying part, which is this acquisition's attitude
        behaviour rather than the instrument's optics. Applied to the **product**, never
        written to the calibration. Without it the bands will not fully co-register; see
        :mod:`spp.refine.scene`.
    stack:
        Also write a single multi-band raster with every band on the common grid.
    gsd_m:
        Output resolution. Defaults to the native sampling, rounded up.
    step:
        Geolocation-lattice spacing, in detector samples.
    """

    def __init__(
        self,
        camera: Camera,
        ephemeris: Ephemeris,
        attitude: Attitude,
        timing: dict[str, LineTiming],
        *,
        terrain: TerrainModel | None = None,
        convention: Convention = RESOLVED_CONVENTION,
        reference_band: str = "PAN",
        refine: bool = True,
        resolve_parities: bool = True,
        correct_absolute: bool = True,
        band_wavelengths: dict[str, float] | None = None,
        reference_search: tuple[tuple[float, float, float, float], object, Path] | None = None,
        scene_correct: bool = True,
        stack: bool = True,
        gsd_m: float | None = None,
        step: int = geoloc_grid.DEFAULT_STEP,
    ) -> None:
        self.camera = camera
        self.ephemeris = ephemeris
        self.attitude = attitude
        self.timing = timing
        self.terrain = terrain or EllipsoidTerrain()
        self.convention = convention
        self.reference_band = reference_band
        self.refine = refine
        self.resolve_parities = resolve_parities
        self.correct_absolute = correct_absolute
        self.band_wavelengths = band_wavelengths or {}
        self.reference_search = reference_search
        self.reference_image: reference.ReferenceImage | None = None
        self.scene_correct = scene_correct
        self.stack = stack
        self.gsd_m = gsd_m
        self.step = step

        self.model = SensorModel(camera, ephemeris, attitude, timing, convention=convention)
        self.native_gsd_m = self._native_gsd()

    # -- public API ---------------------------------------------------------

    def run(self, l1b_dir: Path, output_dir: Path, bands: list[str] | None = None) -> L1CResult:
        """Process every band, and write the products and the quality report."""
        bands = bands or sorted(self.timing)
        paths = {band: l1b_dir / f"{band}.tif" for band in bands}
        missing = [b for b, p in paths.items() if not p.exists()]
        if missing:
            raise FileNotFoundError(
                f"No L1B raster for band(s) {', '.join(missing)} in {l1b_dir}. "
                "Run the L1B stage first (`spp-l1b`)."
            )

        stages = self._stage_plan(bands)
        logger.info("L1C: %d stage(s) over %d band(s) — %s", len(stages), len(bands), ", ".join(bands))
        for number, title in enumerate(stages, start=1):
            logger.info("   %d. %s", number, title)

        step = _Stages(stages)

        step.begin("Resolve the telemetry's mirrors against the terrain")
        convention = self.convention
        parity_report = None
        if self.resolve_parities and isinstance(self.terrain, DEMTerrain):
            convention, parity_report = self._resolve_parities(paths, bands)
            if convention != self.convention:
                self.convention = convention
                self.model = SensorModel(
                    self.camera, self.ephemeris, self.attitude, self.timing,
                    convention=convention,
                )

        qa: dict = {
            "conventions": {
                "resolved": convention.describe(),
                "parities": parity_report or {"note": "not resolved; assumed"},
            },
            "grid": {},
            "terrain": self._terrain_qa(),
            "coregistration": {},
            "los_correction": {},
            "flags": ["no_gnss_lock", "absolute_accuracy_unvalidated"],
        }
        if not (parity_report and parity_report.get("decisive")):
            qa["flags"].append("conventions_assumed")

        step.begin("Correct the absolute pointing against an independent reference")
        if self.correct_absolute and isinstance(self.terrain, DEMTerrain):
            self._correct_absolute(paths, bands, qa)
        else:
            step.skip("no elevation model, or disabled")

        step.begin(f"Self-calibrate the interior orientation against {self.reference_band}")
        model = self.model
        matches: dict[str, list] = {}
        if self.refine and self.reference_band in bands:
            model, matches = self._self_calibrate(paths, bands, qa)
        elif self.refine:
            logger.warning(
                "Reference band %s not among the bands processed; skipping the "
                "self-calibration. The bands will not be co-registered.",
                self.reference_band,
            )
            qa["flags"].append("coregistration_not_achieved")

        step.begin(f"Build the geolocation grids (lattice step {self.step} px)")
        grids = {}
        with rasterio.open(paths[bands[0]]) as probe:
            n_lines, n_columns = probe.height, probe.width
        for band in bands:
            grids[band] = geoloc_grid.build(
                model,
                band,
                n_lines=n_lines,
                n_columns=n_columns,
                terrain=self.terrain,
                step=self.step,
            )

        step.begin("Correct this acquisition's attitude (product only, not the calibration)")
        if self.scene_correct and matches:
            grids = self._scene_correct(model, grids, matches, qa)
        else:
            step.skip("no matches, or disabled")

        step.begin("Choose the common map grid")
        target = TargetGrid.covering(
            {b: g.bounds for b, g in grids.items()},
            native_gsd_m=self.native_gsd_m,
            gsd_m=self.gsd_m,
        )
        qa["grid"] = {
            "crs": target.crs.to_string(),
            "gsd_m": target.gsd_m,
            "native_gsd_m": round(target.native_gsd_m, 3),
            "shape": [target.height, target.width],
            "zone_straddle": target.zone_straddle,
        }
        if target.zone_straddle:
            qa["flags"].append("zone_straddle")

        qa["interpolation"] = {
            "lattice_step_px": self.step,
            "max_error_px": round(
                geoloc_grid.interpolation_error_px(
                    model,
                    grids[bands[0]],
                    bands[0],
                    terrain=self.terrain,
                    gsd_m=self.native_gsd_m,
                ),
                4,
            ),
        }

        step.begin(f"Resample {len(bands)} band(s) onto that grid")
        warper = GeolocWarper()
        products = {}
        for band in bands:
            logger.info("      warping %s", band)
            products[band] = warper.warp(
                paths[band], grids[band], target, output_dir / f"{band}.tif"
            )

        step.begin("Write the co-registered stack")
        stack_path = None
        if self.stack and len(products) > 1:
            stack_path = self._write_stack(products, target, output_dir)
        else:
            step.skip("single band, or disabled")

        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / "qa_report_l1c.json").write_text(json.dumps(qa, indent=2))
        logger.info("Wrote %d band(s) and qa_report_l1c.json to %s", len(products), output_dir)

        return L1CResult(products=products, qa=qa, grid=target, stack=stack_path)

    def _stage_plan(self, bands: list[str]) -> list[str]:
        """What this run is going to do, announced before it does it.

        A log that only says what has happened leaves the reader guessing how much is
        left and which steps were skipped. Stating the plan first makes a skipped stage
        visible as a *gap*, not as an absence nobody notices.
        """
        plan = []
        if self.resolve_parities and isinstance(self.terrain, DEMTerrain):
            plan.append("Resolve the telemetry's mirrors against the terrain")
        if self.correct_absolute and isinstance(self.terrain, DEMTerrain):
            plan.append("Correct the absolute pointing against an independent reference")
        if self.refine and self.reference_band in bands:
            plan.append(f"Self-calibrate the interior orientation against {self.reference_band}")
        plan.append(f"Build the geolocation grids (lattice step {self.step} px)")
        if self.scene_correct:
            plan.append("Correct this acquisition's attitude (product only, not the calibration)")
        plan.append("Choose the common map grid")
        plan.append(f"Resample {len(bands)} band(s) onto that grid")
        if self.stack and len(bands) > 1:
            plan.append("Write the co-registered stack")
        return plan

    def _resolve_parities(self, paths: dict[str, Path], bands: list[str]):
        """Settle the two mirrors that geometry cannot see, against the elevation model.

        Near-infrared is used when available: water is near-black in it, so the
        land/water contrast that carries the signal is strongest.
        """
        band = self._match_band(bands)
        logger.info("      using band %s (water is darkest in the near infrared)", band)
        with rasterio.open(paths[band]) as src:
            image = src.read(1)
        return conventions.resolve_parities(
            self.camera, self.ephemeris, self.attitude, self.timing,
            band=band, image=image, terrain=self.terrain, base=self.convention,
        )

    def _match_band(self, bands: list[str]) -> str:
        """The band used to match against an external reference.

        Near-infrared first: water is near-black in it, so the land/water boundary is
        sharp and vegetation gives texture.
        """
        return next((b for b in ("NIR", "RE3", "RE2", "R") if b in bands), bands[0])

    def _find_reference(self, band: str) -> reference.ReferenceImage | None:
        """Fetch a reference orthoimage **at this band's own wavelength**.

        The wavelength is not a detail. Matching our 665 nm red against a reference's
        842 nm near-infrared correlates two different pictures of the same ground: on the
        reference acquisition it moved the measured bias from 919 m to 283 m and inflated
        the window scatter from 1.6 lattice cells to 9.0. A confident, precise, wrong
        answer -- which is why the reference is searched *after* the band is chosen, not
        before.
        """
        if self.reference_search is None:
            return None
        bounds, acquired, destination = self.reference_search
        wavelength = self.band_wavelengths.get(band)
        if wavelength is None:
            logger.warning(
                "No central wavelength known for band %s; cannot match it to a reference "
                "band. Falling back to the terrain's coastline.",
                band,
            )
            return None
        return reference.find(
            bounds, acquired, destination, wavelength_nm=int(round(wavelength))
        )

    def _correct_absolute(self, paths: dict[str, Path], bands: list[str], qa: dict) -> None:
        """Correct the pointing bias against the terrain -- the only external reference here.

        Without a satellite-navigation lock this bias is expected, and on the reference
        acquisition it is about a kilometre. Nothing else in this level can see it: the
        bands are co-registered onto *each other*, and the delivered footprint was derived
        from the same telemetry the model consumes. The terrain's coastline is the one
        reference in the run that the telemetry did not produce.
        """
        band = self._match_band(bands)
        self.reference_image = self._find_reference(band)
        logger.info("      matching band %s against the reference", band)
        with rasterio.open(paths[band]) as src:
            image = src.read(1)

        ground_speed = self._ground_speed()

        # Prefer a reference orthoimage: it matches *texture*, and pins the product to a
        # pixel. The terrain's coastline is a real reference too, and a coarse one -- it
        # can only be as sharp as a shoreline is, which on the reference acquisition left
        # a 58 m discrepancy against the orthoimage. Fall back to it, and say which was
        # used, rather than quietly reporting one number for two different measurements.
        correction = None
        source = "terrain coastline"
        if self.reference_image is not None:
            correction = absolute.estimate_against_reference(
                self.model,
                band,
                image,
                str(self.reference_image.path),
                self.terrain,
                ground_speed_m_s=ground_speed,
            )
            if correction.fitted:
                source = (
                    f"reference orthoimage ({self.reference_image.acquired}, "
                    f"{self.reference_image.days_from_target} day(s) away, "
                    f"{self.reference_image.cloud_cover:.0f}% cloud)"
                )
            else:
                logger.warning(
                    "The reference orthoimage did not settle the geolocation (%s); "
                    "falling back to the terrain's coastline.",
                    correction.reason,
                )
                correction = None

        if correction is None:
            correction = absolute.estimate(
                self.model, band, image, self.terrain, ground_speed_m_s=ground_speed
            )

        qa["absolute"] = {
            "reference": source,
            "offset_before_m": round(correction.offset_before_m, 1),
            "roll_urad": round(correction.roll_rad * 1e6, 1),
            "pitch_urad": round(correction.pitch_rad * 1e6, 1),
            "equivalent_clock_offset_s": round(correction.equivalent_clock_offset_s, 4),
            "peak": round(correction.peak, 4),
            "fitted": correction.fitted,
            "reason": correction.reason,
            "note": (
                "The along-track correction is equally explained by a clock offset; the "
                "two are not separable from one strip. Attributing it to pitch is a "
                "convention, not a measurement."
            ),
        }
        if not correction.fitted:
            return

        self.model = absolute._with_boresight(
            self.model, roll=correction.roll_rad, pitch=correction.pitch_rad
        )
        self.camera = self.model.camera
        # The bias is corrected, but it is corrected against a *coastline in an
        # elevation model*, not against an orthoimage. That is a real reference and a
        # limited one: it is good to the DEM's own resolution and to how sharply the
        # shoreline is defined -- tens of metres, not sub-pixel. Replacing one honest
        # caveat with silence would be worse than the bias.
        if "absolute_accuracy_unvalidated" in qa["flags"]:
            qa["flags"].remove("absolute_accuracy_unvalidated")

        # Only the *weaker* outcome earns a caveat. Corrected against a reference
        # orthoimage, the geolocation has been established, and listing it under "not
        # established" would be the same mistake as burying it -- just in the other
        # direction. The remaining degeneracy is reported with the measurement itself.
        if "orthoimage" not in source:
            qa["flags"].append("absolute_accuracy_terrain_only")

    def _ground_speed(self) -> float:
        """Speed of the sub-satellite point, from the ephemeris."""
        mid = float(np.median(self.ephemeris.times))
        position, velocity = self.ephemeris.interpolate(np.array([mid]))
        radius = float(np.linalg.norm(position[0]))
        return float(np.linalg.norm(velocity[0]) * 6_371_008.8 / radius)

    def _scene_correct(
        self, model: SensorModel, grids: dict, matches: dict, qa: dict
    ) -> dict:
        """Bend each band's geolocation by the part a constant offset cannot explain.

        This is applied to the **product**, not to the calibration. It describes this
        acquisition's pointing, and feeding it back into the instrument model would make
        every future acquisition wrong. See :mod:`spp.refine.scene`.
        """
        logger.info("Correcting the scene-local (attitude) part of the misregistration...")
        qa["scene_correction"] = {
            "note": (
                "Applied to this product only. Describes THIS acquisition's pointing, "
                "not the instrument's optics; never written to the calibration."
            )
        }
        corrected = dict(grids)
        for band, band_matches in matches.items():
            if band not in grids:
                continue
            correction = scene.estimate(
                model,
                band,
                self.reference_band,
                band_matches,
                terrain=self.terrain,
                gsd_m=self.native_gsd_m,
            )
            qa["scene_correction"][band] = {
                "degree": correction.degree,
                "n_windows": correction.n_windows,
                "rms_before_px": round(correction.residual_before_px, 3),
                "rms_after_px": round(correction.residual_after_px, 3),
                "fitted": correction.fitted,
                "reason": correction.reason,
            }
            if correction.fitted:
                lat0 = float(np.nanmean(grids[band].lat))
                corrected[band] = correction.apply(grids[band], lat0=lat0)
                if "coregistration_not_achieved" in qa["flags"]:
                    qa["flags"].remove("coregistration_not_achieved")
        return corrected

    def _write_stack(self, products: dict, target: TargetGrid, output_dir: Path) -> Path:
        """One multi-band raster: every band on the same grid, co-registered.

        Written **by window, all bands together** — not band by band over the whole
        image, which is the obvious way and which silently destroys the data.

        In a *compressed, tiled* GeoTIFF a tile cannot be revisited: once it has been
        written it is compressed and closed. With the default pixel interleave, a tile
        holds every band, so writing band 1 across the whole image flushes every tile —
        and bands 2..n are then dropped on the floor. No error, no warning, and a file
        that opens fine with one band of data and the rest silently NaN. This pipeline
        produced exactly that before the bug was caught.

        Windowing also keeps the write memory-bounded, which the per-band warp is not.
        """
        path = output_dir / "stack.tif"
        names = list(products)
        block = 1024

        profile = {
            "driver": "GTiff",
            "height": target.height,
            "width": target.width,
            "count": len(names),
            "dtype": "float32",
            "crs": target.crs,
            "transform": target.transform,
            "nodata": float("nan"),
            "tiled": True,
            "blockxsize": 512,
            "blockysize": 512,
            "compress": "deflate",
            "predictor": 3,
            "BIGTIFF": "IF_SAFER",
        }

        sources = [rasterio.open(products[band]) for band in names]
        try:
            with rasterio.open(path, "w", **profile) as dst:
                dst.descriptions = tuple(names)
                for row in range(0, target.height, block):
                    height = min(block, target.height - row)
                    window = ((row, row + height), (0, target.width))
                    tile = np.stack([src.read(1, window=window) for src in sources])
                    dst.write(tile, indexes=list(range(1, len(names) + 1)), window=window)
                dst.build_overviews([2, 4, 8, 16, 32], Resampling.average)
        finally:
            for src in sources:
                src.close()

        logger.info("Wrote stack.tif — %s, co-registered on one grid", ", ".join(names))
        return path

    # -- internals ----------------------------------------------------------

    def _native_gsd(self) -> float:
        """Ground sample distance from the ephemeris — never from a data sheet."""
        mid = float(np.median(self.ephemeris.times))
        position, _ = self.ephemeris.interpolate(np.array([mid]))
        if self.convention.ephemeris_frame == "eci":
            position = np.einsum("mij,mj->mi", eci_to_ecef_matrix(np.array([mid])), position)
        altitude = float(ecef_to_geodetic(position)[0, 2])
        gsd = native_gsd(
            altitude, self.camera.intrinsics.pixel_size_mm, self.camera.intrinsics.focal_length_mm
        )
        logger.info("Altitude %.1f km -> native ground sample distance %.2f m", altitude / 1000, gsd)
        return gsd

    def _terrain_qa(self) -> dict:
        if isinstance(self.terrain, DEMTerrain):
            return {
                "source": "digital elevation model",
                "undulation_range_m": [
                    round(float(self.terrain.geoid.undulation.min()), 2),
                    round(float(self.terrain.geoid.undulation.max()), 2),
                ],
                "void_fraction": round(self.terrain.void_fraction, 3),
            }
        return {
            "source": "ellipsoid (no elevation model)",
            "note": "georeferenced, NOT orthorectified: relief is displaced",
        }

    def _self_calibrate(
        self, paths: dict[str, Path], bands: list[str], qa: dict
    ) -> tuple[SensorModel, dict[str, list]]:
        """Estimate the per-band line-of-sight offsets the calibration file omits."""
        logger.info("Self-calibrating the interior orientation against %s...", self.reference_band)

        with rasterio.open(paths[self.reference_band]) as src:
            reference = src.read(1)

        corrections = {}
        all_matches: dict[str, list] = {}
        achieved = True
        for band in bands:
            if band == self.reference_band:
                continue
            with rasterio.open(paths[band]) as src:
                moving = src.read(1)

            matches = matcher.match_windows(reference, moving, window=384, stride=768)
            all_matches[band] = matches
            correction = relative.estimate(
                model=self.model,
                band=band,
                reference_band=self.reference_band,
                matches=matches,
                terrain=self.terrain,
                gsd_m=self.native_gsd_m,
            )
            qa["coregistration"][band] = {
                "n_windows": correction.n_windows,
                "systematic_before_px": round(correction.residual_before_px, 3),
                "systematic_after_px": round(correction.residual_after_px, 3),
                "scatter_px": round(correction.scatter_px, 3),
                "position_correlation": round(correction.position_correlation, 3),
                "fitted": correction.fitted,
                "reason": correction.reason,
            }
            if correction.fitted:
                corrections[band] = ((correction.along_rad,), (correction.across_rad,))
                qa["los_correction"][band] = {
                    "along_rad": correction.along_rad,
                    "across_rad": correction.across_rad,
                }
            else:
                achieved = False

        if not achieved:
            # Refusing to fit is the correct outcome for a band whose residual varies
            # with position -- but the product is then not co-registered, and must say so.
            qa["flags"].append("coregistration_not_achieved")

        if not corrections:
            logger.warning("No band could be calibrated from the imagery.")
            return self.model, all_matches

        return (
            SensorModel(
                self.camera.with_los(corrections),
                self.ephemeris,
                self.attitude,
                self.timing,
                convention=self.convention,
            ),
            all_matches,
        )
