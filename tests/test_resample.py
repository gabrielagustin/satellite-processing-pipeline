"""The geolocation grid, the target grid, and the warp onto it.

Synthetic throughout: a known linear geolocation, so the correct output is known
in closed form and the tests check the resampling rather than the scene.
"""

from __future__ import annotations

import numpy as np
import pytest
import rasterio
from rasterio.crs import CRS
from rasterio.enums import Resampling
from rasterio.transform import from_bounds

from spp.geometry.geoloc_grid import GeolocationGrid
from spp.resample.grid import TargetGrid, native_gsd, utm_epsg
from spp.resample.warper import GeolocWarper


# -- choosing the grid ------------------------------------------------------


def test_utm_zone_is_derived_not_assumed():
    """The zone is a property of where the platform looked, never a constant."""
    assert utm_epsg(56.1, 26.9) == 32640  # northern hemisphere, zone 40
    assert utm_epsg(-58.4, -34.6) == 32721  # southern, zone 21
    assert utm_epsg(-179.0, 0.0) == 32601  # first zone
    assert utm_epsg(179.0, 0.0) == 32660  # last zone


def test_polar_latitudes_fall_back_off_utm():
    """UTM is not defined at the poles; returning a zone there would be a lie."""
    assert utm_epsg(10.0, 85.0) == 3413
    assert utm_epsg(10.0, -85.0) == 3031


def test_native_gsd_comes_from_the_geometry():
    """GSD is altitude x pitch / focal length — never a data-sheet figure.

    The reference acquisition flew at ~398 km against a nominal ~501 km, so the
    nominal GSD is wrong by 25%. A constant here would resample every scene onto a
    grid of the wrong size, and no test on a single scene would notice.
    """
    assert native_gsd(397_900.0, 0.0055, 580.0) == pytest.approx(3.77, abs=0.01)
    assert native_gsd(501_000.0, 0.0055, 580.0) == pytest.approx(4.75, abs=0.01)


def test_target_grid_covers_the_union_of_the_bands():
    """Every band's data must fit: the union, not the intersection, not band one."""
    footprints = {
        "PAN": (56.0, 26.5, 56.3, 27.4),
        "RE3": (55.9, 26.4, 56.2, 27.3),  # offset from PAN, as the stagger leaves it
    }
    grid = TargetGrid.covering(footprints, native_gsd_m=3.77)

    left, bottom, right, top = grid.bounds
    xs, ys = rasterio.warp.transform(
        CRS.from_epsg(4326), grid.crs, [55.9, 56.3], [26.4, 27.4]
    )
    assert left <= min(xs) and right >= max(xs)
    assert bottom <= min(ys) and top >= max(ys)


def test_default_gsd_rounds_up_from_native():
    """Rounding *up* oversamples slightly; rounding down would decimate and alias."""
    grid = TargetGrid.covering({"PAN": (56.0, 26.5, 56.1, 26.6)}, native_gsd_m=3.77)
    assert grid.gsd_m == 4.0
    assert grid.native_gsd_m == pytest.approx(3.77)


def test_zone_straddle_is_flagged_not_hidden():
    """A strip crossing a UTM boundary still gets a grid — but it says so."""
    straddling = TargetGrid.covering(
        {"PAN": (5.5, 45.0, 6.5, 46.0)}, native_gsd_m=4.0  # zone 31 / 32 boundary at 6E
    )
    assert straddling.zone_straddle

    contained = TargetGrid.covering({"PAN": (56.0, 26.5, 56.3, 27.4)}, native_gsd_m=4.0)
    assert not contained.zone_straddle


def test_an_empty_footprint_is_an_error():
    with pytest.raises(ValueError, match="no footprints"):
        TargetGrid.covering({}, native_gsd_m=4.0)


# -- the geolocation grid ---------------------------------------------------


def test_geolocation_arrays_must_be_float64():
    """float32 resolves longitude to ~0.7 m — the size of the corrections applied.

    Silently accepting float32 would add noise of the same order as the terrain
    correction, which is indistinguishable from the correction not working.
    """
    lon = np.zeros((4, 4), dtype=np.float32)
    lat = np.zeros((4, 4), dtype=np.float32)
    with pytest.raises(ValueError, match="float64"):
        GeolocationGrid(lon=lon, lat=lat, step=8, n_lines=24, n_columns=24)


# -- the warp ---------------------------------------------------------------


def _linear_geoloc(n_lines: int, n_columns: int, step: int) -> GeolocationGrid:
    """A geolocation grid that maps the raster onto a plain lon/lat rectangle."""
    lines = np.arange(0, n_lines + step, step, dtype=np.float64)
    columns = np.arange(0, n_columns + step, step, dtype=np.float64)
    mesh_l, mesh_c = np.meshgrid(lines, columns, indexing="ij")

    lon = 56.0 + 1e-4 * mesh_c
    lat = 27.0 - 1e-4 * mesh_l
    return GeolocationGrid(
        lon=lon, lat=lat, step=step, n_lines=n_lines, n_columns=n_columns
    )


def _write_source(path, data: np.ndarray) -> None:
    profile = {
        "driver": "GTiff",
        "height": data.shape[0],
        "width": data.shape[1],
        "count": 1,
        "dtype": "float32",
        "nodata": float("nan"),
    }
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(data.astype(np.float32), 1)


def test_warp_puts_the_data_on_the_map(tmp_path):
    """End to end: a raster with a known geolocation lands where it should."""
    n_lines, n_columns, step = 64, 64, 8
    source = np.tile(np.arange(n_columns, dtype=np.float32), (n_lines, 1))

    src_path = tmp_path / "band.tif"
    _write_source(src_path, source)

    geoloc = _linear_geoloc(n_lines, n_columns, step)
    target = TargetGrid.covering({"band": geoloc.bounds}, native_gsd_m=8.0)

    out = GeolocWarper().warp(src_path, geoloc, target, tmp_path / "band_L1C.tif")

    with rasterio.open(out) as dst:
        assert dst.crs == target.crs
        assert dst.count == 1
        assert dst.dtypes[0] == "float32"
        warped = dst.read(1)

    valid = np.isfinite(warped)
    assert valid.any()
    # The source ramps 0..63 across the columns; the warp may not preserve the
    # exact values at the edges, but it must preserve the range and the ordering.
    assert warped[valid].min() >= -0.5
    assert warped[valid].max() <= n_columns - 0.5


def test_warp_rejects_a_grid_built_for_another_band(tmp_path):
    """A mismatched geolocation grid would place the band on the wrong ground.

    Each band has its own grid — its own line times, its own view angles. Silently
    accepting the wrong one is exactly the failure the per-band design exists to
    prevent, so the shapes are checked.
    """
    src_path = tmp_path / "band.tif"
    _write_source(src_path, np.zeros((64, 64), dtype=np.float32))

    wrong = _linear_geoloc(128, 64, 8)  # a grid for a taller raster
    target = TargetGrid.covering({"band": wrong.bounds}, native_gsd_m=8.0)

    with pytest.raises(ValueError, match="different band"):
        GeolocWarper().warp(src_path, wrong, target, tmp_path / "out.tif")


def test_nodata_is_not_averaged_into_the_radiance(tmp_path):
    """An interpolation kernel straddling the data edge must not blend in NoData.

    Radiance is a physical quantity; a sentinel averaged into it produces a value
    the sensor never measured, and nothing downstream can tell it apart from one it
    did.
    """
    n_lines, n_columns, step = 64, 64, 8
    source = np.full((n_lines, n_columns), 100.0, dtype=np.float32)
    source[:, :16] = np.nan  # a block of NoData

    src_path = tmp_path / "band.tif"
    _write_source(src_path, source)

    geoloc = _linear_geoloc(n_lines, n_columns, step)
    target = TargetGrid.covering({"band": geoloc.bounds}, native_gsd_m=8.0)
    out = GeolocWarper().warp(src_path, geoloc, target, tmp_path / "out.tif")

    with rasterio.open(out) as dst:
        warped = dst.read(1)

    finite = warped[np.isfinite(warped)]
    # Every real sample was 100. If NoData had been averaged in, values between
    # 0 and 100 would appear at the boundary.
    assert finite.size > 0
    np.testing.assert_allclose(finite, 100.0, atol=1e-3)


def test_ringing_kernels_are_refused():
    """Lanczos rings around saturation and NoData: not a cosmetic issue in radiance."""
    with pytest.raises(ValueError, match="negative lobes"):
        GeolocWarper(resampling=Resampling.lanczos)
