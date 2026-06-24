"""Tests for the detector-temperature line-timing model.

These cover the timestamp-based mapping of telemetry samples to along-track
lines (the exact model) and its uniform-spacing fallback, using small synthetic
fixtures only — no proprietary data.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from spp.calibration.l1b_calibrator import L1BCalibrator
from spp.core.acquisition import Acquisition, ImagerConfiguration
from spp.core.band import Band
from spp.core.calibration_parameters import (
    CalibrationParameters,
    RadiometricCalibration,
)
from spp.readers.package_reader import PackageReader


def make_imager_config(line_period: int = 500) -> ImagerConfiguration:
    return ImagerConfiguration(
        line_period_us=line_period,
        spectral_bands=1,
        band_setup=[4],
        band_start_row=[0],
        band_cwl=[665.0],
        scan_direction=1,
        binning_factor=0,
    )


def make_acquisition(temps, sample_lines=None) -> Acquisition:
    return Acquisition(
        scene_id="synthetic",
        bands={},
        calibration=CalibrationParameters(radiometric=[]),
        imager_config=make_imager_config(),
        sensor_temperatures=np.asarray(temps, dtype=np.float64),
        temperature_sample_lines=sample_lines,
    )


# -- reader timing helpers --------------------------------------------------


def test_iso_to_epoch_ms():
    assert PackageReader._iso_to_epoch_ms("1970-01-01T00:00:00Z") == 0.0
    assert PackageReader._iso_to_epoch_ms("1970-01-01T00:00:01Z") == 1000.0
    assert PackageReader._iso_to_epoch_ms("not-a-date") is None


def test_time_sync_anchor_picks_timeformat_entry():
    session = {
        "TimeSync": [
            {"ImagerTime": 100, "PPS": 1},  # bare PPS, ignored
            {"TimeFormat": 1, "ImagerTime": 200, "PlatformTime": 5000},
        ]
    }
    assert PackageReader._time_sync_anchor(session) == (200.0, 5000.0)
    assert PackageReader._time_sync_anchor({"TimeSync": []}) is None


def test_time_sync_offset_from_properties():
    stac = {"properties": {"vendor:calibration:time_sync_offset": "-0.5"}}
    assert PackageReader._time_sync_offset(stac) == -0.5
    assert PackageReader._time_sync_offset({"properties": {}}) == 0.0
    assert PackageReader._time_sync_offset(None) == 0.0


def test_temperature_sample_lines_exact(tmp_path):
    reader = PackageReader(tmp_path)
    cfg = make_imager_config(line_period=500)
    # Anchor: ImagerTime 1_000_000 µs <-> PlatformTime 1000 ms (epoch).
    session = {"TimeSync": [{"TimeFormat": 1, "ImagerTime": 1_000_000, "PlatformTime": 1000}]}
    # Acquisition start = 1000 ms epoch == anchor, so t_line0 = 1_000_000 µs.
    acquired_at = "1970-01-01T00:00:01Z"
    sample_times = np.array([1_000_000.0, 1_000_500.0, 1_001_000.0])
    lines = reader._temperature_sample_lines(session, cfg, acquired_at, sample_times, 0.0)
    np.testing.assert_allclose(lines, [0.0, 1.0, 2.0])


def test_temperature_sample_lines_none_when_timing_missing(tmp_path):
    reader = PackageReader(tmp_path)
    cfg = make_imager_config()
    anchor = {"TimeSync": [{"TimeFormat": 1, "ImagerTime": 0, "PlatformTime": 0}]}
    # No anchor.
    assert reader._temperature_sample_lines({"TimeSync": []}, cfg, "1970-01-01T00:00:01Z", np.array([1.0]), 0.0) is None
    # No sample timestamps.
    assert reader._temperature_sample_lines(anchor, cfg, "1970-01-01T00:00:01Z", np.array([]), 0.0) is None
    # No acquisition start.
    assert reader._temperature_sample_lines(anchor, cfg, None, np.array([1.0]), 0.0) is None


def test_temperature_samples_drop_timing_if_incomplete():
    # One sample lacks ImagerTime -> timestamps discarded, temps kept.
    session = {
        "ImagerTelemetry": [
            {"ImagerTime": 10, "SensorTemperature": -1},
            {"SensorTemperature": 2},
        ]
    }
    temps, times = PackageReader._temperature_samples(session)
    np.testing.assert_allclose(temps, [-1.0, 2.0])
    assert times.size == 0


def test_optional_float_treats_missing_and_null_as_default():
    # Missing key and explicit JSON null both mean "not specified".
    assert PackageReader._optional_float(None) == 0.0
    assert PackageReader._optional_float(None, default=1.5) == 1.5
    # A real value is coerced to float; zero is preserved, not treated as null.
    assert PackageReader._optional_float(2.5) == 2.5
    assert PackageReader._optional_float(0) == 0.0
    # A malformed (non-numeric) value still fails loudly.
    with pytest.raises((TypeError, ValueError)):
        PackageReader._optional_float("not-a-number")


# -- calibrator interpolation ----------------------------------------------


def test_interpolate_uniform_fallback():
    acq = make_acquisition([0.0, 10.0])  # no sample_lines -> uniform
    profile = L1BCalibrator(acq)._interpolate_temperature(11)
    np.testing.assert_allclose(profile, np.linspace(0.0, 10.0, 11))


def test_interpolate_timestamped_positions():
    # Samples sit at lines 0 and 5; lines beyond 5 clamp to the last value.
    acq = make_acquisition([0.0, 10.0], sample_lines=np.array([0.0, 5.0]))
    profile = L1BCalibrator(acq)._interpolate_temperature(11)
    expected = np.concatenate([np.linspace(0.0, 10.0, 6), np.full(5, 10.0)])
    np.testing.assert_allclose(profile, expected)


def test_interpolate_single_sample_is_constant():
    acq = make_acquisition([5.0])
    profile = L1BCalibrator(acq)._interpolate_temperature(7)
    np.testing.assert_allclose(profile, np.full(7, 5.0))


# -- calibrator core: DN -> radiance ---------------------------------------


def make_calibrate_acquisition():
    """A 4-line x 3-column single-band acquisition with a known calibration.

    The temperature profile is pinned to ``[0, 1, 2, 3]`` (two samples at lines
    0 and 3) so the temperature-dependent darkfield is exact and hand-checkable.
    """
    width, height = 3, 4
    rc = RadiometricCalibration(
        band="R",
        start_row=0,
        tdi=4,
        absolute=2.0,
        absolute_offset=100.0,
        non_uniformity=np.array([1.0, 2.0, 0.5]),
        thermal_intercept=np.array([10.0, 0.0, 4.0]),
        thermal_gradient=np.array([1.0, 0.0, 2.0]),
    )
    band = Band(
        name="R",
        band_id=0,
        cwl_nm=665.0,
        path=Path("R.tif"),
        width=width,
        height=height,
        dtype="uint16",
        nodata=0,
    )
    acq = Acquisition(
        scene_id="synthetic",
        bands={"R": band},
        calibration=CalibrationParameters(radiometric=[rc]),
        imager_config=make_imager_config(),
        sensor_temperatures=np.array([0.0, 3.0]),
        temperature_sample_lines=np.array([0.0, 3.0]),
    )
    return acq, rc


def _expected_radiance(rc, temperature, dn):
    """Reference implementation of the documented calibration formula."""
    darkfield = rc.thermal_intercept[None, :] + (
        rc.thermal_gradient[None, :] * temperature[:, None]
    )
    relative = (dn - darkfield) * rc.non_uniformity[None, :]
    return rc.absolute * relative + rc.absolute_offset


def test_calibrate_applies_documented_formula():
    acq, rc = make_calibrate_acquisition()
    cal = L1BCalibrator(acq)
    dn = np.array(
        [[100, 100, 100],
         [200, 150, 120],
         [300, 250, 130],
         [400, 350, 140]],
        dtype=np.uint16,
    )
    out = cal.calibrate("R", dn, line_start=0)

    temperature = cal.per_line_temperature("R")  # [0, 1, 2, 3]
    expected = _expected_radiance(rc, temperature, dn.astype(np.float64))
    np.testing.assert_allclose(out, expected, rtol=1e-5)
    assert out.dtype == np.float32

    # Two pixels pinned by hand against the documented equation:
    #   radiance = absolute * (DN - (intercept + gradient*T)) * non_uniformity
    #            + absolute_offset
    # line 1, col 1: T=1, darkfield = 0 + 0*1 = 0,
    #   relative = (150 - 0)*2 = 300, radiance = 2*300 + 100 = 700.
    assert out[1, 1] == pytest.approx(700.0)
    # line 2, col 0 (temperature-dependent): T=2, darkfield = 10 + 1*2 = 12,
    #   relative = (300 - 12)*1 = 288, radiance = 2*288 + 100 = 676.
    assert out[2, 0] == pytest.approx(676.0)


def test_calibrate_windowed_matches_full_frame():
    # The property that makes windowed streaming exact: line_start indexes the
    # per-line temperature profile, so any split reproduces the full frame.
    acq, _ = make_calibrate_acquisition()
    cal = L1BCalibrator(acq)
    dn = (np.arange(12, dtype=np.uint16) + 1).reshape(4, 3)  # 1..12, no NoData

    full = cal.calibrate("R", dn, line_start=0)
    top = cal.calibrate("R", dn[:2], line_start=0)
    bottom = cal.calibrate("R", dn[2:], line_start=2)

    np.testing.assert_allclose(np.vstack([top, bottom]), full)


def test_calibrate_marks_nodata_pixels():
    # DN equal to the band's NoData (0) becomes NaN; other pixels stay finite.
    acq, _ = make_calibrate_acquisition()
    cal = L1BCalibrator(acq)
    dn = np.array(
        [[0, 5, 5],
         [5, 5, 0],
         [5, 5, 5],
         [5, 5, 5]],
        dtype=np.uint16,
    )
    out = cal.calibrate("R", dn, line_start=0)

    assert np.isnan(out[0, 0])
    assert np.isnan(out[1, 2])
    assert np.isfinite(out[0, 1])


def test_calibrate_rejects_bad_shape_and_out_of_bounds():
    acq, _ = make_calibrate_acquisition()
    cal = L1BCalibrator(acq)
    with pytest.raises(ValueError):
        cal.calibrate("R", np.zeros((4,), dtype=np.uint16))  # not 2-D
    with pytest.raises(ValueError):
        cal.calibrate("R", np.zeros((4, 2), dtype=np.uint16))  # wrong width
    with pytest.raises(ValueError):
        # block [3:5] exceeds the band height of 4 lines
        cal.calibrate("R", np.zeros((2, 3), dtype=np.uint16), line_start=3)
