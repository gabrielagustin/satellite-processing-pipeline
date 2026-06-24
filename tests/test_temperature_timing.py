"""Tests for the detector-temperature line-timing model.

These cover the timestamp-based mapping of telemetry samples to along-track
lines (the exact model) and its uniform-spacing fallback, using small synthetic
fixtures only — no proprietary data.
"""

from __future__ import annotations

import numpy as np
import pytest

from spp.calibration.l1b_calibrator import L1BCalibrator
from spp.core.acquisition import Acquisition, ImagerConfiguration
from spp.core.calibration_parameters import CalibrationParameters
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
