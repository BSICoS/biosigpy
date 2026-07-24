"""Conformance tests for ecg.pantompkins."""

import importlib

import numpy as np
import pytest

from biosigpy.ecg.pantompkins import pantompkins
from conformance import (
    assert_expected_error,
    assert_expected_outputs,
    load_case,
    load_input,
)


PANTOMPKINS_MODULE = importlib.import_module("biosigpy.ecg.pantompkins")


EXPECTED_ERROR_CASES = [
    "ecg.pantompkins.invalid_sampling_frequency_non_positive",
    "ecg.pantompkins.invalid_sampling_frequency_vector",
    "ecg.pantompkins.invalid_sampling_frequency_non_numeric",
    "ecg.pantompkins.invalid_ecg_matrix",
    "ecg.pantompkins.invalid_ecg_non_numeric",
]


def test_positive_conformance() -> None:
    case_definition = load_case("ecg.pantompkins.medicom_mtd_r_wave_times")
    ecg = load_input(case_definition, "ecg")
    sampling_frequency = load_input(case_definition, "sampling_frequency")

    outputs = pantompkins(ecg, sampling_frequency)

    assert set(outputs) == {
        "r_wave_times",
        "ecg_filtered",
        "decg_squared",
        "decg_envelope",
    }
    for output_id in ("ecg_filtered", "decg_squared", "decg_envelope"):
        output = outputs[output_id]
        assert np.issubdtype(output.dtype, np.number)
        assert output.ndim == 1
        assert output.shape == ecg.shape
    assert_expected_outputs(outputs, case_definition)


@pytest.mark.parametrize("case_id", EXPECTED_ERROR_CASES)
def test_expected_error_conformance(case_id: str) -> None:
    case_definition = load_case(case_id)
    ecg = load_input(case_definition, "ecg")
    sampling_frequency = load_input(case_definition, "sampling_frequency")

    assert_expected_error(
        lambda: pantompkins(ecg, sampling_frequency), case_definition
    )


def test_filters_nan_refined_detections_before_r_wave_times(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fake_snap_to_peak(
        ecg: np.ndarray, detections: np.ndarray, window_size: float
    ) -> np.ndarray:
        np.testing.assert_array_equal(detections, np.array([2.0, 5.0]))
        assert window_size == 3.0
        return np.array([np.nan, 5.0])

    def fake_lpd_filter(
        sampling_frequency: float,
        stop_frequency: float,
        pass_frequency: float | None = None,
        order: int | None = None,
    ) -> tuple[np.ndarray, float]:
        assert sampling_frequency == 40.0
        assert stop_frequency == 12.0
        assert pass_frequency is None
        assert order == 4
        return np.array([1.0]), 0.0

    def fake_nan_filtfilt(
        b: np.ndarray,
        a: np.ndarray,
        ecg: np.ndarray,
        max_gap: int = 0,
    ) -> np.ndarray:
        np.testing.assert_array_equal(b, np.array([1.0]))
        np.testing.assert_array_equal(a, np.array([1.0]))
        assert max_gap == 0
        return np.asarray(ecg, dtype=np.float64)

    def fake_nan_filter(
        b: np.ndarray,
        a: np.ndarray,
        ecg: np.ndarray,
        max_gap: int = 0,
    ) -> np.ndarray:
        np.testing.assert_array_equal(b, np.array([1.0]))
        np.testing.assert_array_equal(a, np.array([1.0]))
        assert max_gap == 0
        return np.asarray(ecg, dtype=np.float64)

    monkeypatch.setattr(
        PANTOMPKINS_MODULE.signal,
        "butter",
        lambda *args, **kwargs: (np.array([1.0]), np.array([1.0])),
    )
    monkeypatch.setattr(PANTOMPKINS_MODULE, "nan_filtfilt", fake_nan_filtfilt)
    monkeypatch.setattr(PANTOMPKINS_MODULE, "nan_filter", fake_nan_filter)
    monkeypatch.setattr(
        PANTOMPKINS_MODULE.signal,
        "find_peaks",
        lambda envelope, distance: (np.array([1, 4]), {}),
    )
    monkeypatch.setattr(PANTOMPKINS_MODULE, "lpd_filter", fake_lpd_filter)
    monkeypatch.setattr(PANTOMPKINS_MODULE, "snap_to_peak", fake_snap_to_peak)

    outputs = pantompkins(
        np.arange(10.0),
        40.0,
        snap_to_peak_window_size=3.0,
    )

    np.testing.assert_array_equal(outputs["r_wave_times"], np.array([0.1]))
    np.testing.assert_array_equal(outputs["decg_squared"], np.arange(10.0) ** 2)


def test_nan_gap_is_preserved_as_a_hard_filtering_boundary() -> None:
    sampling_frequency = 100.0
    time = np.arange(0.0, 12.0, 1.0 / sampling_frequency)
    ecg = np.sin(2.0 * np.pi * 1.2 * time)
    gap = slice(550, 650)
    ecg[gap] = np.nan

    outputs = pantompkins(ecg, sampling_frequency)

    outside_gap = np.ones(ecg.shape, dtype=bool)
    outside_gap[gap] = False

    assert np.all(np.isnan(outputs["ecg_filtered"][gap]))
    assert np.all(np.isnan(outputs["decg_squared"][gap]))
    assert np.all(np.isfinite(outputs["ecg_filtered"][outside_gap]))
    assert np.all(np.isfinite(outputs["decg_squared"][outside_gap]))
    assert not np.any(np.isinf(outputs["ecg_filtered"]))
    assert not np.any(np.isinf(outputs["decg_squared"]))
