"""Conformance tests for ecg.pantompkins."""

import numpy as np
import pytest

from biosigpy.ecg.pantompkins import pantompkins
from conformance import (
    assert_expected_error,
    assert_expected_outputs,
    load_case,
    load_input,
)


EXPECTED_ERROR_CASES = [
    "ecg.pantompkins.invalid_sampling_frequency_non_positive",
    "ecg.pantompkins.invalid_sampling_frequency_vector",
    "ecg.pantompkins.invalid_sampling_frequency_non_numeric",
    "ecg.pantompkins.invalid_ecg_matrix",
    "ecg.pantompkins.invalid_ecg_non_numeric",
]


def test_positive_conformance() -> None:
    case_definition = load_case("ecg.pantompkins.edr_signals_001")
    ecg = load_input(case_definition, "ecg")
    sampling_frequency = load_input(case_definition, "sampling_frequency")

    outputs = pantompkins(ecg, sampling_frequency)

    assert set(outputs) == {
        "r_peak_times",
        "ecg_filtered",
        "decg",
        "decg_envelope",
    }
    for output_id in ("ecg_filtered", "decg", "decg_envelope"):
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
