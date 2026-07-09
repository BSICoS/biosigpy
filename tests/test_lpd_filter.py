"""Conformance tests for tools.lpd_filter."""

import numpy as np
import pytest

from biosigpy.tools.lpd_filter import lpd_filter
from conformance import (
    assert_expected_error,
    assert_expected_outputs,
    load_case,
    load_input,
)


VALID_CASES = [
    "tools.lpd_filter.fs256_stop12_order4_coefficients",
    "tools.lpd_filter.explicit_pass_frequency_order4_coefficients",
]

EXPECTED_ERROR_CASES = [
    "tools.lpd_filter.invalid_pass_frequency_not_less_than_stop",
    "tools.lpd_filter.invalid_stop_frequency_at_nyquist",
]


@pytest.mark.parametrize("case_id", VALID_CASES)
def test_positive_conformance(case_id: str) -> None:
    case_definition = load_case(case_id)
    parameters = case_definition["parameters"]

    filter_coefficients, delay = lpd_filter(
        load_input(case_definition, "sampling_frequency"),
        load_input(case_definition, "stop_frequency"),
        parameters.get("pass_frequency"),
        parameters.get("order"),
    )

    assert np.issubdtype(filter_coefficients.dtype, np.number)
    assert filter_coefficients.ndim == 1
    assert_expected_outputs(
        {"filter_coefficients": filter_coefficients, "delay": delay},
        case_definition,
    )


@pytest.mark.parametrize("case_id", EXPECTED_ERROR_CASES)
def test_expected_error_conformance(case_id: str) -> None:
    case_definition = load_case(case_id)
    parameters = case_definition["parameters"]

    assert_expected_error(
        lambda: lpd_filter(
            load_input(case_definition, "sampling_frequency"),
            load_input(case_definition, "stop_frequency"),
            parameters.get("pass_frequency"),
            parameters.get("order"),
        ),
        case_definition,
    )


@pytest.mark.parametrize(
    ("kwargs", "exception_type"),
    [
        ({"sampling_frequency": 0.0}, ValueError),
        ({"stop_frequency": 128.0}, ValueError),
        ({"pass_frequency": 12.0}, ValueError),
        ({"order": 0}, ValueError),
        ({"order": 4.5}, TypeError),
    ],
)
def test_python_validation(
    kwargs: dict[str, object], exception_type: type[Exception]
) -> None:
    arguments = {
        "sampling_frequency": 256.0,
        "stop_frequency": 12.0,
        "pass_frequency": 10.0,
        "order": 4,
    }
    arguments.update(kwargs)

    with pytest.raises(exception_type):
        lpd_filter(**arguments)


def test_odd_explicit_order_rounds_up_to_even() -> None:
    coefficients_odd, delay_odd = lpd_filter(256.0, 12.0, 10.0, order=5)
    coefficients_even, delay_even = lpd_filter(256.0, 12.0, 10.0, order=6)

    np.testing.assert_allclose(coefficients_odd, coefficients_even, rtol=0.0)
    assert delay_odd == delay_even == 3.0


def test_automatic_order_is_explicitly_unsupported() -> None:
    with pytest.raises(NotImplementedError):
        lpd_filter(256.0, 12.0)
