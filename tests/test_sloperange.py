"""Conformance and Python-native tests for ecg.sloperange."""

import numpy as np
import pytest

from biosigpy import sloperange as top_level_sloperange
from biosigpy.ecg import sloperange as ecg_sloperange
from biosigpy.ecg.sloperange import sloperange
from conformance import (
    assert_expected_error,
    assert_expected_outputs,
    case_id,
    cases_for_specification,
    is_expected_error,
    load_input,
)


@pytest.mark.parametrize(
    "case_definition",
    cases_for_specification("ecg.sloperange"),
    ids=case_id,
)
def test_conformance(case_definition: dict[str, object]) -> None:
    decg = load_input(case_definition, "decg")
    r_wave_times = load_input(case_definition, "r_wave_times")
    sampling_frequency = load_input(case_definition, "sampling_frequency")
    run_case = lambda: sloperange(decg, r_wave_times, sampling_frequency)
    if is_expected_error(case_definition):
        assert_expected_error(run_case, case_definition)
        return

    edr = run_case()
    assert edr.ndim == 1
    assert edr.shape == np.asarray(r_wave_times).reshape(-1).shape
    assert_expected_outputs({"edr": edr}, case_definition)


def test_row_and_column_vectors_have_identical_semantics() -> None:
    decg = np.zeros(40)
    decg[[9, 19, 29]] = [3.0, 6.0, 2.0]
    decg[[13, 23, 33]] = [-2.0, -1.0, -4.0]
    r_wave_times = np.array([0.1, 0.2, 0.3])

    expected = sloperange(decg, r_wave_times, 100.0)
    actual = sloperange(
        decg.reshape(1, -1), r_wave_times.reshape(-1, 1), 100.0
    )

    np.testing.assert_array_equal(actual, expected)
    assert actual.shape == (3,)


def test_half_sample_positions_round_like_biosiglib() -> None:
    decg = np.zeros(30)
    decg[13] = 5.0
    decg[15] = -2.0

    edr = sloperange(decg, [0.105], 100.0)

    np.testing.assert_array_equal(edr, np.array([7.0]))


@pytest.mark.parametrize(
    ("arguments", "exception_type"),
    [
        (([0.0], [0.0], 100.0), ValueError),
        (([0.0, np.nan], [0.0], 100.0), ValueError),
        (([0.0, np.inf], [0.0], 100.0), ValueError),
        (([0.0, 1.0], [], 100.0), ValueError),
        (([0.0, 1.0], [np.nan], 100.0), ValueError),
        (([0.0, 1.0], [np.inf], 100.0), ValueError),
        (([[0.0, 1.0], [2.0, 3.0]], [0.0], 100.0), ValueError),
        (([0.0, 1.0], [[0.0, 0.1], [0.2, 0.3]], 100.0), ValueError),
        (([0.0, 1.0], [0.0], 0.0), ValueError),
        (([0.0, 1.0], [0.0], [100.0]), ValueError),
        (([0.0, 1.0], [0.0], "100"), TypeError),
    ],
)
def test_python_validation(
    arguments: tuple[object, object, object],
    exception_type: type[Exception],
) -> None:
    with pytest.raises(exception_type):
        sloperange(*arguments)


def test_public_reexports_reference_the_canonical_function() -> None:
    assert ecg_sloperange is sloperange
    assert top_level_sloperange is sloperange
