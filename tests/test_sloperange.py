"""Conformance and Python-native tests for ecg.sloperange."""

import numpy as np
import pytest

from biosigpy import SlopeRangeResult as top_level_slope_range_result
from biosigpy import sloperange as top_level_sloperange
from biosigpy.ecg import SlopeRangeResult as ecg_slope_range_result
from biosigpy.ecg import sloperange as ecg_sloperange
from biosigpy.ecg.sloperange import SlopeRangeResult, sloperange
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

    result = run_case()
    beat_shape = np.asarray(r_wave_times).reshape(-1).shape
    signal_shape = np.asarray(decg).reshape(-1).shape
    assert result.edr.shape == beat_shape
    assert result.upslope_max_positions.shape == beat_shape
    assert result.downslope_min_positions.shape == beat_shape
    assert result.upslopes.shape == signal_shape
    assert result.downslopes.shape == signal_shape
    assert all(output.ndim == 1 for output in result)
    assert_expected_outputs(result._asdict(), case_definition)


def test_row_and_column_vectors_have_identical_semantics() -> None:
    decg = np.zeros(40)
    decg[[9, 19, 29]] = [3.0, 6.0, 2.0]
    decg[[13, 23, 33]] = [-2.0, -1.0, -4.0]
    r_wave_times = np.array([0.1, 0.2, 0.3])

    expected = sloperange(decg, r_wave_times, 100.0)
    actual = sloperange(
        decg.reshape(1, -1), r_wave_times.reshape(-1, 1), 100.0
    )

    for field_name in SlopeRangeResult._fields:
        np.testing.assert_array_equal(
            getattr(actual, field_name), getattr(expected, field_name)
        )
    assert actual.edr.shape == (3,)
    assert actual.upslopes.shape == (40,)


def test_half_sample_positions_round_like_biosiglib() -> None:
    decg = np.zeros(30)
    decg[13] = 5.0
    decg[15] = -2.0

    result = sloperange(decg, [0.105], 100.0)

    np.testing.assert_array_equal(result.edr, np.array([7.0]))
    np.testing.assert_array_equal(
        result.upslope_max_positions, np.array([13.0])
    )
    np.testing.assert_array_equal(
        result.downslope_min_positions, np.array([15.0])
    )


def test_named_result_unpacks_in_documented_order() -> None:
    result = sloperange(np.zeros(40), [0.1], 100.0)
    edr, upslopes, downslopes, upmaxpos, downminpos = result

    assert isinstance(result, SlopeRangeResult)
    assert SlopeRangeResult._fields == (
        "edr",
        "upslopes",
        "downslopes",
        "upslope_max_positions",
        "downslope_min_positions",
    )
    assert edr is result.edr
    assert upslopes is result.upslopes
    assert downslopes is result.downslopes
    assert upmaxpos is result.upslope_max_positions
    assert downminpos is result.downslope_min_positions


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
    assert ecg_slope_range_result is SlopeRangeResult
    assert top_level_slope_range_result is SlopeRangeResult
