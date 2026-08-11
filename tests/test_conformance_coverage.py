"""Regression tests for shared conformance harness behavior."""

from types import SimpleNamespace

import pytest

from conformance import assert_complete_case_coverage, assert_expected_warnings


def test_newly_discovered_uncollected_case_fails() -> None:
    discovered = {
        "hrv.tdmetrics.existing_case",
        "hrv.tdmetrics.new_case",
    }
    collected = {"hrv.tdmetrics.existing_case"}

    with pytest.raises(RuntimeError, match="hrv.tdmetrics.new_case"):
        assert_complete_case_coverage(discovered, collected)


def test_expected_warnings_compare_ids_and_affected_sets() -> None:
    case_definition = {
        "expected_warnings": [
            {"id": "example_warning", "affected_ids": ["first", "second"]}
        ]
    }
    actual = [
        SimpleNamespace(
            warning_id="example_warning",
            affected_ids=("second", "first"),
        )
    ]

    assert_expected_warnings(actual, case_definition)


def test_duplicate_canonical_warning_id_fails() -> None:
    case_definition = {
        "expected_warnings": [
            {"id": "example_warning", "affected_ids": ["output"]}
        ]
    }
    warning = SimpleNamespace(
        warning_id="example_warning",
        affected_ids=("output",),
    )

    with pytest.raises(AssertionError, match="once per canonical id"):
        assert_expected_warnings([warning, warning], case_definition)
