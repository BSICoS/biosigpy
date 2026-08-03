"""Regression tests for the conformance-case coverage gate."""

import pytest

from conformance import assert_complete_case_coverage


def test_newly_discovered_uncollected_case_fails() -> None:
    discovered = {
        "hrv.tdmetrics.existing_case",
        "hrv.tdmetrics.new_case",
    }
    collected = {"hrv.tdmetrics.existing_case"}

    with pytest.raises(RuntimeError, match="hrv.tdmetrics.new_case"):
        assert_complete_case_coverage(discovered, collected)
