"""Pytest integration for complete Biosiglib conformance-case coverage."""

from collections.abc import Mapping
from typing import Any

import matplotlib
import pytest

from conformance import assert_complete_case_coverage, discover_cases


# Render plotting and debug paths without opening GUI windows during tests.
matplotlib.use("Agg", force=True)


def pytest_collection_modifyitems(
    config: pytest.Config, items: list[pytest.Item]
) -> None:
    """Fail default full-suite collection when a discovered case is uncollected."""

    if config.invocation_params.args:
        return

    discovered_case_ids = {case["id"] for case in discover_cases()}
    collected_case_ids: set[str] = set()
    for item in items:
        callspec = getattr(item, "callspec", None)
        if callspec is None:
            continue
        for parameter in callspec.params.values():
            if _is_conformance_case(parameter):
                collected_case_ids.add(str(parameter["id"]))

    try:
        assert_complete_case_coverage(discovered_case_ids, collected_case_ids)
    except RuntimeError as error:
        raise pytest.UsageError(str(error)) from error


def _is_conformance_case(parameter: Any) -> bool:
    return (
        isinstance(parameter, Mapping)
        and "id" in parameter
        and "specification_id" in parameter
        and ("expected_outputs" in parameter or "expected_error" in parameter)
    )
