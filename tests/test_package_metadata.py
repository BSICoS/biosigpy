"""Tests for public package version metadata."""

from importlib.metadata import version

import biosigpy


def test_public_version_matches_installed_metadata() -> None:
    assert biosigpy.__version__ == version("biosigpy")


def test_initial_release_is_pre_one() -> None:
    major, minor, patch = (int(part) for part in biosigpy.__version__.split("."))
    assert (major, minor, patch) >= (0, 1, 0)
    assert major == 0
