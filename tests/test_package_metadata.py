"""Tests for public package version metadata."""

from importlib.metadata import version
import re

from packaging.version import Version

import biosigpy


def test_public_version_matches_installed_metadata() -> None:
    assert Version(biosigpy.__version__) == Version(version("biosigpy"))


def test_public_version_is_stable_or_beta_semver() -> None:
    assert re.fullmatch(
        r"(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)(-beta)?",
        biosigpy.__version__,
    )


def test_beta_spelling_matches_python_package_metadata() -> None:
    assert Version("1.0.0-beta") == Version("1.0.0b0")
    assert Version("1.0.0-beta") < Version("1.0.0")
