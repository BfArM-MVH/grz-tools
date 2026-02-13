"""Tests for VersionInfo schema and validation."""

from datetime import date

import pytest
from grz_common.models.version import VersionInfo
from packaging.version import Version


def test_versioninfo_full_parsing():
    vi = VersionInfo(
        minimal_version="1.6.0",
        recommended_version="1.6",
        max_version="1.6.0",
        enforced_from="2026-03-01",
    )

    assert vi.minimal_version == Version("1.6.0")
    assert vi.recommended_version == Version("1.6")
    assert vi.max_version == Version("1.6.0")
    assert vi.enforced_from == date(2026, 3, 1)


def test_recommended_cannot_be_lower_than_minimal():
    with pytest.raises(ValueError):
        VersionInfo(
            minimal_version="1.5.0",
            recommended_version="1.4.0",
            max_version="1.5.0",
            enforced_from="2026-01-01",
        )


def test_max_cannot_be_lower_than_recommended():
    with pytest.raises(ValueError):
        VersionInfo(
            minimal_version="1.5.0",
            recommended_version="1.5.1",
            max_version="1.5.0",
            enforced_from="2026-01-01",
        )
