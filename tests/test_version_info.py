"""Tests for VersionInfo schema and validation."""

from datetime import UTC, datetime

import pytest
from grz_common.models.version import VersionInfo
from packaging.version import Version


def test_versioninfo_full_parsing():
    vi = VersionInfo(
        minimal_version="1.6.0",
        recommended_version="1.6",
        max_version="1.6.0",
        enforced_from="2026-03-01T00:00:00",
    )

    assert vi.minimal_version == Version("1.6.0")
    assert vi.recommended_version == Version("1.6")
    assert vi.max_version == Version("1.6.0")
    assert vi.enforced_from == datetime(2026, 3, 1, tzinfo=UTC)


def test_versioninfo_minimal_only():
    """Only required fields should be sufficient."""
    vi = VersionInfo(
        minimal_version="1.6.0",
        enforced_from="2026-03-01T00:00:00",
    )

    assert vi.minimal_version == Version("1.6.0")
    assert vi.recommended_version is None
    assert vi.max_version is None
    assert vi.enforced_from == datetime(2026, 3, 1, tzinfo=UTC)


def test_versioninfo_with_recommended_only():
    """max_version can be absent when recommended_version is present."""
    vi = VersionInfo(
        minimal_version="1.6.0",
        recommended_version="1.6.1",
        enforced_from="2026-03-01T00:00:00",
    )

    assert vi.recommended_version == Version("1.6.1")
    assert vi.max_version is None


def test_versioninfo_with_max_only():
    """recommended_version can be absent when max_version is present."""
    vi = VersionInfo(
        minimal_version="1.6.0",
        max_version="1.6.1",
        enforced_from="2026-03-01T00:00:00",
    )

    assert vi.recommended_version is None
    assert vi.max_version == Version("1.6.1")


def test_recommended_cannot_be_lower_than_minimal_when_provided():
    """Version order check still applies when recommended_version is present."""
    with pytest.raises(ValueError):
        VersionInfo(
            minimal_version="1.5.0",
            recommended_version="1.4.0",
            enforced_from="2026-01-01T00:00:00",
        )


def test_max_cannot_be_lower_than_recommended_when_both_provided():
    """Version order check still applies when both optional fields are present."""
    with pytest.raises(ValueError):
        VersionInfo(
            minimal_version="1.5.0",
            recommended_version="1.5.1",
            max_version="1.5.0",
            enforced_from="2026-01-01T00:00:00",
        )
