"""Tests for VersionFile container behavior."""

import pytest
from grz_common.models.version import VersionFile, VersionInfo

def test_versionfile_accepts_policy_list():
    vf = VersionFile(
        grzcli_version=[
            {
                "minimal_version": "1.6.0",
                "recommended_version": "1.6",
                "max_version": "1.6.0",
                "enforced_from": "2026-03-01",
            }
        ]
    )

    assert len(vf.grzcli_version) == 1
    assert isinstance(vf.grzcli_version[0], VersionInfo)

def test_versionfile_requires_at_least_one_policy():
    with pytest.raises(ValueError):
        VersionFile(grzcli_version=[])

def test_multiple_policies_parsed():
    vf = VersionFile(
        grzcli_version=[
            {
                "minimal_version": "1.6.0",
                "recommended_version": "1.6",
                "max_version": "1.6.0",
                "enforced_from": "2026-03-01",
            },
            {
                "minimal_version": "1.7.0",
                "recommended_version": "1.7",
                "max_version": "1.7.0",
                "enforced_from": "2026-06-01",
            },
        ]
    )

    assert len(vf.grzcli_version) == 2

