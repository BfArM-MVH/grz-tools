from datetime import date
from unittest.mock import patch

import pytest
from grz_cli.utils.version_check import check_version_and_exit_if_needed
from packaging.version import Version


class DummyVersionFile:
    def __init__(self, minimal_version, recommended_version, max_version, enforced_from):
        self.schema_version = 1
        self.minimal_version = Version(minimal_version)
        self.recommended_version = Version(recommended_version)
        self.max_version = Version(max_version)
        self.enforced_from = enforced_from


def test_too_old_after_enforcement():
    """The version is old --> sys.exit"""
    # Use PAST date for enforcement
    vf = DummyVersionFile("1.4.0", "1.5.0", "1.5.1", date(2020, 1, 1))

    with (
        patch("grz_cli.utils.version_check.VersionFile.from_s3", return_value=vf),
        patch("grz_cli.utils.version_check.version", return_value="1.3.0"),
        patch("sys.exit", side_effect=SystemExit) as mock_exit,
    ):
        with pytest.raises(SystemExit):
            check_version_and_exit_if_needed(None)

        mock_exit.assert_called_once_with(1)


def test_too_old_before_enforcement(caplog):
    # Use a FUTURE date
    vf = DummyVersionFile("1.4.0", "1.5.0", "1.5.1", date(2099, 1, 1))

    with (
        patch("grz_cli.utils.version_check.VersionFile.from_s3", return_value=vf),
        patch("grz_cli.utils.version_check.version", return_value="1.3.0"),
    ):
        check_version_and_exit_if_needed(None)

    assert "unsupported" in caplog.text
    assert "2099-01-01" in caplog.text
    assert "1.4.0" in caplog.text


def test_outdated_but_supported(caplog):
    """The version is behind the recommended version --> raise warning"""
    vf = DummyVersionFile("1.4.0", "1.5.0", "1.5.1", date(2026, 1, 1))

    with (
        patch("grz_cli.utils.version_check.VersionFile.from_s3", return_value=vf),
        patch("grz_cli.utils.version_check.version", return_value="1.4.9"),
    ):
        check_version_and_exit_if_needed(None)
        assert "recommended version" in caplog.text


def test_version_ok(caplog):
    """The version is fine"""
    vf = DummyVersionFile("1.4.0", "1.5.0", "1.5.1", date(2026, 1, 1))

    with (
        patch("grz_cli.utils.version_check.VersionFile.from_s3", return_value=vf),
        patch("grz_cli.utils.version_check.version", return_value="1.5.0"),
    ):
        check_version_and_exit_if_needed(None)
        assert "up to date" in caplog.text or "supported and tested range" in caplog.text


def test_newer_than_tested(caplog):
    vf = DummyVersionFile("1.4.0", "1.5.0", "1.5.0", date(2026, 1, 1))  # max_tested = 1.5.0

    with (
        patch("grz_cli.utils.version_check.VersionFile.from_s3", return_value=vf),
        patch("grz_cli.utils.version_check.version", return_value="1.6.0"),
    ):
        check_version_and_exit_if_needed(None)
        assert "newer than the latest tested version" in caplog.text
