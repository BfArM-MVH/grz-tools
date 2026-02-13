import pytest
from packaging.version import Version
from grz_cli.utils.version_check import check_version_and_exit_if_needed

from datetime import date
from unittest.mock import patch

class DummyVersionInfo:
    def __init__(self, minimal_version, recommended_version, max_version, enforced_from):
        self.minimal_version = Version(minimal_version)
        self.recommended_version = Version(recommended_version)
        self.max_version = Version(max_version)
        self.enforced_from = enforced_from

class DummyVersionFile:
    def __init__(self, policies):
        self.schema_version = 1
        self.grzcli_version = policies


def test_too_old_after_enforcement():
    """The version is old --> sys.exit"""
    policy = DummyVersionInfo("1.4.0", "1.5.0", "1.5.1", date(2020, 1, 1))
    vf = DummyVersionFile([policy])

    with (
        patch("grz_cli.utils.version_check.VersionFile.from_s3", return_value=vf),
        patch("grz_cli.utils.version_check.version", return_value="1.3.0"),
        patch("sys.exit", side_effect=SystemExit) as mock_exit,
    ):
        with pytest.raises(SystemExit):
            check_version_and_exit_if_needed(None)

        mock_exit.assert_called_once_with(1)

def test_outdated_but_supported(caplog):
    """The version is behind the recommended version --> raise warning"""
    policy = DummyVersionInfo("1.4.0", "1.5.0", "1.5.1", date(2020, 1, 1))
    vf = DummyVersionFile([policy])

    with (
        patch("grz_cli.utils.version_check.VersionFile.from_s3", return_value=vf),
        patch("grz_cli.utils.version_check.version", return_value="1.4.9"),
    ):
        check_version_and_exit_if_needed(None)

    assert "recommended version" in caplog.text


def test_version_ok(caplog):
    """The version is fine"""
    policy = DummyVersionInfo("1.4.0", "1.5.0", "1.5.1", date(2020, 1, 1))
    vf = DummyVersionFile([policy])

    with (
        patch("grz_cli.utils.version_check.VersionFile.from_s3", return_value=vf),
        patch("grz_cli.utils.version_check.version", return_value="1.5.0"),
    ):
        check_version_and_exit_if_needed(None)

    assert "supported and tested range" in caplog.text

def test_newer_than_tested():
    """The version is more recent than the tested version --> sys.exit"""
    policy = DummyVersionInfo("1.4.0", "1.5.0", "1.5.0", date(2020, 1, 1))
    vf = DummyVersionFile([policy])

    with (
        patch("grz_cli.utils.version_check.VersionFile.from_s3", return_value=vf),
        patch("grz_cli.utils.version_check.version", return_value="1.6.0"),
        patch("sys.exit", side_effect=SystemExit) as mock_exit,
    ):
        with pytest.raises(SystemExit):
            check_version_and_exit_if_needed(None)

        mock_exit.assert_called_once_with(1)

def test_selects_latest_active_policy():
    """Implicitly verifies the older policy was selected"""
    past_policy = DummyVersionInfo("1.0.0", "1.1.0", "2.0.0", date(2020, 1, 1))
    future_policy = DummyVersionInfo("2.0.0", "2.1.0", "3.0.0", date(2099, 1, 1))

    vf = DummyVersionFile([past_policy, future_policy])

    with (
        patch("grz_cli.utils.version_check.VersionFile.from_s3", return_value=vf),
        patch("grz_cli.utils.version_check.version", return_value="1.5.0"),
    ):
        check_version_and_exit_if_needed(None)

