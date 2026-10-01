"""Tests for the private key that ``grzctl decrypt`` takes with ``--archive`` or ``--private-key-path``.

``test_inbox_tracking`` covers the default, the key of the inbox that the submission came from.
"""

from pathlib import Path
from unittest.mock import MagicMock, patch

import click.testing
import crypt4gh.keys
import crypt4gh.keys.c4gh
import grzctl.cli
import pytest
import yaml
from grz_common.exceptions import ConfigurationError
from grz_db.models.submission import FailureReasonEnum, SubmissionDb, SubmissionStateEnum, SubmissionStateLog
from grzctl.models.config import GrzctlConfig

SUBMITTER_ID = "260914050"
SUBMISSION_ID = "260914050_2025-09-15_c64603a7"
PASSPHRASE = "decrypt-key-passphrase"


@pytest.fixture
def db(migrated_database_config: GrzctlConfig) -> SubmissionDb:
    """A SubmissionDb on the migrated database, without an author (adds submissions and reads states only)."""
    return SubmissionDb(db_url=migrated_database_config.db.database_url, author=None)


def _generate_key(tmp_path: Path, name: str, passphrase: str | None = None) -> Path:
    """Write a crypt4gh key pair under *tmp_path* and return the path of its private key."""
    private_key_path = tmp_path / f"{name}.sec"
    crypt4gh.keys.c4gh.generate(
        private_key_path, tmp_path / f"{name}.pub", passphrase.encode() if passphrase else None, comment=None
    )
    return private_key_path


def _raw_key(private_key_path: Path, passphrase: str | None = None) -> bytes:
    return crypt4gh.keys.get_private_key(private_key_path, lambda: passphrase)


def _config_path(config: GrzctlConfig, tmp_path: Path, inboxes=("inbox",), **archive_keys: dict) -> Path:
    """A config file whose submitter :data:`SUBMITTER_ID` has the given inboxes, all with the key of the DB author.

    Each keyword names an archive, and its value holds the private key fields of that archive.
    """
    data = config.model_dump(mode="json", exclude_none=True, context={"reveal_secrets": True})
    author_key_path = data["db"]["author"]["private_key_path"]
    data["leistungserbringer"] = {
        SUBMITTER_ID: {"inbox_buckets": {name: {"private_key_path": author_key_path} for name in inboxes}}
    }
    for archive, key_fields in archive_keys.items():
        data["archives"][archive].update(key_fields)
    config_path = tmp_path / "config.decrypt.yaml"
    config_path.write_text(yaml.safe_dump(data))
    return config_path


def _decrypt_worker() -> MagicMock:
    encrypted = MagicMock()
    encrypted.submission_id = SUBMISSION_ID
    encrypted.metadata.content.submission.submitter_id = SUBMITTER_ID
    worker = MagicMock()
    worker.parse_encrypted_submission.return_value = encrypted
    return worker


def _decrypt(config_path: Path, tmp_path: Path, worker: MagicMock, *args: str) -> click.testing.Result:
    """Run ``grzctl decrypt`` with the mocked *worker* and the extra *args*."""
    submission_dir = tmp_path / "submission"
    submission_dir.mkdir(exist_ok=True)
    with patch("grzctl.commands.decrypt.Worker", return_value=worker):
        return click.testing.CliRunner().invoke(
            grzctl.cli.build_cli(),
            ["--config", str(config_path), "decrypt", "--submission-dir", str(submission_dir), *args],
        )


def _decrypted_with(worker: MagicMock) -> bytes:
    return worker.decrypt.call_args.kwargs["recipient_private_key"].private_bytes_raw()


def _latest_state(db: SubmissionDb) -> SubmissionStateLog:
    submission = db.get_submission(SUBMISSION_ID)
    assert submission is not None
    latest_state = submission.get_latest_state()
    assert latest_state is not None
    return latest_state


@pytest.mark.parametrize(("option", "archive"), [("consented", "consented"), ("non-consented", "non_consented")])
def test_archive_decrypts_with_the_key_of_that_archive(
    offline_config: GrzctlConfig, tmp_path: Path, no_prompt, option: str, archive: str
):
    key_paths = {name: _generate_key(tmp_path, name, PASSPHRASE) for name in ("consented", "non_consented")}
    config_path = _config_path(
        offline_config,
        tmp_path,
        **{
            name: {"private_key_path": str(key_path), "private_key_passphrase": PASSPHRASE}
            for name, key_path in key_paths.items()
        },
    )
    worker = _decrypt_worker()

    result = _decrypt(config_path, tmp_path, worker, "--archive", option, "--no-update-db")

    assert result.exit_code == 0, result.stderr
    assert _decrypted_with(worker) == _raw_key(key_paths[archive], PASSPHRASE)


def test_private_key_path_wins_over_the_key_of_the_archive(offline_config: GrzctlConfig, tmp_path: Path, no_prompt):
    archive_key_path = _generate_key(tmp_path, "consented")
    other_key_path = _generate_key(tmp_path, "other")
    config_path = _config_path(offline_config, tmp_path, consented={"private_key_path": str(archive_key_path)})
    worker = _decrypt_worker()

    result = _decrypt(
        config_path,
        tmp_path,
        worker,
        "--archive",
        "consented",
        "--private-key-path",
        str(other_key_path),
        "--no-update-db",
    )

    assert result.exit_code == 0, result.stderr
    assert _decrypted_with(worker) == _raw_key(other_key_path)


def test_archive_without_a_private_key_records_a_configuration_error(
    migrated_database_config: GrzctlConfig, tmp_path: Path, db: SubmissionDb
):
    config_path = _config_path(migrated_database_config, tmp_path)
    db.add_submission(SUBMISSION_ID)
    worker = _decrypt_worker()

    result = _decrypt(config_path, tmp_path, worker, "--archive", "non-consented")

    assert isinstance(result.exception, ConfigurationError), result.output
    assert "Pass --private-key-path, or set archives.non_consented.private_key_path." in str(result.exception)
    worker.decrypt.assert_not_called()
    latest_state = _latest_state(db)
    assert latest_state.state == SubmissionStateEnum.ERROR
    assert latest_state.failure_reason == FailureReasonEnum.CONFIGURATION_ERROR


def test_private_key_path_decrypts_without_resolving_an_inbox(
    migrated_database_config: GrzctlConfig, tmp_path: Path, db: SubmissionDb, monkeypatch
):
    """Two inboxes and no recorded inbox leave the inbox unknown, which only fails if grzctl resolves it.

    The passphrase of the key file comes from ``C4GH_PASSPHRASE``.
    """
    monkeypatch.setenv("C4GH_PASSPHRASE", PASSPHRASE)
    monkeypatch.setattr("grz_common.utils.crypt.getpass", MagicMock(side_effect=AssertionError("no prompt")))
    key_path = _generate_key(tmp_path, "grz", PASSPHRASE)
    config_path = _config_path(migrated_database_config, tmp_path, inboxes=("inbox-a", "inbox-b"))
    db.add_submission(SUBMISSION_ID)
    worker = _decrypt_worker()

    result = _decrypt(config_path, tmp_path, worker, "--private-key-path", str(key_path))

    assert result.exit_code == 0, result.stderr
    assert _decrypted_with(worker) == _raw_key(key_path, PASSPHRASE)
    assert _latest_state(db).state == SubmissionStateEnum.DECRYPTED


@pytest.mark.parametrize("option", ["--archive", "--private-key-path"])
def test_option_with_inbox_fails_before_the_database_is_touched(
    offline_config: GrzctlConfig, tmp_path: Path, option: str
):
    value = "consented" if option == "--archive" else str(_generate_key(tmp_path, "grz"))
    config_path = _config_path(offline_config, tmp_path)
    worker = _decrypt_worker()

    with patch("grzctl.commands.decrypt.DbContext") as db_context:
        result = _decrypt(config_path, tmp_path, worker, option, value, "--inbox", "inbox")

    assert result.exit_code == 2, result.output
    assert f"{option} and --inbox are mutually exclusive." in result.stderr
    db_context.assert_not_called()
    worker.decrypt.assert_not_called()
