"""Tests for the grzctl ``encrypt`` command."""

import logging
from base64 import b64decode
from pathlib import Path
from unittest.mock import MagicMock, patch

import click.testing
import crypt4gh.keys
import crypt4gh.keys.c4gh
import grz_common.exceptions as grzexc
import grzctl.cli
import pytest
import yaml
from grz_db.models.submission import FailureReasonEnum, SubmissionDb, SubmissionStateEnum, SubmissionStateLog
from grzctl.models.config import GrzctlConfig

from .conftest import PRUEFBERICHT

SUBMITTER_ID = "260914050"
SUBMISSION_ID = "260914050_2025-09-15_c64603a7"
PASSPHRASE = "inbox-key-passphrase"


@pytest.fixture
def key_paths(tmp_path) -> dict[str, Path]:
    """Two crypt4gh private keys, encrypted with ``PASSPHRASE``, by inbox name."""
    paths = {}
    for name in ("inbox-a", "inbox-b"):
        paths[name] = tmp_path / f"{name}.sec"
        crypt4gh.keys.c4gh.generate(paths[name], tmp_path / f"{name}.pub", PASSPHRASE.encode(), comment=None)
    return paths


@pytest.fixture
def db(migrated_database_config: GrzctlConfig) -> SubmissionDb:
    """A SubmissionDb on the migrated database, without an author (adds submissions and records inboxes only)."""
    return SubmissionDb(db_url=migrated_database_config.db.database_url, author=None)


def _raw_key(private_key_path: Path) -> bytes:
    return crypt4gh.keys.get_private_key(private_key_path, lambda: PASSPHRASE)


def _inboxes(key_paths: dict[str, Path]) -> dict[str, dict[str, str]]:
    """Config entries for inboxes that use the keys of *key_paths*, with their passphrase."""
    return {
        name: {"private_key_path": str(path), "private_key_passphrase": PASSPHRASE} for name, path in key_paths.items()
    }


def _write_config(
    tmp_path: Path,
    public_key: str,
    unread_file: str,
    inboxes: dict[str, dict[str, str]],
    database_config: GrzctlConfig | None = None,
) -> Path:
    """Write a config whose submitter :data:`SUBMITTER_ID` has the given inboxes.

    The consented archive has *public_key* inline.
    The ``db`` section is that of *database_config*, else an in-memory database.
    """
    if database_config is None:
        db_section = {"database_url": "sqlite:///:memory:", "author": {"name": "test"}}
    else:
        db_section = database_config.model_dump(mode="json", exclude_none=True, context={"reveal_secrets": True})["db"]
    config = {
        "leistungserbringer": {SUBMITTER_ID: {"inbox_buckets": inboxes}},
        "archives": {
            "consented": {"s3": {"bucket": "consented"}, "public_key": public_key},
            "non_consented": {"s3": {"bucket": "non_consented"}, "public_key_path": unread_file},
        },
        "db": db_section,
        "pruefbericht": PRUEFBERICHT,
        "identifiers": {"grz": "GRZT00000"},
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config))
    return config_path


def _submission_dir(tmp_path: Path) -> Path:
    submission_dir = tmp_path / "submission"
    for sub in ("metadata", "files", "logs", "encrypted_files"):
        (submission_dir / sub).mkdir(parents=True)
    return submission_dir


def _invoke_encrypt(
    config_path: Path, submission_dir: Path, mock_worker_cls: MagicMock, *args: str
) -> click.testing.Result:
    """Run ``grzctl encrypt`` with the extra *args* against a mocked ``Worker`` for a consented submission."""
    mock_worker = mock_worker_cls.return_value
    mock_submission = mock_worker.parse_submission.return_value
    mock_submission.metadata.content.submission_id = SUBMISSION_ID
    mock_submission.metadata.content.submission.submitter_id = SUBMITTER_ID
    mock_submission.metadata.content.consents_to_research.return_value = True

    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    return runner.invoke(cli, ["--config", str(config_path), "encrypt", "--submission-dir", str(submission_dir), *args])


def _signing_key(mock_worker_cls: MagicMock):
    return mock_worker_cls.return_value.encrypt.call_args.kwargs["submitter_private_key"]


def _latest_state(db: SubmissionDb) -> SubmissionStateLog | None:
    submission = db.get_submission(SUBMISSION_ID)
    assert submission is not None
    return submission.get_latest_state()


def test_encrypt_uses_an_inline_archive_public_key(tmp_path, key_paths, crypt4gh_public_key, unread_file, no_prompt):
    """The consented archive's public key may be given inline instead of as a file path."""
    config_path = _write_config(tmp_path, crypt4gh_public_key, unread_file, _inboxes({"main": key_paths["inbox-a"]}))

    with patch("grzctl.commands.encrypt.Worker") as mock_worker_cls:
        result = _invoke_encrypt(config_path, _submission_dir(tmp_path), mock_worker_cls, "--no-update-db")

    assert result.exit_code == 0, result.output
    mock_worker_cls.return_value.encrypt.assert_called_once()
    recipient_public_key = mock_worker_cls.return_value.encrypt.call_args.kwargs["recipient_public_key"]
    assert recipient_public_key.public_bytes_raw() == b64decode(crypt4gh_public_key.splitlines()[1])


def test_encrypt_signs_with_the_key_of_the_only_inbox(tmp_path, key_paths, crypt4gh_public_key, unread_file, no_prompt):
    """The key of the submitter's only inbox reaches ``Worker.encrypt``, decrypted with its configured passphrase."""
    config_path = _write_config(tmp_path, crypt4gh_public_key, unread_file, _inboxes({"main": key_paths["inbox-a"]}))

    with patch("grzctl.commands.encrypt.Worker") as mock_worker_cls:
        result = _invoke_encrypt(config_path, _submission_dir(tmp_path), mock_worker_cls, "--no-update-db")

    assert result.exit_code == 0, result.output
    assert _signing_key(mock_worker_cls).private_bytes_raw() == _raw_key(key_paths["inbox-a"])


def test_encrypt_signs_with_the_key_of_the_recorded_inbox(
    tmp_path, key_paths, migrated_database_config, db, crypt4gh_public_key, unread_file, no_prompt
):
    """With several inboxes, the inbox recorded in the database names the signing key."""
    config_path = _write_config(
        tmp_path, crypt4gh_public_key, unread_file, _inboxes(key_paths), database_config=migrated_database_config
    )
    db.add_submission(SUBMISSION_ID)
    db.set_submission_inbox(SUBMISSION_ID, "inbox-b")

    with patch("grzctl.commands.encrypt.Worker") as mock_worker_cls:
        result = _invoke_encrypt(config_path, _submission_dir(tmp_path), mock_worker_cls, "--update-db")

    assert result.exit_code == 0, result.output
    assert _signing_key(mock_worker_cls).private_bytes_raw() == _raw_key(key_paths["inbox-b"])
    latest_state = _latest_state(db)
    assert latest_state is not None
    assert latest_state.state == SubmissionStateEnum.ENCRYPTED


@pytest.mark.parametrize(
    ("update_db_flag", "recorded_inbox"),
    [("--update-db", None), ("--no-update-db", "inbox-b")],
    ids=["nothing-recorded", "database-not-read"],
)
def test_encrypt_signs_with_a_random_key_if_no_inbox_resolves(
    tmp_path,
    key_paths,
    migrated_database_config,
    db,
    crypt4gh_public_key,
    unread_file,
    caplog,
    update_db_flag: str,
    recorded_inbox: str | None,
):
    """Several inboxes and no recorded inbox leave the inbox unknown, so grz-common signs with a random key.

    grzctl reads the recorded inbox only with ``--update-db``.
    """
    config_path = _write_config(
        tmp_path, crypt4gh_public_key, unread_file, _inboxes(key_paths), database_config=migrated_database_config
    )
    db.add_submission(SUBMISSION_ID)
    if recorded_inbox is not None:
        db.set_submission_inbox(SUBMISSION_ID, recorded_inbox)

    with (
        patch("grzctl.commands.encrypt.Worker") as mock_worker_cls,
        caplog.at_level(logging.WARNING, logger="grzctl.commands.encrypt"),
    ):
        result = _invoke_encrypt(config_path, _submission_dir(tmp_path), mock_worker_cls, update_db_flag)

    assert result.exit_code == 0, result.output
    assert _signing_key(mock_worker_cls) is None
    warnings = [r.getMessage() for r in caplog.records if r.name == "grzctl.commands.encrypt"]
    assert len(warnings) == 1, warnings
    assert f"No inbox resolves for submission {SUBMISSION_ID}" in warnings[0]
    assert "signed with a random key" in warnings[0]
    assert "grzctl db backfill" in warnings[0]


def test_encrypt_fails_for_a_recorded_inbox_missing_from_the_config(
    tmp_path, key_paths, migrated_database_config, db, crypt4gh_public_key, unread_file
):
    """A recorded inbox that the config no longer names is a configuration error, as in ``grzctl decrypt``."""
    config_path = _write_config(
        tmp_path, crypt4gh_public_key, unread_file, _inboxes(key_paths), database_config=migrated_database_config
    )
    db.add_submission(SUBMISSION_ID)
    db.set_submission_inbox(SUBMISSION_ID, "removed-inbox")

    with patch("grzctl.commands.encrypt.Worker") as mock_worker_cls:
        result = _invoke_encrypt(config_path, _submission_dir(tmp_path), mock_worker_cls, "--update-db")

    assert isinstance(result.exception, grzexc.ConfigurationError), result.output
    assert "Inbox 'removed-inbox' not configured" in str(result.exception)
    mock_worker_cls.return_value.encrypt.assert_not_called()
    latest_state = _latest_state(db)
    assert latest_state is not None
    assert latest_state.state == SubmissionStateEnum.ERROR
    assert latest_state.failure_reason == FailureReasonEnum.CONFIGURATION_ERROR


def test_encrypt_fails_if_the_inbox_key_cannot_be_loaded(
    tmp_path, migrated_database_config, db, crypt4gh_public_key, unread_file
):
    """An inbox key that cannot be loaded is a configuration error, which the GRZ has to fix."""
    not_a_key_path = tmp_path / "not_a_key.sec"
    not_a_key_path.write_text("not a key")
    config_path = _write_config(
        tmp_path,
        crypt4gh_public_key,
        unread_file,
        {"main": {"private_key_path": str(not_a_key_path)}},
        database_config=migrated_database_config,
    )
    db.add_submission(SUBMISSION_ID)

    with patch("grzctl.commands.encrypt.Worker") as mock_worker_cls:
        result = _invoke_encrypt(config_path, _submission_dir(tmp_path), mock_worker_cls, "--update-db")

    assert isinstance(result.exception, grzexc.ConfigurationError), result.output
    assert f"Secret key {not_a_key_path} cannot be read" in str(result.exception)
    mock_worker_cls.return_value.encrypt.assert_not_called()
    latest_state = _latest_state(db)
    assert latest_state is not None
    assert latest_state.failure_reason == FailureReasonEnum.CONFIGURATION_ERROR
