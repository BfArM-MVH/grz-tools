"""Tests for recording and using the inbox of a submission.

``grzctl download``, ``grzctl db sync-from-inbox`` and ``grzctl db submission populate --inbox``
record the inbox a submission came from; ``grzctl decrypt`` then picks the key of that inbox.
"""

from datetime import UTC, datetime
from pathlib import Path
from unittest.mock import MagicMock, patch

import click.testing
import crypt4gh.keys
import crypt4gh.keys.c4gh
import grzctl.cli
import pytest
import yaml
from grz_common.exceptions import ConfigurationError
from grz_common.workers.download import InboxSubmissionState, InboxSubmissionSummary
from grz_db.models.submission import FailureReasonEnum, SubmissionDb, SubmissionStateEnum, SubmissionStateLog
from grz_pydantic_models.submission.metadata import GrzSubmissionMetadata
from grzctl.models.config import GrzctlConfig

SUBMITTER_ID = "260914050"
SUBMISSION_ID = "260914050_2025-09-15_c64603a7"
INBOX = "inbox"


@pytest.fixture
def db(migrated_database_config: GrzctlConfig) -> SubmissionDb:
    """A SubmissionDb on the migrated database, without an author (reads and writes inboxes only)."""
    return SubmissionDb(db_url=migrated_database_config.db.database_url, author=None)


def _config_with_inbox(
    migrated_database_config: GrzctlConfig,
    tmp_path: Path,
    inboxes=("inbox",),
    key_paths: dict[str, Path] | None = None,
) -> Path:
    """A config file whose submitter :data:`SUBMITTER_ID` has the given inboxes.

    Each inbox uses the private key that *key_paths* names for it, else the key of the DB author.
    """
    data = migrated_database_config.model_dump(mode="json", exclude_none=True, context={"reveal_secrets": True})
    author_key_path = data["db"]["author"]["private_key_path"]
    key_paths = key_paths or {}
    data["leistungserbringer"] = {
        SUBMITTER_ID: {
            "inbox_buckets": {name: {"private_key_path": str(key_paths.get(name, author_key_path))} for name in inboxes}
        }
    }
    config_path = tmp_path / "config.inbox.yaml"
    config_path.write_text(yaml.safe_dump(data))
    return config_path


def _invoke(*args: str) -> click.testing.Result:
    runner = click.testing.CliRunner()
    return runner.invoke(grzctl.cli.build_cli(), list(args))


def test_sync_from_inbox_records_inbox(migrated_database_config: GrzctlConfig, tmp_path: Path, db: SubmissionDb):
    """Every submission the scan finds is recorded with the inbox it was found in."""
    old = datetime.fromisoformat("1970-01-01T00:00:00+00:00")
    now = datetime.now(UTC)
    found = [
        InboxSubmissionSummary(
            submission_id=f"{SUBMITTER_ID}_2025-01-0{i}_0000000{i}",
            state=InboxSubmissionState.COMPLETE,
            oldest_upload=old,
            newest_upload=now,
            total_size_bytes=100,
        )
        for i in (1, 2)
    ]
    db.add_submission(found[1].submission_id)
    config_path = _config_with_inbox(migrated_database_config, tmp_path)

    with patch("grzctl.commands.db.cli.query_submissions", return_value=found):
        result = _invoke(
            "--config",
            str(config_path),
            "db",
            "sync-from-inbox",
            "--submitter-id",
            SUBMITTER_ID,
            "--inbox",
            INBOX,
        )

    assert result.exit_code == 0, result.stderr
    for summary in found:
        recorded = db.get_submission(summary.submission_id)
        assert recorded is not None
        assert recorded.inbox == INBOX


def test_sync_from_inbox_defaults_to_the_only_inbox(
    migrated_database_config: GrzctlConfig, tmp_path: Path, db: SubmissionDb
):
    """Without --inbox, a submitter with one inbox is unambiguous and records it."""
    summary = InboxSubmissionSummary(
        submission_id=f"{SUBMITTER_ID}_2025-01-01_00000001",
        state=InboxSubmissionState.COMPLETE,
        oldest_upload=datetime.fromisoformat("1970-01-01T00:00:00+00:00"),
        newest_upload=datetime.now(UTC),
        total_size_bytes=100,
    )
    config_path = _config_with_inbox(migrated_database_config, tmp_path)

    with patch("grzctl.commands.db.cli.query_submissions", return_value=[summary]):
        result = _invoke("--config", str(config_path), "db", "sync-from-inbox", "--submitter-id", SUBMITTER_ID)

    assert result.exit_code == 0, result.stderr
    recorded = db.get_submission(summary.submission_id)
    assert recorded is not None
    assert recorded.inbox == INBOX


def test_populate_with_inbox_records_inbox(
    migrated_database_config: GrzctlConfig, tmp_path: Path, test_metadata_path: Path, db: SubmissionDb
):
    """``db submission populate --inbox`` records the inbox, even when the diff has nothing to write."""
    metadata = GrzSubmissionMetadata.model_validate_json(test_metadata_path.read_text())
    db.add_submission(metadata.submission_id)
    config_path = _config_with_inbox(migrated_database_config, tmp_path)

    result = _invoke(
        "--config",
        str(config_path),
        "db",
        "submission",
        "populate",
        metadata.submission_id,
        str(test_metadata_path),
        "--submission-date",
        "2025-09-15",
        "--inbox",
        INBOX,
        "--no-confirm",
    )

    assert result.exit_code == 0, result.stderr
    recorded = db.get_submission(metadata.submission_id)
    assert recorded is not None
    assert recorded.inbox == INBOX


def test_populate_records_the_inbox_even_when_nothing_is_written(
    migrated_database_config: GrzctlConfig, tmp_path: Path, test_metadata_path: Path, db: SubmissionDb
):
    """A populate whose diff is empty (already up to date) still records an explicit --inbox."""
    metadata = GrzSubmissionMetadata.model_validate_json(test_metadata_path.read_text())
    db.add_submission(metadata.submission_id)
    db.populate(metadata.submission_id, metadata, metadata.submission.submission_date)
    config_path = _config_with_inbox(migrated_database_config, tmp_path)

    result = _invoke(
        "--config",
        str(config_path),
        "db",
        "submission",
        "populate",
        metadata.submission_id,
        str(test_metadata_path),
        "--submission-date",
        "2025-09-15",
        "--inbox",
        INBOX,
        "--no-confirm",
    )

    assert result.exit_code == 0, result.stderr
    assert db.get_submission(metadata.submission_id).inbox == INBOX


def test_populate_without_inbox_records_nothing(
    migrated_database_config: GrzctlConfig, tmp_path: Path, test_metadata_path: Path, db: SubmissionDb
):
    """Without an explicit --inbox, populate does not guess the inbox."""
    metadata = GrzSubmissionMetadata.model_validate_json(test_metadata_path.read_text())
    db.add_submission(metadata.submission_id)
    config_path = _config_with_inbox(migrated_database_config, tmp_path)

    result = _invoke(
        "--config",
        str(config_path),
        "db",
        "submission",
        "populate",
        metadata.submission_id,
        str(test_metadata_path),
        "--submission-date",
        "2025-09-15",
        "--no-confirm",
    )

    assert result.exit_code == 0, result.stderr
    recorded = db.get_submission(metadata.submission_id)
    assert recorded is not None
    assert recorded.inbox is None


def test_download_records_inbox(migrated_database_config: GrzctlConfig, tmp_path: Path, db: SubmissionDb):
    """A successful download records the inbox it came from, independent of populate."""
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    db.add_submission(SUBMISSION_ID)
    config_path = _config_with_inbox(migrated_database_config, tmp_path)

    context = MagicMock()
    context.__enter__.return_value.db = db
    with (
        patch("grzctl.commands.download.DbContext", return_value=context) as db_context,
        patch("grzctl.commands.download.Worker") as worker_cls,
    ):
        worker = worker_cls.return_value
        result = _invoke(
            "--config",
            str(config_path),
            "download",
            "--submission-id",
            SUBMISSION_ID,
            "--output-dir",
            str(output_dir),
            "--inbox",
            INBOX,
            "--no-populate",
        )

    assert result.exit_code == 0, result.stderr
    db_context.assert_called_once()
    recorded = db.get_submission(SUBMISSION_ID)
    assert recorded is not None
    assert recorded.inbox == INBOX
    worker.download.assert_called_once()


def test_download_without_update_db_touches_no_database(migrated_database_config: GrzctlConfig, tmp_path: Path):
    """With --no-update-db, download neither reads the recorded inbox nor records one, even with --populate."""
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    config_path = _config_with_inbox(migrated_database_config, tmp_path)

    no_database = AssertionError("no database may be opened")
    with (
        patch("grzctl.commands.download.get_submission_db_instance", side_effect=no_database),
        patch("grzctl.dbcontext.get_submission_db_instance", side_effect=no_database),
        patch("grzctl.commands.download.Worker") as worker_cls,
    ):
        result = _invoke(
            "--config",
            str(config_path),
            "download",
            "--submission-id",
            SUBMISSION_ID,
            "--output-dir",
            str(output_dir),
            "--no-update-db",
            "--populate",
        )

    assert result.exit_code == 0, result.stderr
    worker_cls.return_value.download.assert_called_once()


def _decrypt_worker() -> MagicMock:
    encrypted = MagicMock()
    encrypted.submission_id = SUBMISSION_ID
    encrypted.metadata.content.submission.submitter_id = SUBMITTER_ID
    worker = MagicMock()
    worker.parse_encrypted_submission.return_value = encrypted
    return worker


def _submission_dir(tmp_path: Path, db: SubmissionDb, inbox: str | None) -> Path:
    """A submission directory plus a database that records *inbox* (or nothing) for it."""
    db.add_submission(SUBMISSION_ID)
    if inbox is not None:
        db.set_submission_inbox(SUBMISSION_ID, inbox)
    submission_dir = tmp_path / "submission"
    submission_dir.mkdir()
    return submission_dir


def _decrypt(config_path: Path, submission_dir: Path, worker: MagicMock) -> click.testing.Result:
    """Run ``grzctl decrypt`` without --inbox, with the mocked *worker* and the real ``DbContext``."""
    with patch("grzctl.commands.decrypt.Worker", return_value=worker):
        return _invoke("--config", str(config_path), "decrypt", "--submission-dir", str(submission_dir))


def _latest_state(db: SubmissionDb) -> SubmissionStateLog:
    submission = db.get_submission(SUBMISSION_ID)
    assert submission is not None
    latest_state = submission.get_latest_state()
    assert latest_state is not None
    return latest_state


def test_decrypt_uses_the_key_of_the_recorded_inbox(
    migrated_database_config: GrzctlConfig, tmp_path: Path, db: SubmissionDb
):
    """A recorded inbox names its own key, so two inboxes with different keys still decrypt."""
    key_paths = {}
    for name in ("inbox-a", "inbox-b"):
        key_paths[name] = tmp_path / f"{name}.sec"
        crypt4gh.keys.c4gh.generate(key_paths[name], tmp_path / f"{name}.pub", None, comment=None)
    config_path = _config_with_inbox(
        migrated_database_config, tmp_path, inboxes=("inbox-a", "inbox-b"), key_paths=key_paths
    )
    submission_dir = _submission_dir(tmp_path, db, inbox="inbox-b")
    worker = _decrypt_worker()

    result = _decrypt(config_path, submission_dir, worker)

    assert result.exit_code == 0, result.stderr
    private_key = worker.decrypt.call_args.kwargs["recipient_private_key"]
    assert private_key.private_bytes_raw() == crypt4gh.keys.get_private_key(key_paths["inbox-b"], None)
    assert _latest_state(db).state == SubmissionStateEnum.DECRYPTED


def test_decrypt_without_a_recorded_inbox_fails_for_several_inboxes(
    migrated_database_config: GrzctlConfig, tmp_path: Path, db: SubmissionDb
):
    """Without --inbox and a recorded inbox, several inboxes leave the key unknown."""
    config_path = _config_with_inbox(migrated_database_config, tmp_path, inboxes=("inbox-a", "inbox-b"))
    submission_dir = _submission_dir(tmp_path, db, inbox=None)
    worker = _decrypt_worker()

    result = _decrypt(config_path, submission_dir, worker)

    assert isinstance(result.exception, ConfigurationError), result.output
    worker.decrypt.assert_not_called()
    assert "Pass --inbox, or record the inbox with 'grzctl db backfill'." in str(result.exception)
    latest_state = _latest_state(db)
    assert latest_state.state == SubmissionStateEnum.ERROR
    assert latest_state.failure_reason == FailureReasonEnum.CONFIGURATION_ERROR


def test_decrypt_fails_for_a_recorded_inbox_missing_from_the_config(
    migrated_database_config: GrzctlConfig, tmp_path: Path, db: SubmissionDb
):
    """A recorded inbox that the config no longer names is a configuration error."""
    config_path = _config_with_inbox(migrated_database_config, tmp_path)
    submission_dir = _submission_dir(tmp_path, db, inbox="removed-inbox")
    worker = _decrypt_worker()

    result = _decrypt(config_path, submission_dir, worker)

    assert isinstance(result.exception, ConfigurationError), result.output
    worker.decrypt.assert_not_called()
    assert "Inbox 'removed-inbox' not configured" in str(result.exception)
    latest_state = _latest_state(db)
    assert latest_state.state == SubmissionStateEnum.ERROR
    assert latest_state.failure_reason == FailureReasonEnum.CONFIGURATION_ERROR
