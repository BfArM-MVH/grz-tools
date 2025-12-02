"""
Tests for the 'grzctl db sync-from-inbox' subcommand.
"""

from datetime import UTC, datetime
from unittest.mock import patch

import click.testing
import grzctl.cli
import yaml
from grz_common.workers.download import InboxSubmissionState, InboxSubmissionSummary
from grz_db.models.author import Author
from grz_db.models.submission import SubmissionDb, SubmissionStateEnum
from grzctl.models.config import DbConfig


def test_sync_from_inbox(blank_database_config_path, tmp_path):
    """
    Test that sync-from-inbox correctly updates the database based on inbox state.
    """
    with open(blank_database_config_path, encoding="utf-8") as f:
        config_data = yaml.safe_load(f)

    config_data["s3"] = {
        "bucket": "test-bucket",
        "endpoint_url": "http://localhost:9000",
        "access_key_id": "minioadmin",
        "secret_access_key": "minioadmin",
    }

    with open(blank_database_config_path, "w", encoding="utf-8") as f:
        yaml.dump(config_data, f)

    old = datetime.fromisoformat("1970-01-01")
    now = datetime.now(UTC)

    # complete in inbox + new in db → uploaded
    sub_new_complete = InboxSubmissionSummary(
        submission_id="123456789_1970-01-01_00000000",
        state=InboxSubmissionState.COMPLETE,
        oldest_upload=old,
        newest_upload=now,
        total_size_bytes=100,
    )

    # incomplete in inbox + new in db → uploading
    sub_new_incomplete = InboxSubmissionSummary(
        submission_id="123456789_1970-01-01_11111111",
        state=InboxSubmissionState.INCOMPLETE,
        oldest_upload=old,
        newest_upload=now,
        total_size_bytes=100,
    )

    # complete in inbox + uploading in db → uploaded
    sub_update_needed = InboxSubmissionSummary(
        submission_id="123456789_1970-01-01_22222222",
        state=InboxSubmissionState.COMPLETE,
        oldest_upload=old,
        newest_upload=now,
        total_size_bytes=100,
    )

    # uploaded in db + incomplete in inbox (mismatch!) → no change
    sub_no_change = InboxSubmissionSummary(
        submission_id="123456789_1970-01-01_33333333",
        state=InboxSubmissionState.INCOMPLETE,
        oldest_upload=old,
        newest_upload=now,
        total_size_bytes=100,
    )

    mock_s3_submissions = [
        sub_new_complete,
        sub_new_incomplete,
        sub_update_needed,
        sub_no_change,
    ]

    db_config = DbConfig.from_path(blank_database_config_path)

    with open(config_data["db"]["author"]["private_key_path"], "rb") as f:
        pk_bytes = f.read()
    author = Author(name="alice", private_key_bytes=pk_bytes, private_key_passphrase="")

    db = SubmissionDb(db_url=db_config.db.database_url, author=author)

    db.add_submission("123456789_1970-01-01_22222222")
    db.update_submission_state("123456789_1970-01-01_22222222", SubmissionStateEnum.UPLOADING)

    db.add_submission("123456789_1970-01-01_33333333")
    db.update_submission_state("123456789_1970-01-01_33333333", SubmissionStateEnum.UPLOADED)

    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()

    with patch("grzctl.commands.db.cli.query_submissions", return_value=mock_s3_submissions) as mock_query:
        result = runner.invoke(cli, ["db", "--config-file", str(blank_database_config_path), "sync-from-inbox"])

        assert result.exit_code == 0, result.stderr
        mock_query.assert_called_once()

    s1 = db.get_submission("123456789_1970-01-01_00000000")
    assert s1 is not None
    assert s1.get_latest_state().state == SubmissionStateEnum.UPLOADED

    s2 = db.get_submission("123456789_1970-01-01_11111111")
    assert s2 is not None
    assert s2.get_latest_state().state == SubmissionStateEnum.UPLOADING

    s3 = db.get_submission("123456789_1970-01-01_22222222")
    assert s3.get_latest_state().state == SubmissionStateEnum.UPLOADED

    s4 = db.get_submission("123456789_1970-01-01_33333333")
    assert s4.get_latest_state().state == SubmissionStateEnum.UPLOADED
