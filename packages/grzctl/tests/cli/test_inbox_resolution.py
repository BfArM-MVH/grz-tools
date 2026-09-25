"""Tests for resolving the inbox of a submission.

Resolution order: explicit ``--inbox``, the inbox recorded in the database, the
submitter's only inbox, and (when asked) the sole inbox that S3 still holds the
submission in.
"""

import boto3
import click
import grzctl.commands.inbox_resolution as inbox_resolution
import pytest
from grzctl.models.config import GrzctlConfig
from moto import mock_aws

SUBMITTER_ID = "260914050"
SUBMISSION_ID = "260914050_2024-01-01_abcdef01"
REGION = "us-east-1"


@pytest.fixture(autouse=True)
def _clear_inbox_listing_cache() -> None:
    """The listing cache lives for the process; a test must not inherit previous listings."""
    inbox_resolution._inbox_listing_cache.clear()


def _config(tmp_path, unread_file: str, inbox_names: list[str]) -> GrzctlConfig:
    archives = {
        name: {"s3": {"bucket": name}, "public_key_path": unread_file} for name in ("consented", "non_consented")
    }
    return GrzctlConfig.from_configuration(
        {
            "leistungserbringer": {
                SUBMITTER_ID: {
                    "inbox_buckets": {inbox_name: {"private_key_path": unread_file} for inbox_name in inbox_names}
                }
            },
            "archives": {**archives, "signing_key_path": unread_file},
            "db": {"database_url": f"sqlite:///{tmp_path / 'unused.sqlite'}", "author": {"name": "test"}},
            "pruefbericht": {},
            "identifiers": {"grz": "GRZK00007"},
        }
    )


def _db(recorded_inbox: str | None = None, *, error: Exception | None = None):
    """A mocked submission database answering ``get_submission`` with *recorded_inbox*."""
    submission = type("Submission", (), {"inbox": recorded_inbox})() if recorded_inbox is not None else None

    class Db:
        def get_submission(self, submission_id: str):
            if error is not None:
                raise error
            return submission

    return Db()


def _one_inbox_listing(submission_id: str, inbox_name: str):
    """Mock ``query_submissions`` so only *inbox_name* still holds *submission_id*."""

    def _listing(s3_options, show_cleaned: bool):
        summary = type("Summary", (), {"submission_id": submission_id})()
        return [summary] if s3_options.bucket == inbox_name else []

    return _listing


def test_an_explicit_inbox_wins_over_every_source(tmp_path, unread_file):
    config = _config(tmp_path, unread_file, ["inbox-a", "inbox-b"])
    db = _db(recorded_inbox="inbox-b")

    resolved = inbox_resolution.resolve_inbox(
        config, submitter_id=SUBMITTER_ID, submission_id=SUBMISSION_ID, inbox_name="inbox-a", db_service=db, scan=True
    )

    assert resolved == "inbox-a"


def test_the_database_record_wins_over_the_only_inbox(tmp_path, unread_file):
    """A recorded inbox is authoritative, even when the submitter now has just one inbox."""
    config = _config(tmp_path, unread_file, ["inbox-b"])

    resolved = inbox_resolution.resolve_inbox(
        config, submitter_id=SUBMITTER_ID, submission_id=SUBMISSION_ID, db_service=_db(recorded_inbox="inbox-b")
    )

    assert resolved == "inbox-b"


def test_the_database_record_wins_over_the_s3_scan(tmp_path, unread_file):
    config = _config(tmp_path, unread_file, ["inbox-a", "inbox-b"])
    listing = _one_inbox_listing(SUBMISSION_ID, "inbox-a")

    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.setattr(inbox_resolution, "query_submissions", listing)
        resolved = inbox_resolution.resolve_inbox(
            config,
            submitter_id=SUBMITTER_ID,
            submission_id=SUBMISSION_ID,
            db_service=_db(recorded_inbox="inbox-b"),
            scan=True,
        )

    assert resolved == "inbox-b"


def test_the_only_inbox_falls_back_when_nothing_is_recorded(tmp_path, unread_file):
    config = _config(tmp_path, unread_file, ["inbox"])

    resolved = inbox_resolution.resolve_inbox(config, submitter_id=SUBMITTER_ID, submission_id=SUBMISSION_ID)

    assert resolved == "inbox"


def test_scan_finds_the_sole_inbox_holding_the_submission(tmp_path, unread_file):
    config = _config(tmp_path, unread_file, ["inbox-a", "inbox-b"])
    listing = _one_inbox_listing(SUBMISSION_ID, "inbox-a")

    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.setattr(inbox_resolution, "query_submissions", listing)
        resolved = inbox_resolution.resolve_inbox(
            config, submitter_id=SUBMITTER_ID, submission_id=SUBMISSION_ID, scan=True
        )

    assert resolved == "inbox-a"


def test_scan_is_ambiguous_when_several_inboxes_hold_the_submission(tmp_path, unread_file):
    config = _config(tmp_path, unread_file, ["inbox-a", "inbox-b"])
    summary = type("Summary", (), {"submission_id": SUBMISSION_ID})()

    def listing(s3_options, show_cleaned: bool):
        return [summary]

    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.setattr(inbox_resolution, "query_submissions", listing)
        resolved = inbox_resolution.resolve_inbox(
            config, submitter_id=SUBMITTER_ID, submission_id=SUBMISSION_ID, scan=True
        )

    assert resolved is None


def test_scan_finds_no_inbox_without_a_marker(tmp_path, unread_file):
    config = _config(tmp_path, unread_file, ["inbox-a", "inbox-b"])
    listing = _one_inbox_listing("other_2024-01-01_abcdef02", "inbox-a")

    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.setattr(inbox_resolution, "query_submissions", listing)
        resolved = inbox_resolution.resolve_inbox(
            config, submitter_id=SUBMITTER_ID, submission_id=SUBMISSION_ID, scan=True
        )

    assert resolved is None


def test_a_database_that_errors_falls_through(tmp_path, unread_file):
    """A database that does not answer (not configured, out of sync) must not block resolution."""
    config = _config(tmp_path, unread_file, ["inbox"])

    resolved = inbox_resolution.resolve_inbox(
        config, submitter_id=SUBMITTER_ID, submission_id=SUBMISSION_ID, db_service=_db(error=RuntimeError("boom"))
    )

    assert resolved == "inbox"


@pytest.mark.parametrize("scan", [False, True])
def test_several_inboxes_without_any_source_resolve_to_none(tmp_path, unread_file, scan: bool):
    config = _config(tmp_path, unread_file, ["inbox-a", "inbox-b"])

    with pytest.MonkeyPatch.context() as monkeypatch:
        if scan:
            monkeypatch.setattr(inbox_resolution, "query_submissions", lambda s3_options, show_cleaned: [])
        resolved = inbox_resolution.resolve_inbox(
            config, submitter_id=SUBMITTER_ID, submission_id=SUBMISSION_ID, scan=scan
        )

    assert resolved is None


@mock_aws
def test_scan_inbox_uses_the_metadata_marker_still_kept_in_the_inbox(tmp_path, unread_file):
    """``grzctl clean`` leaves a metadata.json marker behind, so cleaning must not hide the inbox."""
    inboxes = {"inbox-a": "bucket-a", "inbox-b": "bucket-b"}
    config = _config(tmp_path, unread_file, list(inboxes))
    # re-point the inboxes at the moto buckets
    data = config.model_dump(mode="json", exclude_none=True, context={"reveal_secrets": True})
    for inbox_name, bucket in inboxes.items():
        data["leistungserbringer"][SUBMITTER_ID]["inbox_buckets"][inbox_name].update(
            {"bucket": bucket, "region_name": REGION}
        )
    config = GrzctlConfig.from_configuration(data)

    for bucket in inboxes.values():
        boto3.client("s3", region_name=REGION).create_bucket(Bucket=bucket)
    s3_client = boto3.client("s3", region_name=REGION)
    s3_client.put_object(Bucket="bucket-a", Key=f"{SUBMISSION_ID}/metadata/metadata.json", Body=b"")

    resolved = inbox_resolution.scan_inbox(config, SUBMITTER_ID, SUBMISSION_ID)

    assert resolved == "inbox-a"


def test_require_inbox_returns_the_resolved_inbox(tmp_path, unread_file):
    config = _config(tmp_path, unread_file, ["inbox"])

    assert inbox_resolution.require_inbox(config, submitter_id=SUBMITTER_ID, submission_id=SUBMISSION_ID) == "inbox"


def test_require_inbox_aborts_when_nothing_resolves(tmp_path, unread_file):
    config = _config(tmp_path, unread_file, ["inbox-a", "inbox-b"])

    with pytest.raises(click.ClickException, match="several inboxes \\(inbox-a, inbox-b\\)"):
        inbox_resolution.require_inbox(config, submitter_id=SUBMITTER_ID, submission_id=SUBMISSION_ID)


def test_require_inbox_aborts_without_a_submitter_entry(tmp_path, unread_file):
    config = _config(tmp_path, unread_file, ["inbox"])

    with pytest.raises(click.ClickException, match="no inbox in the configuration"):
        inbox_resolution.require_inbox(config, submitter_id="999999999", submission_id=SUBMISSION_ID)


def test_require_inbox_repeats_the_hint(tmp_path, unread_file):
    config = _config(tmp_path, unread_file, ["inbox-a", "inbox-b"])

    with pytest.raises(click.ClickException, match="Record the inbox by downloading it once"):
        inbox_resolution.require_inbox(
            config,
            submitter_id=SUBMITTER_ID,
            submission_id=SUBMISSION_ID,
            hint="Record the inbox by downloading it once.",
        )


def test_require_inbox_raises_a_usage_error_when_asked(tmp_path, unread_file):
    config = _config(tmp_path, unread_file, ["inbox-a", "inbox-b"])

    with pytest.raises(click.UsageError):
        inbox_resolution.require_inbox(config, submitter_id=SUBMITTER_ID, exc_type=click.UsageError)


def test_db_inbox_returns_the_recorded_inbox_or_none(tmp_path, unread_file):
    recorded = inbox_resolution.db_inbox(_db(recorded_inbox="inbox-b"), SUBMISSION_ID)
    assert recorded == "inbox-b"

    errored = inbox_resolution.db_inbox(_db(error=RuntimeError("boom")), SUBMISSION_ID)
    assert errored is None

    unrecorded = inbox_resolution.db_inbox(_db(recorded_inbox=""), SUBMISSION_ID)
    assert unrecorded is None
