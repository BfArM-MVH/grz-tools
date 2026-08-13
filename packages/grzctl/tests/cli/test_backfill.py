# Backfill unit tests — see .vbw-planning/phases/02-implement-grzctl-db-backfill-command/02-02-PLAN.md
"""Unit tests for `_backfill_submission`.

The tests target `_backfill_submission` directly with a real SubmissionDb on every supported
backend and a moto-mocked S3. Live in `grzctl/tests/` (not `grz-db/tests/`) because `grzctl`
is not a dev/test dependency of `grz-db`, but `moto[s3]`, `pytest-postgresql`, and
`grz-pydantic-models-testing` are all in `grzctl`'s [test] dependency group.
"""

import datetime
import importlib.resources
import json
from collections.abc import Iterator
from typing import Any

import boto3
import pytest
from grz_db.models.submission import Submission, SubmissionDb
from grz_pydantic_models.submission.metadata import GrzSubmissionMetadata
from grz_pydantic_models_testing.example_metadata import grzctl as grzctl_metadata
from grzctl.commands.db.cli import _backfill_submission, _BackfillResult
from moto import mock_aws

BUCKET = "test-backfill-bucket"
REGION = "us-east-1"
IGNORE_FIELDS = {"submission_uploaded_date", "tan_g", "pseudonym"}
DIFFERENT_TAN_G = "b" * 64
DIFFERENT_PSEUDONYM = "different-pseudonym"
DIFFERENT_DATE = datetime.date(1999, 1, 1)


@pytest.fixture(scope="session")
def metadata() -> GrzSubmissionMetadata:
    """Load the wes_tumor_germline v1.2.1 example shipped with grz-pydantic-models-testing."""
    path = importlib.resources.files(grzctl_metadata).joinpath("metadata.json")
    with path.open() as fh:
        return GrzSubmissionMetadata(**json.load(fh))


@pytest.fixture(scope="session")
def submission_id(metadata: GrzSubmissionMetadata) -> str:
    return metadata.submission_id


@pytest.fixture
def s3_client_mock() -> Iterator[Any]:
    """A moto-backed S3 client with a pre-created test bucket."""
    with mock_aws():
        client = boto3.client("s3", region_name=REGION)
        client.create_bucket(Bucket=BUCKET)
        yield client


def _put_metadata(s3_client: Any, submission_id: str, metadata: GrzSubmissionMetadata) -> None:
    s3_client.put_object(
        Bucket=BUCKET,
        Key=f"{submission_id}/metadata/metadata.json",
        Body=metadata.model_dump_json(by_alias=True).encode("utf-8"),
    )


def _populate_full_row(db: SubmissionDb, submission_id: str, metadata: GrzSubmissionMetadata) -> Submission:
    """Persist a fully-populated row by running the same diff/commit path the production code uses."""
    db.add_submission(submission_id)
    submission_diff, donors_diff = db.diff(submission_id, metadata, submission_uploaded_date=None)
    db.commit_changes(submission_id, submission_diff, donors_diff)
    return db.get_submission(submission_id)


def test_backfill_submission_happy_path(
    db: SubmissionDb, s3_client_mock: Any, metadata: GrzSubmissionMetadata, submission_id: str
) -> None:
    """A NULL row plus valid metadata.json in S3 yields one update with size + redacted metadata persisted."""
    current = db.add_submission(submission_id)
    assert current.submission_size is None
    assert current.submission_metadata is None
    _put_metadata(s3_client_mock, submission_id, metadata)

    result = _backfill_submission(
        current_submission=current,
        s3_client=s3_client_mock,
        bucket=BUCKET,
        db_service=db,
        dry_run=False,
        force=False,
        ignore_fields=IGNORE_FIELDS,
    )

    assert result == _BackfillResult.UPDATED
    persisted = db.get_submission(submission_id)
    assert persisted.submission_size == metadata.get_submission_size()
    assert persisted.submission_metadata == metadata.to_redacted_dict()


def test_backfill_submission_returns_not_found_when_metadata_missing_in_s3(
    db: SubmissionDb, s3_client_mock: Any, metadata: GrzSubmissionMetadata, submission_id: str
) -> None:
    """When metadata.json does not exist in S3, the submission is recorded as NOT_FOUND and the row stays NULL."""
    current = db.add_submission(submission_id)

    result = _backfill_submission(
        current_submission=current,
        s3_client=s3_client_mock,
        bucket=BUCKET,
        db_service=db,
        dry_run=False,
        force=False,
        ignore_fields=IGNORE_FIELDS,
    )

    assert result == _BackfillResult.NOT_FOUND
    persisted = db.get_submission(submission_id)
    assert persisted.submission_size is None
    assert persisted.submission_metadata is None


def test_backfill_submission_records_error_on_invalid_json(
    db: SubmissionDb, s3_client_mock: Any, submission_id: str
) -> None:
    """A metadata.json that fails model_validate_json is recorded under errors and the row stays NULL."""
    current = db.add_submission(submission_id)
    s3_client_mock.put_object(
        Bucket=BUCKET,
        Key=f"{submission_id}/metadata/metadata.json",
        Body=b"{not json",
    )

    result = _backfill_submission(
        current_submission=current,
        s3_client=s3_client_mock,
        bucket=BUCKET,
        db_service=db,
        dry_run=False,
        force=False,
        ignore_fields=IGNORE_FIELDS,
    )

    assert result == _BackfillResult.ERROR
    persisted = db.get_submission(submission_id)
    assert persisted.submission_size is None
    assert persisted.submission_metadata is None


def test_backfill_submission_returns_up_to_date_when_no_pending_diff(
    db: SubmissionDb, s3_client_mock: Any, metadata: GrzSubmissionMetadata, submission_id: str
) -> None:
    """With force=True against an already-in-sync row, the diff path runs and reports UP_TO_DATE (no updates)."""
    current = _populate_full_row(db, submission_id, metadata)
    _put_metadata(s3_client_mock, submission_id, metadata)

    result = _backfill_submission(
        current_submission=current,
        s3_client=s3_client_mock,
        bucket=BUCKET,
        db_service=db,
        dry_run=False,
        force=True,
        ignore_fields=IGNORE_FIELDS,
    )

    assert result == _BackfillResult.UP_TO_DATE


def test_backfill_submission_holds_back_a_destructive_change_without_force(
    db: SubmissionDb, s3_client_mock: Any, metadata: GrzSubmissionMetadata, submission_id: str
) -> None:
    """Without --force, a differing non-NULL field is left alone while the missing ones are filled.

    Populating a field that was NULL destroys nothing, and refusing to do it because some other
    column disagrees would defeat what this command is for.
    """
    current = db.add_submission(submission_id)
    current.submission_size = 1  # non-NULL value that will differ from metadata
    db.update_submission(current)
    current = db.get_submission(submission_id)
    assert current.submission_metadata is None, "the fixture must leave something additive to write"
    _put_metadata(s3_client_mock, submission_id, metadata)

    result = _backfill_submission(
        current_submission=current,
        s3_client=s3_client_mock,
        bucket=BUCKET,
        db_service=db,
        dry_run=False,
        force=False,
        ignore_fields=IGNORE_FIELDS,
    )

    assert result == _BackfillResult.UPDATED
    persisted = db.get_submission(submission_id)
    assert persisted.submission_size == 1, "an existing value must not be overwritten"
    assert persisted.submission_metadata == metadata.to_redacted_dict(), "a NULL field must still be filled"


def test_backfill_submission_returns_would_overwrite_when_only_overwrites_are_pending(
    db: SubmissionDb, s3_client_mock: Any, metadata: GrzSubmissionMetadata, submission_id: str
) -> None:
    """When every pending change is an overwrite and none is permitted, nothing is written."""
    db.add_submission(submission_id)
    _put_metadata(s3_client_mock, submission_id, metadata)

    # bring the row up to date first, so the only difference afterwards is the one below
    _backfill_submission(
        current_submission=db.get_submission(submission_id),
        s3_client=s3_client_mock,
        bucket=BUCKET,
        db_service=db,
        dry_run=False,
        force=True,
        ignore_fields=IGNORE_FIELDS,
    )
    current = db.get_submission(submission_id)
    current.submission_size = 1
    db.update_submission(current)

    result = _backfill_submission(
        current_submission=db.get_submission(submission_id),
        s3_client=s3_client_mock,
        bucket=BUCKET,
        db_service=db,
        dry_run=False,
        force=False,
        ignore_fields=IGNORE_FIELDS,
    )

    assert result == _BackfillResult.WOULD_OVERWRITE
    assert db.get_submission(submission_id).submission_size == 1


def test_backfill_submission_force_applies_destructive_changes(
    db: SubmissionDb, s3_client_mock: Any, metadata: GrzSubmissionMetadata, submission_id: str
) -> None:
    """With --force, a submission whose non-NULL field differs from S3 is updated even though it overwrites data."""
    current = db.add_submission(submission_id)
    current.submission_size = 1  # non-NULL value that will differ from metadata
    db.update_submission(current)
    current = db.get_submission(submission_id)
    _put_metadata(s3_client_mock, submission_id, metadata)

    result = _backfill_submission(
        current_submission=current,
        s3_client=s3_client_mock,
        bucket=BUCKET,
        db_service=db,
        dry_run=False,
        force=True,
        ignore_fields=IGNORE_FIELDS,
    )

    assert result == _BackfillResult.UPDATED
    persisted = db.get_submission(submission_id)
    assert persisted.submission_size == metadata.get_submission_size()
    assert persisted.submission_metadata == metadata.to_redacted_dict()


def test_backfill_submission_force_does_not_overwrite_ignore_fields(
    db: SubmissionDb, s3_client_mock: Any, metadata: GrzSubmissionMetadata, submission_id: str
) -> None:
    """Even with --force, the hard-coded ignore_fields are not overwritten by re-derived metadata."""
    current = db.add_submission(submission_id)
    current.submission_uploaded_date = DIFFERENT_DATE
    current.tan_g = DIFFERENT_TAN_G
    current.pseudonym = DIFFERENT_PSEUDONYM
    db.update_submission(current)
    current = db.get_submission(submission_id)
    assert metadata.submission.tan_g != DIFFERENT_TAN_G
    assert metadata.submission.local_case_id != DIFFERENT_PSEUDONYM
    assert metadata.submission.submission_date != DIFFERENT_DATE
    _put_metadata(s3_client_mock, submission_id, metadata)

    result = _backfill_submission(
        current_submission=current,
        s3_client=s3_client_mock,
        bucket=BUCKET,
        db_service=db,
        dry_run=False,
        force=True,
        ignore_fields=IGNORE_FIELDS,
    )

    assert result == _BackfillResult.UPDATED
    persisted = db.get_submission(submission_id)
    assert persisted.submission_uploaded_date == DIFFERENT_DATE
    assert persisted.tan_g == DIFFERENT_TAN_G
    assert persisted.pseudonym == DIFFERENT_PSEUDONYM
    assert persisted.submission_size == metadata.get_submission_size()
    assert persisted.submission_metadata is not None
    assert persisted.submission_metadata == metadata.to_redacted_dict()


def test_backfill_submission_reads_a_consent_datetime_without_a_timezone(
    db: SubmissionDb, s3_client_mock: Any, metadata: GrzSubmissionMetadata, submission_id: str
) -> None:
    """A submission accepted before FHIR's timezone rule was enforced is read, not skipped.

    The models keep a consent dateTime as submitted, so such a document parses as it stands and
    backfill needs no repair step. What it stores is what the submitter wrote.
    """
    current = db.add_submission(submission_id)
    raw = json.loads(metadata.model_dump_json(by_alias=True))
    scope = raw["donors"][0]["researchConsents"][0]["scope"]
    scope["dateTime"] = "2020-09-01T14:37:22"  # as an older submission would have stated it
    s3_client_mock.put_object(
        Bucket=BUCKET,
        Key=f"{submission_id}/metadata/metadata.json",
        Body=json.dumps(raw).encode("utf-8"),
    )

    result = _backfill_submission(
        current_submission=current,
        s3_client=s3_client_mock,
        bucket=BUCKET,
        db_service=db,
        dry_run=False,
        force=False,
        ignore_fields=IGNORE_FIELDS,
    )

    assert result == _BackfillResult.UPDATED
    stored = db.get_submission(submission_id).submission_metadata
    stored_scope = stored["donors"][0]["researchConsents"][0]["scope"]
    assert stored_scope["dateTime"] == "2020-09-01T14:37:22", "the submitted value must be stored as written"


def test_backfill_submission_allow_overwrite_writes_only_the_named_field(
    db: SubmissionDb, s3_client_mock: Any, metadata: GrzSubmissionMetadata, submission_id: str
) -> None:
    """--allow-overwrite overwrites the field it names and holds every other overwrite back.

    Refreshing stored submission_metadata is the motivating case: it has to be rewritten from S3
    without --force also reverting a manually corrected value in some other column.
    """
    current = db.add_submission(submission_id)
    current.submission_size = 1  # non-NULL, differs from metadata, and is NOT allowed to change
    current.submission_metadata = {"stale": True}  # non-NULL, differs, and IS allowed to change
    db.update_submission(current)
    current = db.get_submission(submission_id)
    _put_metadata(s3_client_mock, submission_id, metadata)

    result = _backfill_submission(
        current_submission=current,
        s3_client=s3_client_mock,
        bucket=BUCKET,
        db_service=db,
        dry_run=False,
        force=False,
        ignore_fields=IGNORE_FIELDS,
        allow_overwrite=frozenset({"submission_metadata"}),
    )

    assert result == _BackfillResult.UPDATED
    persisted = db.get_submission(submission_id)
    assert persisted.submission_metadata == metadata.to_redacted_dict(), "the named field must be refreshed"
    assert persisted.submission_size == 1, "an overwrite outside --allow-overwrite must be held back"


def test_backfill_submission_allow_overwrite_reports_would_overwrite_when_nothing_is_writable(
    db: SubmissionDb, s3_client_mock: Any, metadata: GrzSubmissionMetadata, submission_id: str
) -> None:
    """When every pending change is an overwrite the allow-list does not cover, nothing is written."""
    db.add_submission(submission_id)
    _put_metadata(s3_client_mock, submission_id, metadata)

    # Bring the row fully up to date first, so the only pending change afterwards is the one below.
    # A freshly added submission still has many NULL columns, and filling those is additive.
    _backfill_submission(
        current_submission=db.get_submission(submission_id),
        s3_client=s3_client_mock,
        bucket=BUCKET,
        db_service=db,
        dry_run=False,
        force=True,
        ignore_fields=IGNORE_FIELDS,
    )

    current = db.get_submission(submission_id)
    current.submission_size = 1  # now the only difference, and not in the allow-list
    db.update_submission(current)

    result = _backfill_submission(
        current_submission=db.get_submission(submission_id),
        s3_client=s3_client_mock,
        bucket=BUCKET,
        db_service=db,
        dry_run=False,
        force=False,
        ignore_fields=IGNORE_FIELDS,
        allow_overwrite=frozenset({"submission_metadata"}),
    )

    assert result == _BackfillResult.WOULD_OVERWRITE
    assert db.get_submission(submission_id).submission_size == 1
