# Backfill unit tests — see .vbw-planning/phases/02-implement-grzctl-db-backfill-command/02-02-PLAN.md
"""Unit tests for `_backfill_submission`.

The tests target `_backfill_submission` directly with a real SubmissionDb on SQLite
and a moto-mocked S3. Live in `grzctl/tests/` (not `grz-db/tests/`) because `grzctl`
is not a dev/test dependency of `grz-db`, but `moto[s3]`, `pytest-postgresql`, and
`grz-pydantic-models-testing` are all in `grzctl`'s [test] dependency group.
"""

import datetime
import importlib.resources
import json
from collections.abc import Iterator
from typing import Any
from unittest.mock import MagicMock

import boto3
import pytest
from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives.asymmetric import ed25519
from grz_db.models.author import Author
from grz_db.models.submission import Submission, SubmissionDb
from grz_pydantic_models.submission.metadata import GrzSubmissionMetadata
from grz_pydantic_models_testing.example_metadata import grzctl as grzctl_metadata
from grzctl.commands.db.cli import _BackfillResult, _backfill_submission
from moto import mock_aws

BUCKET = "test-backfill-bucket"
REGION = "us-east-1"
IGNORE_FIELDS = {"submission_date", "tan_g", "local_case_id"}
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
def test_author() -> Author:
    key = ed25519.Ed25519PrivateKey.generate()
    private_key_bytes = key.private_bytes(
        encoding=serialization.Encoding.PEM,
        format=serialization.PrivateFormat.OpenSSH,
        encryption_algorithm=serialization.NoEncryption(),
    )
    return Author(name="alice", private_key_bytes=private_key_bytes, private_key_passphrase="")


@pytest.fixture
def db(tmp_path, test_author: Author) -> SubmissionDb:
    db_url = f"sqlite:///{tmp_path / 'backfill.db'}"
    submission_db = SubmissionDb(db_url=db_url, author=test_author, debug=False)
    submission_db.initialize_schema()
    return submission_db


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
    submission_diff, donors_diff = db.diff(submission_id, metadata, submission_date=None)
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
    result = _BackfillResult()

    _backfill_submission(
        current_submission=current,
        s3_client=s3_client_mock,
        bucket=BUCKET,
        db_service=db,
        dry_run=False,
        force=False,
        ignore_fields=IGNORE_FIELDS,
        result=result,
    )

    assert result.updated == 1
    assert result.skipped == 0
    assert result.errors == []
    persisted = db.get_submission(submission_id)
    assert persisted.submission_size == metadata.get_submission_size()
    assert persisted.submission_metadata == metadata.to_redacted_dict()


def test_backfill_submission_skips_when_metadata_missing_in_s3(
    db: SubmissionDb, s3_client_mock: Any, metadata: GrzSubmissionMetadata, submission_id: str
) -> None:
    """When metadata.json does not exist in S3, the submission is recorded as skipped and the row stays NULL."""
    current = db.add_submission(submission_id)
    result = _BackfillResult()

    _backfill_submission(
        current_submission=current,
        s3_client=s3_client_mock,
        bucket=BUCKET,
        db_service=db,
        dry_run=False,
        force=False,
        ignore_fields=IGNORE_FIELDS,
        result=result,
    )

    assert result.updated == 0
    assert result.skipped == 1
    assert result.errors == []
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
    result = _BackfillResult()

    _backfill_submission(
        current_submission=current,
        s3_client=s3_client_mock,
        bucket=BUCKET,
        db_service=db,
        dry_run=False,
        force=False,
        ignore_fields=IGNORE_FIELDS,
        result=result,
    )

    assert result.updated == 0
    assert result.skipped == 0
    assert len(result.errors) == 1
    failed_id, _exc = result.errors[0]
    assert failed_id == submission_id
    persisted = db.get_submission(submission_id)
    assert persisted.submission_size is None
    assert persisted.submission_metadata is None


def test_backfill_submission_skips_when_no_pending_diff(
    db: SubmissionDb, s3_client_mock: Any, metadata: GrzSubmissionMetadata, submission_id: str
) -> None:
    """With force=True against an already-in-sync row, the diff path runs and reports skipped (no updates)."""
    current = _populate_full_row(db, submission_id, metadata)
    _put_metadata(s3_client_mock, submission_id, metadata)
    result = _BackfillResult()

    _backfill_submission(
        current_submission=current,
        s3_client=s3_client_mock,
        bucket=BUCKET,
        db_service=db,
        dry_run=False,
        force=True,
        ignore_fields=IGNORE_FIELDS,
        result=result,
    )

    assert result.updated == 0
    assert result.skipped == 1
    assert result.errors == []


def test_backfill_submission_pre_filter_skips_when_both_columns_set_and_not_force(
    db: SubmissionDb, s3_client_mock: Any, metadata: GrzSubmissionMetadata, submission_id: str
) -> None:
    """Without --force, an already-populated row is short-circuited before any S3 access."""
    current = _populate_full_row(db, submission_id, metadata)
    spy_client = MagicMock(wraps=s3_client_mock)
    result = _BackfillResult()

    _backfill_submission(
        current_submission=current,
        s3_client=spy_client,
        bucket=BUCKET,
        db_service=db,
        dry_run=False,
        force=False,
        ignore_fields=IGNORE_FIELDS,
        result=result,
    )

    assert result.updated == 0
    assert result.skipped == 1
    assert result.errors == []
    assert spy_client.get_object.call_count == 0


def test_backfill_submission_force_does_not_overwrite_ignore_fields(
    db: SubmissionDb, s3_client_mock: Any, metadata: GrzSubmissionMetadata, submission_id: str
) -> None:
    """Even with --force, the hard-coded ignore_fields are not overwritten by re-derived metadata."""
    current = db.add_submission(submission_id)
    current.submission_date = DIFFERENT_DATE
    current.tan_g = DIFFERENT_TAN_G
    current.pseudonym = DIFFERENT_PSEUDONYM
    db.update_submission(current)
    current = db.get_submission(submission_id)
    assert metadata.submission.tan_g != DIFFERENT_TAN_G
    assert metadata.submission.local_case_id != DIFFERENT_PSEUDONYM
    assert metadata.submission.submission_date != DIFFERENT_DATE
    _put_metadata(s3_client_mock, submission_id, metadata)
    result = _BackfillResult()

    _backfill_submission(
        current_submission=current,
        s3_client=s3_client_mock,
        bucket=BUCKET,
        db_service=db,
        dry_run=False,
        force=True,
        ignore_fields=IGNORE_FIELDS,
        result=result,
    )

    assert result.errors == []
    persisted = db.get_submission(submission_id)
    assert persisted.submission_date == DIFFERENT_DATE
    assert persisted.tan_g == DIFFERENT_TAN_G
    assert persisted.pseudonym == DIFFERENT_PSEUDONYM
    assert persisted.submission_size == metadata.get_submission_size()
    assert persisted.submission_metadata is not None
    assert persisted.submission_metadata == metadata.to_redacted_dict()
