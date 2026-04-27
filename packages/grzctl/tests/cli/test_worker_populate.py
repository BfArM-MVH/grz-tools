"""
Tests for Worker.populate

These are integration-style tests that exercise Worker.populate directly
(not via the grzctl CLI). The S3 client is mocked; the database uses the
same SQLite/PostgreSQL fixtures as the CLI tests.
"""

import json
from datetime import UTC, date, datetime
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest
from grz_common.models.s3 import S3Options
from grz_common.workers.worker import Worker
from grz_db.models.submission import SubmissionDb
from grz_db.models.submission.diff import DiffState, DonorDiff, DonorsDiffCollection, SubmissionDiffCollection
from grz_pydantic_models.submission.metadata import GrzSubmissionMetadata
from grzctl.models.config import DbConfig

# S3 endpoint is the same across every test - declare it once.
# model_construct skips Pydantic validation; endpoint_url is never
# accessed because init_s3_client is fully mocked in every test.
_S3_OPTIONS = S3Options.model_construct(endpoint_url=None, bucket="test-bucket")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_populate_worker(tmp_path: Path, blank_database_config_path: Path, metadata_raw: dict):
    """Return a (Worker, SubmissionDb) pair with *metadata_raw* written to disk."""
    metadata_dir = tmp_path / "metadata"
    metadata_dir.mkdir(exist_ok=True)
    (metadata_dir / "metadata.json").write_text(json.dumps(metadata_raw))

    worker = Worker(
        metadata_dir=metadata_dir,
        files_dir=tmp_path / "files",
        log_dir=tmp_path / "logs",
        encrypted_files_dir=tmp_path / "encrypted",
    )
    config = DbConfig.from_path(blank_database_config_path)
    db = SubmissionDb(db_url=config.db.database_url, author=None)
    return worker, db


def _fake_s3_client(submission_date: date = date(2025, 9, 15)) -> MagicMock:
    """Mock boto3 client whose head_object returns a fixed LastModified date."""
    mock = MagicMock()
    mock.head_object.return_value = {
        "LastModified": datetime(submission_date.year, submission_date.month, submission_date.day, tzinfo=UTC)
    }
    return mock


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def worker_ctx(tmp_path: Path, blank_database_config_path: Path, test_metadata_path: Path) -> SimpleNamespace:
    """Worker + SubmissionDb wired up for populate tests.

    The submission is registered in the database (``db.add_submission``) but not yet
    populated, so every field starts as NULL. Returns a ``SimpleNamespace`` with:

    - ``worker`` - :class:`Worker` instance pointing at ``tmp_path``
    - ``db`` - :class:`SubmissionDb` instance
    - ``metadata`` - parsed :class:`GrzSubmissionMetadata`
    - ``metadata_raw`` - mutable raw dict; write it back to disk to change what
      ``worker.populate`` will read next
    - ``submission_id`` - str
    """
    metadata_raw = json.loads(test_metadata_path.read_text())
    worker, db = _make_populate_worker(tmp_path, blank_database_config_path, metadata_raw)
    metadata = GrzSubmissionMetadata.model_validate_json(json.dumps(metadata_raw))
    submission_id = metadata.submission_id
    db.add_submission(submission_id)
    return SimpleNamespace(
        worker=worker, db=db, metadata=metadata, metadata_raw=metadata_raw, submission_id=submission_id
    )


@pytest.fixture
def mock_s3():
    """Patch ``init_s3_client`` so no real S3 connection is attempted."""
    with patch("grz_common.workers.worker.init_s3_client", return_value=_fake_s3_client()):
        yield


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_worker_populate_no_raise_on_additive_changes(worker_ctx: SimpleNamespace, mock_s3):
    """No RuntimeError is raised when all changes are additive (no updates or deletions).

    A freshly-added submission has NULL for every field, so the first populate creates
    only NEW/ADDED diffs. Both guards (submission and donors) are skipped regardless of
    force, and commit_changes is called to persist the data.
    """
    ctx = worker_ctx
    # force=False — must not raise because nothing is being updated or deleted
    ctx.worker.populate(_S3_OPTIONS, ctx.db, ctx.submission_id, force=False)

    submission = ctx.db.get_submission(ctx.submission_id)
    assert submission.submitter_id == ctx.metadata.submission.submitter_id
    assert submission.pseudonym == ctx.metadata.submission.local_case_id
    assert len(ctx.db.get_donors(ctx.submission_id)) == len(ctx.metadata.donors)


def test_worker_populate_no_raise_when_already_up_to_date(worker_ctx: SimpleNamespace, mock_s3):
    """No RuntimeError is raised and no changes are written when the database is already
    in sync with the metadata (all diffs are UNCHANGED).
    """
    ctx = worker_ctx
    ctx.worker.populate(_S3_OPTIONS, ctx.db, ctx.submission_id, force=False)

    # Second call with identical metadata — all diffs UNCHANGED, nothing to commit
    ctx.worker.populate(_S3_OPTIONS, ctx.db, ctx.submission_id, force=False)

    # Data is intact after the no-op second call
    submission = ctx.db.get_submission(ctx.submission_id)
    assert submission.submitter_id == ctx.metadata.submission.submitter_id
    assert len(ctx.db.get_donors(ctx.submission_id)) == len(ctx.metadata.donors)


def test_worker_populate_raises_without_force_on_submission_update(worker_ctx: SimpleNamespace, mock_s3):
    """RuntimeError is raised when a submission-level field would be updated and force=False."""
    ctx = worker_ctx
    ctx.worker.populate(_S3_OPTIONS, ctx.db, ctx.submission_id, force=True)

    # Changing submitterId creates an UPDATE diff for the submitter_id field
    ctx.metadata_raw["submission"]["submitterId"] = "999999999"
    (ctx.worker.metadata_dir / "metadata.json").write_text(json.dumps(ctx.metadata_raw))

    with pytest.raises(RuntimeError, match="submission data"):
        ctx.worker.populate(_S3_OPTIONS, ctx.db, ctx.submission_id, force=False)


def test_worker_populate_raises_without_force_on_donor_deletion(worker_ctx: SimpleNamespace, mock_s3):
    """RuntimeError is raised when donors would be deleted and force=False.

    ``db.diff`` is mocked to return a scenario where the submission diff is clean
    but one donor is deleted, isolating the donors guard on lines 402-405 of worker.py.
    Any real change to donor data would also update submission_metadata, causing the
    submission guard to fire first; mocking lets us target the donors guard directly.
    """
    ctx = worker_ctx
    clean_submission_diff = SubmissionDiffCollection()
    donors_diff_with_deletion = DonorsDiffCollection()
    donors_diff_with_deletion.deleted.append(
        DonorDiff(before=MagicMock(), after=None, state=DiffState.DELETED, pseudonym="deleted_donor")
    )

    with patch.object(ctx.db, "diff", return_value=(clean_submission_diff, donors_diff_with_deletion)):
        with pytest.raises(RuntimeError, match="donors"):
            ctx.worker.populate(_S3_OPTIONS, ctx.db, ctx.submission_id, force=False)


def test_worker_populate_force_commits_destructive_changes(worker_ctx: SimpleNamespace, mock_s3):
    """force=True allows submission-level updates and donor renames to be committed.

    A donor's pseudonym is renamed (old pseudonym deleted, new pseudonym added) to create
    a donors_diff.deleted entry without reducing the total donor count below the duo minimum.
    """
    ctx = worker_ctx
    assert len(ctx.metadata_raw["donors"]) >= 2, "Test requires at least 2 donors in the fixture metadata"

    original_pseudonym = ctx.metadata_raw["donors"][1]["donorPseudonym"]
    renamed_pseudonym = "renamed_donor_for_force_test"
    new_submitter_id = "999999999"

    ctx.worker.populate(_S3_OPTIONS, ctx.db, ctx.submission_id, force=True)

    # Rename one donor (DELETE old pseudonym + ADD new pseudonym) and update submitter ID.
    # With 2 donors remaining the duo-study constraint is still satisfied.
    ctx.metadata_raw["submission"]["submitterId"] = new_submitter_id
    ctx.metadata_raw["donors"][1]["donorPseudonym"] = renamed_pseudonym
    (ctx.worker.metadata_dir / "metadata.json").write_text(json.dumps(ctx.metadata_raw))

    ctx.worker.populate(_S3_OPTIONS, ctx.db, ctx.submission_id, force=True)

    submission = ctx.db.get_submission(ctx.submission_id)
    assert submission.submitter_id == new_submitter_id
    donor_pseudonyms = {d.pseudonym for d in ctx.db.get_donors(ctx.submission_id)}
    assert original_pseudonym not in donor_pseudonyms, "Old donor pseudonym should be removed"
    assert renamed_pseudonym in donor_pseudonyms, "Renamed donor pseudonym should be present"
    assert len(donor_pseudonyms) == 2
