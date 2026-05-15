"""
Tests for SubmissionDb.populate.

Exercises the populate orchestration on grz-db directly. The S3 last-modified
date and the parsed metadata are passed in as arguments, so neither S3 nor
filesystem I/O is involved in these tests.
"""

import json
from datetime import date
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest
from grz_db.models.submission import SubmissionDb
from grz_db.models.submission.diff import DiffState, DonorDiff, DonorsDiffCollection, SubmissionDiffCollection
from grz_pydantic_models.submission.metadata import GrzSubmissionMetadata
from grzctl.models.config import DbConfig

SUBMISSION_DATE = date(2025, 9, 15)


def _parse(metadata_raw: dict) -> GrzSubmissionMetadata:
    return GrzSubmissionMetadata.model_validate_json(json.dumps(metadata_raw))


@pytest.fixture
def db_ctx(blank_database_config_path: Path, test_metadata_path: Path) -> SimpleNamespace:
    """SubmissionDb + parsed metadata wired up for populate tests.

    The submission is registered in the database (``db.add_submission``) but not
    yet populated, so every field starts as NULL.
    """
    metadata_raw = json.loads(test_metadata_path.read_text())
    metadata = _parse(metadata_raw)
    submission_id = metadata.submission_id
    config = DbConfig.from_path(blank_database_config_path)
    db = SubmissionDb(db_url=config.db.database_url, author=None)
    db.add_submission(submission_id)
    return SimpleNamespace(db=db, metadata=metadata, metadata_raw=metadata_raw, submission_id=submission_id)


def test_populate_no_raise_on_additive_changes(db_ctx: SimpleNamespace):
    """No RuntimeError when all changes are additive (no updates or deletions).

    A freshly-added submission has NULL for every field, so the first populate
    creates only NEW/ADDED diffs and both destructive-change guards are skipped
    regardless of force.
    """
    ctx = db_ctx
    ctx.db.populate(ctx.submission_id, ctx.metadata, SUBMISSION_DATE, force=False)

    submission = ctx.db.get_submission(ctx.submission_id)
    assert submission.submitter_id == ctx.metadata.submission.submitter_id
    assert submission.pseudonym == ctx.metadata.submission.local_case_id
    assert len(ctx.db.get_donors(ctx.submission_id)) == len(ctx.metadata.donors)


def test_populate_no_raise_when_already_up_to_date(db_ctx: SimpleNamespace):
    """No RuntimeError and no changes are written when the database is already
    in sync with the metadata (all diffs UNCHANGED).
    """
    ctx = db_ctx
    ctx.db.populate(ctx.submission_id, ctx.metadata, SUBMISSION_DATE, force=False)

    # Second call with identical metadata: all diffs UNCHANGED, nothing to commit.
    ctx.db.populate(ctx.submission_id, ctx.metadata, SUBMISSION_DATE, force=False)

    submission = ctx.db.get_submission(ctx.submission_id)
    assert submission.submitter_id == ctx.metadata.submission.submitter_id
    assert len(ctx.db.get_donors(ctx.submission_id)) == len(ctx.metadata.donors)


def test_populate_raises_without_force_on_submission_update(db_ctx: SimpleNamespace):
    """RuntimeError when a submission-level field would be updated and force=False."""
    ctx = db_ctx
    ctx.db.populate(ctx.submission_id, ctx.metadata, SUBMISSION_DATE, force=True)

    # Changing submitterId creates an UPDATE diff for the submitter_id field.
    ctx.metadata_raw["submission"]["submitterId"] = "999999999"
    mutated_metadata = _parse(ctx.metadata_raw)

    with pytest.raises(RuntimeError, match="submission data"):
        ctx.db.populate(ctx.submission_id, mutated_metadata, SUBMISSION_DATE, force=False)


def test_populate_raises_without_force_on_donor_deletion(db_ctx: SimpleNamespace):
    """RuntimeError when donors would be deleted and force=False.

    ``db.diff`` is mocked so the submission diff is clean but one donor is
    deleted, isolating the donors guard.
    """
    ctx = db_ctx
    clean_submission_diff = SubmissionDiffCollection()
    donors_diff_with_deletion = DonorsDiffCollection()
    donors_diff_with_deletion.deleted.append(
        DonorDiff(before=MagicMock(), after=None, state=DiffState.DELETED, pseudonym="deleted_donor")
    )

    with patch.object(ctx.db, "diff", return_value=(clean_submission_diff, donors_diff_with_deletion)):
        with pytest.raises(RuntimeError, match="donors"):
            ctx.db.populate(ctx.submission_id, ctx.metadata, SUBMISSION_DATE, force=False)


def test_populate_force_commits_destructive_changes(db_ctx: SimpleNamespace):
    """force=True allows submission-level updates and donor renames to be committed.

    A donor's pseudonym is renamed (old pseudonym deleted, new pseudonym added)
    to create a ``donors_diff.deleted`` entry without reducing the total donor
    count below the duo minimum.
    """
    ctx = db_ctx
    assert len(ctx.metadata_raw["donors"]) >= 2, "Test requires at least 2 donors in the fixture metadata"

    original_pseudonym = ctx.metadata_raw["donors"][1]["donorPseudonym"]
    renamed_pseudonym = "renamed_donor_for_force_test"
    new_submitter_id = "999999999"

    ctx.db.populate(ctx.submission_id, ctx.metadata, SUBMISSION_DATE, force=True)

    # Rename one donor (DELETE old pseudonym + ADD new pseudonym) and update submitter ID.
    ctx.metadata_raw["submission"]["submitterId"] = new_submitter_id
    ctx.metadata_raw["donors"][1]["donorPseudonym"] = renamed_pseudonym
    mutated_metadata = _parse(ctx.metadata_raw)

    ctx.db.populate(ctx.submission_id, mutated_metadata, SUBMISSION_DATE, force=True)

    submission = ctx.db.get_submission(ctx.submission_id)
    assert submission.submitter_id == new_submitter_id
    donor_pseudonyms = {d.pseudonym for d in ctx.db.get_donors(ctx.submission_id)}
    assert original_pseudonym not in donor_pseudonyms, "Old donor pseudonym should be removed"
    assert renamed_pseudonym in donor_pseudonyms, "Renamed donor pseudonym should be present"
    assert len(donor_pseudonyms) == 2
