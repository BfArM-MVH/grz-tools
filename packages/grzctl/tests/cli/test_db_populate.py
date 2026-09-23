"""
Tests for SubmissionDb.populate and for ``grzctl db submission populate``.

The grz-db half exercises the populate orchestration directly. The S3 last-modified
date and the parsed metadata are passed in as arguments, so neither S3 nor
filesystem I/O is involved in those tests. The command half runs the CLI against the
same database, since the two refuse a destructive change on their own paths. The
command reads the upload date from a moto inbox unless --submission-date is given.
"""

import json
from collections.abc import Iterator
from datetime import date
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import boto3
import click.testing
import grzctl.cli
import pytest
import yaml
from grz_db.errors import SubmissionNotFoundError
from grz_db.models.submission import DONORS_KEY, SubmissionDb
from grz_db.models.submission.diff import DiffState, DonorDiff, SubmissionChangeSet
from grz_pydantic_models.submission.metadata import REDACTED_LOCAL_CASE_ID, REDACTED_TAN, GrzSubmissionMetadata
from grzctl.models.config import GrzctlConfig
from moto import mock_aws

SUBMISSION_DATE = date(2025, 9, 15)
INBOX_BUCKET = "inbox"
REGION = "us-east-1"


def _parse(metadata_raw: dict) -> GrzSubmissionMetadata:
    return GrzSubmissionMetadata.model_validate_json(json.dumps(metadata_raw))


def _db_ctx(config_path: Path, test_metadata_path: Path, *, register: bool) -> SimpleNamespace:
    metadata_raw = json.loads(test_metadata_path.read_text())
    metadata = _parse(metadata_raw)
    submission_id = metadata.submission_id
    config = GrzctlConfig.from_path(config_path)
    db = SubmissionDb(db_url=config.db.database_url, author=None)
    if register:
        db.add_submission(submission_id)
    return SimpleNamespace(db=db, metadata=metadata, metadata_raw=metadata_raw, submission_id=submission_id)


@pytest.fixture
def db_ctx(migrated_database_config_path: Path, test_metadata_path: Path) -> SimpleNamespace:
    """SubmissionDb + parsed metadata wired up for populate tests.

    The submission is registered in the database (``db.add_submission``) but not
    yet populated, so every field starts as NULL.
    """
    return _db_ctx(migrated_database_config_path, test_metadata_path, register=True)


@pytest.fixture
def unregistered_db_ctx(migrated_database_config_path: Path, test_metadata_path: Path) -> SimpleNamespace:
    """Like :func:`db_ctx`, but the submission was never added, as after a fresh download."""
    return _db_ctx(migrated_database_config_path, test_metadata_path, register=False)


def test_populate_registers_an_unknown_submission_when_asked(unregistered_db_ctx: SimpleNamespace):
    """``download`` populates submissions it has just fetched, which are not yet in the database."""
    ctx = unregistered_db_ctx

    ctx.db.populate(ctx.submission_id, ctx.metadata, SUBMISSION_DATE, on_missing="create")

    submission = ctx.db.get_submission(ctx.submission_id)
    assert submission is not None
    assert submission.submitter_id == ctx.metadata.submission.submitter_id
    assert len(ctx.db.get_donors(ctx.submission_id)) == len(ctx.metadata.donors)


def test_populate_refuses_an_unknown_submission_by_default(unregistered_db_ctx: SimpleNamespace):
    """Defaulting to an error keeps a mistyped ID from quietly creating a second submission."""
    ctx = unregistered_db_ctx

    with pytest.raises(SubmissionNotFoundError):
        ctx.db.populate(ctx.submission_id, ctx.metadata, SUBMISSION_DATE)

    assert ctx.db.get_submission(ctx.submission_id) is None


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
    assert submission.local_case_id == ctx.metadata.submission.local_case_id
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

    ``db.diff`` is mocked so the field diffs are clean but one donor is
    deleted, isolating the donor part of the destructive guard.
    """
    ctx = db_ctx
    changes = SubmissionChangeSet()
    changes.donors.deleted.append(
        DonorDiff(before=MagicMock(), after=None, state=DiffState.DELETED, pseudonym="deleted_donor")
    )

    with patch.object(ctx.db, "diff", return_value=changes):
        with pytest.raises(RuntimeError, match="donor 'deleted_donor'"):
            ctx.db.populate(ctx.submission_id, ctx.metadata, SUBMISSION_DATE, force=False)


def _flip_a_donors_mv_consent(ctx: SimpleNamespace) -> bool:
    """Make one stored donor differ from the metadata, and return the value the metadata still has."""
    donor = ctx.db.get_donors(ctx.submission_id)[0]
    from_metadata = donor.mv_consented
    donor.mv_consented = not from_metadata
    ctx.db.update_donor(donor)
    return from_metadata


def test_populate_leaves_the_donors_alone_when_they_are_ignored(db_ctx: SimpleNamespace):
    """Ignoring the donors is how a caller writes the rest without deciding the donor question."""
    ctx = db_ctx
    ctx.db.populate(ctx.submission_id, ctx.metadata, SUBMISSION_DATE, force=True)
    from_metadata = _flip_a_donors_mv_consent(ctx)

    ctx.db.populate(ctx.submission_id, ctx.metadata, SUBMISSION_DATE, ignore_fields={DONORS_KEY})

    assert ctx.db.get_donors(ctx.submission_id)[0].mv_consented == (not from_metadata)


def test_populate_writes_a_donor_overwrite_that_allow_overwrite_names(db_ctx: SimpleNamespace):
    """Naming the donors is how a caller settles it the other way, without permitting everything."""
    ctx = db_ctx
    ctx.db.populate(ctx.submission_id, ctx.metadata, SUBMISSION_DATE, force=True)
    from_metadata = _flip_a_donors_mv_consent(ctx)

    ctx.db.populate(ctx.submission_id, ctx.metadata, SUBMISSION_DATE, allow_overwrite={DONORS_KEY})

    assert ctx.db.get_donors(ctx.submission_id)[0].mv_consented == from_metadata


def test_populate_raises_on_redacted_tan_g(db_ctx: SimpleNamespace):
    """ValueError when ``tan_g`` is redacted and ``"tan_g"`` is not in ``ignore_fields``."""
    ctx = db_ctx
    ctx.metadata_raw["submission"]["tanG"] = REDACTED_TAN
    redacted_metadata = _parse(ctx.metadata_raw)

    with pytest.raises(ValueError, match="redacted tan_g"):
        ctx.db.populate(ctx.submission_id, redacted_metadata, SUBMISSION_DATE, force=True)


def test_populate_skips_redacted_tan_g_when_ignored(db_ctx: SimpleNamespace):
    """``ignore_fields={"tan_g"}`` bypasses the redacted-TAN guard."""
    ctx = db_ctx
    ctx.metadata_raw["submission"]["tanG"] = REDACTED_TAN
    redacted_metadata = _parse(ctx.metadata_raw)

    ctx.db.populate(
        ctx.submission_id,
        redacted_metadata,
        SUBMISSION_DATE,
        force=True,
        ignore_fields={"tan_g"},
    )

    submission = ctx.db.get_submission(ctx.submission_id)
    assert submission.local_case_id == ctx.metadata.submission.local_case_id


def test_populate_raises_on_redacted_local_case_id(db_ctx: SimpleNamespace):
    """ValueError when ``local_case_id`` is redacted and ``"local_case_id"`` is not in ``ignore_fields``."""
    ctx = db_ctx
    ctx.metadata_raw["submission"]["localCaseId"] = REDACTED_LOCAL_CASE_ID
    redacted_metadata = _parse(ctx.metadata_raw)

    with pytest.raises(ValueError, match="local_case_id"):
        ctx.db.populate(ctx.submission_id, redacted_metadata, SUBMISSION_DATE, force=True)


def test_populate_skips_redacted_local_case_id_when_ignored(db_ctx: SimpleNamespace):
    """``ignore_fields={"local_case_id"}`` bypasses the missing/redacted local_case_id guard."""
    ctx = db_ctx
    ctx.metadata_raw["submission"]["localCaseId"] = REDACTED_LOCAL_CASE_ID
    redacted_metadata = _parse(ctx.metadata_raw)

    ctx.db.populate(
        ctx.submission_id,
        redacted_metadata,
        SUBMISSION_DATE,
        force=True,
        ignore_fields={"local_case_id"},
    )

    submission = ctx.db.get_submission(ctx.submission_id)
    assert submission.submitter_id == ctx.metadata.submission.submitter_id


def test_populate_forwards_ignore_fields_to_diff(db_ctx: SimpleNamespace):
    """``ignore_fields`` flows through to :meth:`SubmissionDb.diff`."""
    ctx = db_ctx

    with patch.object(ctx.db, "diff", wraps=ctx.db.diff) as diff_spy:
        ctx.db.populate(
            ctx.submission_id,
            ctx.metadata,
            SUBMISSION_DATE,
            force=False,
            ignore_fields={"submitter_id"},
        )

    diff_spy.assert_called_once()
    _, kwargs = diff_spy.call_args
    assert kwargs["ignore_fields"] == {"submitter_id"}


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


def _write_metadata(tmp_path: Path, metadata_raw: dict) -> Path:
    """Write *metadata_raw* where the command can read it."""
    path = tmp_path / "mutated.metadata.json"
    path.write_text(json.dumps(metadata_raw))
    return path


def _invoke_populate(
    config_path: Path,
    submission_id: str,
    metadata_path: Path,
    *args: str,
    submission_date: date | None = SUBMISSION_DATE,
):
    """Invoke ``grzctl db submission populate`` against the database *config_path* names.

    With *submission_date* set to ``None``, the command looks for the upload date in the inbox.
    """
    date_args = ["--submission-date", submission_date.isoformat()] if submission_date is not None else []
    runner = click.testing.CliRunner()
    return runner.invoke(
        grzctl.cli.build_cli(),
        [
            "--config",
            str(config_path),
            "db",
            "submission",
            "populate",
            submission_id,
            str(metadata_path),
            *date_args,
            *args,
        ],
    )


@pytest.fixture
def inbox_config_path(
    migrated_database_config: GrzctlConfig, tmp_path: Path, test_metadata_path: Path
) -> Iterator[Path]:
    """A config whose only inbox for the example submitter is an empty moto bucket, which lives as long as the test."""
    submitter_id = GrzSubmissionMetadata.model_validate_json(test_metadata_path.read_text()).submission.submitter_id
    data = migrated_database_config.model_dump(mode="json", exclude_none=True, context={"reveal_secrets": True})
    data["leistungserbringer"] = {
        submitter_id: {"inbox_buckets": {INBOX_BUCKET: {"private_key_path": "/dev/null", "region_name": REGION}}}
    }
    config_path = tmp_path / "config.inbox.yaml"
    config_path.write_text(yaml.safe_dump(data))
    with mock_aws():
        boto3.client("s3", region_name=REGION).create_bucket(Bucket=INBOX_BUCKET)
        yield config_path


def _put_inbox_metadata(submission_id: str, body: bytes) -> date:
    """Write metadata.json into the inbox, and return the date that S3 records for it."""
    s3_client = boto3.client("s3", region_name=REGION)
    key = f"{submission_id}/metadata/metadata.json"
    s3_client.put_object(Bucket=INBOX_BUCKET, Key=key, Body=body)
    return s3_client.head_object(Bucket=INBOX_BUCKET, Key=key)["LastModified"].date()


def test_populate_command_takes_the_upload_date_from_the_inbox(
    db_ctx: SimpleNamespace, inbox_config_path: Path, test_metadata_path: Path
):
    """The date is when metadata.json arrived in the inbox, not the submissionDate it contains."""
    ctx = db_ctx
    uploaded = _put_inbox_metadata(ctx.submission_id, test_metadata_path.read_bytes())
    assert uploaded != ctx.metadata.submission.submission_date

    result = _invoke_populate(
        inbox_config_path, ctx.submission_id, test_metadata_path, "--no-confirm", submission_date=None
    )

    assert result.exit_code == 0, result.stderr
    assert ctx.db.get_submission(ctx.submission_id).submission_uploaded_date == uploaded


@pytest.mark.parametrize(
    ("body", "marker"),
    [(b"", "cleaned"), (b"{}", "cleaning"), (None, None)],
    ids=["cleaned", "being-cleaned", "missing"],
)
def test_populate_command_needs_a_date_when_the_inbox_has_none(
    db_ctx: SimpleNamespace,
    inbox_config_path: Path,
    test_metadata_path: Path,
    body: bytes | None,
    marker: str | None,
):
    """``grzctl clean`` leaves an empty metadata.json, whose LastModified is the time of cleaning."""
    ctx = db_ctx
    if body is not None:
        _put_inbox_metadata(ctx.submission_id, body)
    if marker is not None:
        s3_client = boto3.client("s3", region_name=REGION)
        s3_client.put_object(Bucket=INBOX_BUCKET, Key=f"{ctx.submission_id}/{marker}", Body=b"")

    result = _invoke_populate(
        inbox_config_path, ctx.submission_id, test_metadata_path, "--no-confirm", submission_date=None
    )

    assert result.exit_code != 0
    assert "Pass --submission-date" in result.stderr
    assert ctx.db.get_submission(ctx.submission_id).local_case_id is None, "nothing is written"


def test_populate_command_needs_a_date_even_when_one_is_stored(
    db_ctx: SimpleNamespace, migrated_database_config_path: Path, test_metadata_path: Path
):
    """Without an inbox for the submitter, a re-populate needs --submission-date again."""
    ctx = db_ctx
    ctx.db.populate(ctx.submission_id, ctx.metadata, SUBMISSION_DATE, force=True)

    result = _invoke_populate(
        migrated_database_config_path, ctx.submission_id, test_metadata_path, "--no-confirm", submission_date=None
    )

    assert result.exit_code != 0
    assert "has no inbox in the configuration" in result.stderr
    assert ctx.db.get_submission(ctx.submission_id).submission_uploaded_date == SUBMISSION_DATE


def test_populate_command_refuses_an_overwrite_and_writes_nothing(
    db_ctx: SimpleNamespace, migrated_database_config_path: Path, tmp_path: Path
):
    """A value the database already holds stops the command, and the rest is not written either."""
    ctx = db_ctx
    ctx.db.populate(ctx.submission_id, ctx.metadata, SUBMISSION_DATE, force=True, ignore_fields={"disease_type"})
    stored_submitter_id = ctx.metadata.submission.submitter_id

    ctx.metadata_raw["submission"]["submitterId"] = "999999999"
    metadata_path = _write_metadata(tmp_path, ctx.metadata_raw)

    result = _invoke_populate(migrated_database_config_path, ctx.submission_id, metadata_path, "--no-confirm")

    assert result.exit_code != 0
    assert "Refusing to overwrite or remove" in result.stderr
    submission = ctx.db.get_submission(ctx.submission_id)
    assert submission.submitter_id == stored_submitter_id
    assert submission.disease_type is None, "the NULL the run could have filled is not written either"


def test_populate_command_force_writes_the_overwrite(
    db_ctx: SimpleNamespace, migrated_database_config_path: Path, tmp_path: Path
):
    ctx = db_ctx
    ctx.db.populate(ctx.submission_id, ctx.metadata, SUBMISSION_DATE, force=True)

    ctx.metadata_raw["submission"]["submitterId"] = "999999999"
    metadata_path = _write_metadata(tmp_path, ctx.metadata_raw)

    result = _invoke_populate(
        migrated_database_config_path, ctx.submission_id, metadata_path, "--force", "--no-confirm"
    )

    assert result.exit_code == 0, result.stderr
    assert ctx.db.get_submission(ctx.submission_id).submitter_id == "999999999"


def test_populate_command_allow_overwrite_writes_the_fields_it_names(
    db_ctx: SimpleNamespace, migrated_database_config_path: Path, tmp_path: Path
):
    """The stored metadata.json dump carries submitterId too, so both have to be named."""
    ctx = db_ctx
    ctx.db.populate(ctx.submission_id, ctx.metadata, SUBMISSION_DATE, force=True)

    ctx.metadata_raw["submission"]["submitterId"] = "999999999"
    metadata_path = _write_metadata(tmp_path, ctx.metadata_raw)

    result = _invoke_populate(
        migrated_database_config_path,
        ctx.submission_id,
        metadata_path,
        "--allow-overwrite",
        "submitter_id",
        "--allow-overwrite",
        "submission_metadata",
        "--no-confirm",
    )

    assert result.exit_code == 0, result.stderr
    assert ctx.db.get_submission(ctx.submission_id).submitter_id == "999999999"


def test_populate_command_allow_overwrite_releases_donors(
    db_ctx: SimpleNamespace, migrated_database_config_path: Path, tmp_path: Path
):
    """``donors`` releases a rename, which deletes one donor row and adds another."""
    ctx = db_ctx
    ctx.db.populate(ctx.submission_id, ctx.metadata, SUBMISSION_DATE, force=True)
    original_pseudonym = ctx.metadata_raw["donors"][1]["donorPseudonym"]

    ctx.metadata_raw["donors"][1]["donorPseudonym"] = "renamed_donor"
    metadata_path = _write_metadata(tmp_path, ctx.metadata_raw)

    result = _invoke_populate(
        migrated_database_config_path,
        ctx.submission_id,
        metadata_path,
        "--allow-overwrite",
        "donors",
        "--allow-overwrite",
        "submission_metadata",
        "--no-confirm",
    )

    assert result.exit_code == 0, result.stderr
    pseudonyms = {d.pseudonym for d in ctx.db.get_donors(ctx.submission_id)}
    assert "renamed_donor" in pseudonyms
    assert original_pseudonym not in pseudonyms
