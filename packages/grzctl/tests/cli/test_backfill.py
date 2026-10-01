# Backfill unit tests, see .vbw-planning/phases/02-implement-grzctl-db-backfill-command/02-02-PLAN.md
"""Tests for `grzctl db backfill`.

Most tests target `_backfill_submission` directly with a real SubmissionDb on every supported
backend. The tests at the end run the command against a moto-mocked S3 with both archives.
Live in `grzctl/tests/` (not `grz-db/tests/`) because `grzctl`
is not a dev/test dependency of `grz-db`, but `moto[s3]`, `pytest-postgresql`, and
`grz-pydantic-models-testing` are all in `grzctl`'s [test] dependency group.
"""

import datetime
import importlib.resources
import json
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import boto3
import click.testing
import grzctl.cli
import pytest
import sqlalchemy
from grz_db.models.submission import DONORS_KEY, Submission, SubmissionDb
from grz_pydantic_models.submission.metadata import GrzSubmissionMetadata
from grz_pydantic_models_testing.example_metadata import grzctl as grzctl_metadata
from grzctl.commands.db.cli import (
    _backfill_submission,
    _BackfillOutcome,
    _BackfillResult,
    _fetch_metadata_json_from_archives,
)
from moto import mock_aws

ARCHIVE_BUCKETS = ("consented", "non_consented")
"""The archive bucket names that ``tests/cli/conftest.py`` configures."""
REGION = "us-east-1"
DIFFERENT_TAN_G = "b" * 64
DIFFERENT_LOCAL_CASE_ID = "different-local-case-id"
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
    """A moto-backed S3 client without buckets, so each test decides which archives exist."""
    with mock_aws():
        yield boto3.client("s3", region_name=REGION)


def _metadata_json(metadata: GrzSubmissionMetadata) -> str:
    """Serialize *metadata* unredacted, unlike archival, which redacts first.

    For the redacted shape backfill actually reads from the archive, pass ``_archived(metadata)``.
    """
    return json.dumps(metadata.get_raw_dict())


def _put_metadata(s3_client: Any, bucket: str, submission_id: str, metadata: GrzSubmissionMetadata) -> None:
    s3_client.put_object(
        Bucket=bucket,
        Key=f"{submission_id}/metadata/metadata.json",
        Body=_metadata_json(metadata).encode("utf-8"),
    )


def _populate_full_row(db: SubmissionDb, submission_id: str, metadata: GrzSubmissionMetadata) -> Submission:
    """Persist a fully-populated row by running the same diff/commit path the production code uses."""
    db.add_submission(submission_id)
    db.commit_changes(submission_id, db.diff(submission_id, metadata, submission_uploaded_date=None))
    return db.get_submission(submission_id)


def test_backfill_submission_happy_path(db: SubmissionDb, metadata: GrzSubmissionMetadata, submission_id: str) -> None:
    """A NULL row plus valid metadata.json in S3 yields one update with size + redacted metadata persisted."""
    current = db.add_submission(submission_id)
    assert current.submission_size is None
    assert current.submission_metadata is None

    result = _backfill_submission(
        current_submission=current,
        raw_json=_metadata_json(metadata),
        db_service=db,
        dry_run=False,
        force=False,
        ignore_fields=set(),
    )

    assert result.status == _BackfillResult.UPDATED
    persisted = db.get_submission(submission_id)
    assert persisted.submission_size == metadata.get_submission_size()
    assert persisted.submission_metadata == metadata.to_redacted_dict()


def test_backfill_submission_records_error_on_invalid_json(db: SubmissionDb, submission_id: str) -> None:
    """A metadata.json that fails model_validate_json is recorded under errors and the row stays NULL."""
    current = db.add_submission(submission_id)

    result = _backfill_submission(
        current_submission=current,
        raw_json="{not json",
        db_service=db,
        dry_run=False,
        force=False,
        ignore_fields=set(),
    )

    assert result.status == _BackfillResult.ERROR
    persisted = db.get_submission(submission_id)
    assert persisted.submission_size is None
    assert persisted.submission_metadata is None


def test_backfill_submission_returns_up_to_date_when_no_pending_diff(
    db: SubmissionDb, metadata: GrzSubmissionMetadata, submission_id: str
) -> None:
    """With force=True against an already-in-sync row, the diff path runs and reports UP_TO_DATE (no updates)."""
    current = _populate_full_row(db, submission_id, metadata)

    result = _backfill_submission(
        current_submission=current,
        raw_json=_metadata_json(metadata),
        db_service=db,
        dry_run=False,
        force=True,
        ignore_fields=set(),
    )

    assert result.status == _BackfillResult.UP_TO_DATE


def test_backfill_submission_skips_a_submission_with_an_undeclared_overwrite(
    db: SubmissionDb, metadata: GrzSubmissionMetadata, submission_id: str
) -> None:
    """An overwrite nobody asked for skips the submission whole, NULLs it could have filled included.

    A row the run updated in part is harder to reason about later than one it left alone, and the
    next run fills those NULLs once the overwrite is settled.
    """
    current = db.add_submission(submission_id)
    current.submission_size = 1  # non-NULL value that will differ from metadata
    db.update_submission(current)
    current = db.get_submission(submission_id)
    assert current.submission_metadata is None, "the fixture must leave something additive to write"

    result = _backfill_submission(
        current_submission=current,
        raw_json=_metadata_json(metadata),
        db_service=db,
        dry_run=False,
        force=False,
        ignore_fields=set(),
    )

    assert result.status == _BackfillResult.WOULD_OVERWRITE
    persisted = db.get_submission(submission_id)
    assert persisted.submission_size == 1, "an existing value must not be overwritten"
    assert persisted.submission_metadata is None, "and the rest of the submission waits with it"


def test_backfill_submission_returns_would_overwrite_when_only_overwrites_are_pending(
    db: SubmissionDb, metadata: GrzSubmissionMetadata, submission_id: str
) -> None:
    """When every pending change is an overwrite and none is permitted, nothing is written."""
    db.add_submission(submission_id)

    # bring the row up to date first, so the only difference afterwards is the one below
    _backfill_submission(
        current_submission=db.get_submission(submission_id),
        raw_json=_metadata_json(metadata),
        db_service=db,
        dry_run=False,
        force=True,
        ignore_fields=set(),
    )
    current = db.get_submission(submission_id)
    current.submission_size = 1
    db.update_submission(current)

    result = _backfill_submission(
        current_submission=db.get_submission(submission_id),
        raw_json=_metadata_json(metadata),
        db_service=db,
        dry_run=False,
        force=False,
        ignore_fields=set(),
    )

    assert result.status == _BackfillResult.WOULD_OVERWRITE
    assert db.get_submission(submission_id).submission_size == 1


def test_backfill_submission_force_applies_destructive_changes(
    db: SubmissionDb, metadata: GrzSubmissionMetadata, submission_id: str
) -> None:
    """With --force, a submission whose non-NULL field differs from S3 is updated even though it overwrites data."""
    current = db.add_submission(submission_id)
    current.submission_size = 1  # non-NULL value that will differ from metadata
    db.update_submission(current)
    current = db.get_submission(submission_id)

    result = _backfill_submission(
        current_submission=current,
        raw_json=_metadata_json(metadata),
        db_service=db,
        dry_run=False,
        force=True,
        ignore_fields=set(),
    )

    assert result.status == _BackfillResult.UPDATED
    persisted = db.get_submission(submission_id)
    assert persisted.submission_size == metadata.get_submission_size()
    assert persisted.submission_metadata == metadata.to_redacted_dict()


def test_backfill_force_reconciles_against_an_unredacted_copy(
    db: SubmissionDb, metadata: GrzSubmissionMetadata, submission_id: str
) -> None:
    """An unredacted copy is authoritative, so --force reconciles the database against it.

    The upload date is the exception: metadata.json does not carry one, so diff() excludes it.
    """
    current = db.add_submission(submission_id)
    current.submission_uploaded_date = DIFFERENT_DATE
    current.tan_g = DIFFERENT_TAN_G
    current.local_case_id = DIFFERENT_LOCAL_CASE_ID
    db.update_submission(current)
    current = db.get_submission(submission_id)
    assert metadata.submission.tan_g != DIFFERENT_TAN_G
    assert metadata.submission.local_case_id != DIFFERENT_LOCAL_CASE_ID
    assert metadata.submission.submission_date != DIFFERENT_DATE

    result = _backfill_submission(
        current_submission=current,
        raw_json=_metadata_json(metadata),
        db_service=db,
        dry_run=False,
        force=True,
        ignore_fields=set(),
    )

    assert result.status == _BackfillResult.UPDATED
    persisted = db.get_submission(submission_id)
    assert persisted.submission_uploaded_date == DIFFERENT_DATE
    assert persisted.tan_g == metadata.submission.tan_g
    assert persisted.local_case_id == metadata.submission.local_case_id
    assert persisted.submission_size == metadata.get_submission_size()
    assert persisted.submission_metadata is not None
    assert persisted.submission_metadata == metadata.to_redacted_dict()


def test_backfill_leaves_a_missing_upload_date_alone(
    db: SubmissionDb, metadata: GrzSubmissionMetadata, submission_id: str
) -> None:
    """A NULL upload date stays NULL, even with --force.

    The column records when the upload finished; the metadata's submission_date is what the
    submitter declared, and it drives the reporting windows, so a declared date must not stand in
    for a real one. diff() makes that call: given no date it falls back to the metadata's and then
    ignores the field, so backfill passes the stored value straight through rather than deciding
    again.
    """
    current = db.add_submission(submission_id)
    assert current.submission_uploaded_date is None
    assert metadata.submission.submission_date is not None

    result = _backfill_submission(
        current_submission=current,
        raw_json=_metadata_json(metadata),
        db_service=db,
        dry_run=False,
        force=True,
        ignore_fields=set(),
    )

    assert result.status == _BackfillResult.UPDATED
    assert db.get_submission(submission_id).submission_uploaded_date is None


def test_backfill_does_not_write_an_upload_date_from_a_stale_snapshot(
    db: SubmissionDb, metadata: GrzSubmissionMetadata, submission_id: str
) -> None:
    """A row corrected mid-run keeps the correction, not the value the run started with.

    Backfill builds its candidate list once and then spends the run fetching from S3, so the
    Submission it holds is a snapshot: diff() re-reads the row from the database, and the two can
    disagree if an operator writes in between. Offering the snapshot's date would write it back over
    the correction, and over a deliberate clear it would not even need --force, since filling a NULL
    is additive.
    """
    db.add_submission(submission_id)
    db.modify_submission(submission_id, "submission_uploaded_date", DIFFERENT_DATE)
    stale = db.get_submission(submission_id)
    assert stale.submission_uploaded_date == DIFFERENT_DATE

    # the operator clears the column while the run is in flight
    db.modify_submission(submission_id, "submission_uploaded_date", None)

    result = _backfill_submission(
        current_submission=stale,
        raw_json=_metadata_json(metadata),
        db_service=db,
        dry_run=False,
        force=True,
        ignore_fields=set(),
    )

    assert result.status == _BackfillResult.UPDATED
    assert db.get_submission(submission_id).submission_uploaded_date is None


def test_backfill_submission_reads_a_consent_datetime_without_a_timezone(
    db: SubmissionDb, metadata: GrzSubmissionMetadata, submission_id: str
) -> None:
    """A submission accepted before FHIR's timezone rule was enforced is read, not skipped.

    The models keep a consent dateTime as submitted, so such a document parses as it stands and
    backfill needs no repair step. What it stores is what the submitter wrote.
    """
    current = db.add_submission(submission_id)
    raw = json.loads(metadata.model_dump_json(by_alias=True))
    scope = raw["donors"][0]["researchConsents"][0]["scope"]
    scope["dateTime"] = "2020-09-01T14:37:22"  # as an older submission would have stated it

    result = _backfill_submission(
        current_submission=current,
        raw_json=json.dumps(raw),
        db_service=db,
        dry_run=False,
        force=False,
        ignore_fields=set(),
    )

    assert result.status == _BackfillResult.UPDATED
    stored = db.get_submission(submission_id).submission_metadata
    stored_scope = stored["donors"][0]["researchConsents"][0]["scope"]
    assert stored_scope["dateTime"] == "2020-09-01T14:37:22", "the submitted value must be stored as written"


def test_backfill_submission_allow_overwrite_writes_the_field_it_names(
    db: SubmissionDb, metadata: GrzSubmissionMetadata, submission_id: str
) -> None:
    """Once --allow-overwrite covers the overwrite, the submission is written whole.

    Refreshing stored submission_metadata is the motivating case: it has to be rewritten from S3
    without --force also permitting every other overwrite the same diff carries.
    """
    current = db.add_submission(submission_id)
    current.submission_metadata = {"stale": True}  # non-NULL, differs, and IS allowed to change
    db.update_submission(current)
    current = db.get_submission(submission_id)

    result = _backfill_submission(
        current_submission=current,
        raw_json=_metadata_json(metadata),
        db_service=db,
        dry_run=False,
        force=False,
        ignore_fields=set(),
        allow_overwrite=frozenset({"submission_metadata"}),
    )

    assert result.status == _BackfillResult.UPDATED
    persisted = db.get_submission(submission_id)
    assert persisted.submission_metadata == metadata.to_redacted_dict(), "the named field must be refreshed"
    assert persisted.submission_size == metadata.get_submission_size(), "the NULLs are filled with it"


def test_backfill_submission_allow_overwrite_reports_would_overwrite_when_nothing_is_writable(
    db: SubmissionDb, metadata: GrzSubmissionMetadata, submission_id: str
) -> None:
    """When every pending change is an overwrite the allow-list does not cover, nothing is written."""
    db.add_submission(submission_id)

    # Bring the row fully up to date first, so the only pending change afterwards is the one below.
    # A freshly added submission still has many NULL columns, and filling those is additive.
    _backfill_submission(
        current_submission=db.get_submission(submission_id),
        raw_json=_metadata_json(metadata),
        db_service=db,
        dry_run=False,
        force=True,
        ignore_fields=set(),
    )

    current = db.get_submission(submission_id)
    current.submission_size = 1  # now the only difference, and not in the allow-list
    db.update_submission(current)

    result = _backfill_submission(
        current_submission=db.get_submission(submission_id),
        raw_json=_metadata_json(metadata),
        db_service=db,
        dry_run=False,
        force=False,
        ignore_fields=set(),
        allow_overwrite=frozenset({"submission_metadata"}),
    )

    assert result.status == _BackfillResult.WOULD_OVERWRITE
    assert db.get_submission(submission_id).submission_size == 1


def test_backfill_never_overwrites_stored_values_with_placeholders(
    db: SubmissionDb, metadata: GrzSubmissionMetadata, submission_id: str
) -> None:
    """A redacted archive copy is restored from the row first, so the stored values survive."""
    db.add_submission(submission_id)
    db.modify_submission(submission_id, "tan_g", DIFFERENT_TAN_G)
    db.modify_submission(submission_id, "local_case_id", DIFFERENT_LOCAL_CASE_ID)
    current = db.get_submission(submission_id)

    assert _run_backfill(db, _archived(metadata), current).status == _BackfillResult.UPDATED

    persisted = db.get_submission(submission_id)
    assert persisted.tan_g == DIFFERENT_TAN_G
    assert persisted.local_case_id == DIFFERENT_LOCAL_CASE_ID
    # and the restored local case ID is what keyed the case, not the placeholder
    assert [case.local_case_id for case, _count in db.list_cases()] == [DIFFERENT_LOCAL_CASE_ID]


def _as_initial(metadata: GrzSubmissionMetadata) -> GrzSubmissionMetadata:
    raw = json.loads(metadata.model_dump_json(by_alias=True))
    raw["submission"]["submissionType"] = "initial"
    return GrzSubmissionMetadata.model_validate(raw)


def test_backfill_submission_links_case_by_default(db: SubmissionDb, metadata: GrzSubmissionMetadata) -> None:
    """Backfill links submissions to cases by default (skip via --ignore-field case_id)."""
    initial_metadata = _as_initial(metadata)
    sid = initial_metadata.submission_id
    db.add_submission(sid)
    db.commit_changes(sid, db.diff(sid, initial_metadata, submission_uploaded_date=None, ignore_fields={"case_id"}))
    current = db.get_submission(sid)
    assert current.case_id is None

    result = _backfill_submission(
        current_submission=current,
        raw_json=_metadata_json(initial_metadata),
        db_service=db,
        dry_run=False,
        force=False,
        ignore_fields={"case_id"},
    )
    assert result.status == _BackfillResult.UP_TO_DATE
    assert db.get_submission(sid).case_id is None

    result = _backfill_submission(
        current_submission=db.get_submission(sid),
        raw_json=_metadata_json(initial_metadata),
        db_service=db,
        dry_run=False,
        force=False,
        ignore_fields=set(),
    )
    assert result.status == _BackfillResult.UPDATED
    linked = db.get_submission(sid)
    assert linked.case_id is not None
    # this copy is unredacted, so it is authoritative and fills the NULL local case ID too
    assert linked.local_case_id == initial_metadata.submission.local_case_id


def _archived(metadata: GrzSubmissionMetadata) -> GrzSubmissionMetadata:
    """The copy archival uploads, in its older spelling: tanG zeroed and localCaseId emptied.

    The archive still holds objects written before ``REDACTED_LOCAL_CASE_ID`` replaced this
    spelling. Their case key is a placeholder shared by every submission of that submitter.
    """
    raw = metadata.to_redacted_dict()
    raw["submission"]["submissionType"] = "initial"
    # the older archival spelling, which is still what most objects in the archive carry
    raw["submission"]["localCaseId"] = ""
    return GrzSubmissionMetadata.model_validate(raw)


def _run_backfill(db: SubmissionDb, metadata: GrzSubmissionMetadata, current: Submission) -> _BackfillOutcome:
    return _backfill_submission(
        current_submission=current,
        raw_json=_metadata_json(metadata),
        db_service=db,
        dry_run=False,
        force=False,
        ignore_fields=set(),
    )


def test_backfill_keys_cases_on_the_stored_local_case_id(db: SubmissionDb, metadata: GrzSubmissionMetadata) -> None:
    """Two patients whose archived metadata both read localCaseId "" must not share a case."""
    archived = _archived(metadata)
    submitter = metadata.submission.submitter_id
    rows = []
    for sid, tan_g, local_case_id in (
        (f"{submitter}_2024-01-01_aaaaaaa1", "a" * 64, "patient-A"),
        (f"{submitter}_2024-01-02_aaaaaaa2", "b" * 64, "patient-B"),
    ):
        db.add_submission(sid)
        db.modify_submission(sid, "tan_g", tan_g)
        db.modify_submission(sid, "local_case_id", local_case_id)
        db.modify_submission(sid, "submission_type", "initial")
        rows.append(db.get_submission(sid))

    for row in rows:
        assert _run_backfill(db, archived, row).status == _BackfillResult.UPDATED

    assert {(case.submitter_id, case.local_case_id) for case, _count in db.list_cases()} == {
        (submitter, "patient-A"),
        (submitter, "patient-B"),
    }
    linked = {row.id: db.get_submission(row.id).case_id for row in rows}
    assert None not in linked.values()
    assert len(set(linked.values())) == 2


def test_backfill_without_a_stored_local_case_id_skips_the_case_link(
    db: SubmissionDb, metadata: GrzSubmissionMetadata, submission_id: str
) -> None:
    """Nothing to restore from, so the placeholders are ignored rather than written or keyed on."""
    current = db.add_submission(submission_id)

    assert _run_backfill(db, _archived(metadata), current).status == _BackfillResult.UPDATED

    persisted = db.get_submission(submission_id)
    assert persisted.submission_size == metadata.get_submission_size()
    assert persisted.case_id is None
    assert persisted.local_case_id is None
    assert persisted.tan_g is None
    assert db.list_cases() == []


def _second_case_for_the_same_key(db: SubmissionDb, submitter_id: str, local_case_id: str) -> None:
    """Give one key a second case, which ``ux_cases_submitter_local_case`` forbids.

    Only reachable by writing around the application, which is exactly the state backfill's
    unresolvable-link handling exists to survive: the index is dropped first so the row can
    be inserted at all.
    """
    with db.transaction() as session:
        session.execute(sqlalchemy.text("DROP INDEX ux_cases_submitter_local_case"))
        session.execute(
            sqlalchemy.text("INSERT INTO cases (submitter_id, local_case_id) VALUES (:submitter, :local_case)"),
            {"submitter": submitter_id, "local_case": local_case_id},
        )
        session.commit()


def test_backfill_writes_everything_but_the_link_when_the_case_key_is_ambiguous(
    db: SubmissionDb, metadata: GrzSubmissionMetadata
) -> None:
    """An ambiguous key needs an operator to merge the cases; the submission still gets recorded.

    Discarding the change set would leave submission_size, submission_metadata and the donors
    unwritten, which are what the Prüfbericht is built from.
    """
    submitter = metadata.submission.submitter_id
    db.create_case(submitter, "patient-A")
    _second_case_for_the_same_key(db, submitter, "patient-A")

    sid = f"{submitter}_2024-01-01_aaaaaaa1"
    db.add_submission(sid)
    db.modify_submission(sid, "local_case_id", "patient-A")
    current = db.get_submission(sid)

    # link_unresolved=True: the key is ambiguous, not merely missing. See _BackfillOutcome.
    assert _run_backfill(db, _archived(metadata), current) == _BackfillOutcome(_BackfillResult.UPDATED, True)

    persisted = db.get_submission(sid)
    assert persisted.case_id is None
    assert persisted.submission_size == metadata.get_submission_size()
    assert persisted.submission_metadata is not None
    assert db.get_donors(sid)


def test_backfill_links_once_the_ambiguity_is_gone(db: SubmissionDb, metadata: GrzSubmissionMetadata) -> None:
    """Re-running after an operator deletes the spurious duplicate case completes the link."""
    submitter = metadata.submission.submitter_id
    kept = db.create_case(submitter, "patient-A")
    _second_case_for_the_same_key(db, submitter, "patient-A")
    spare = next(case for case, _n in db.list_cases() if case.id != kept.id)

    sid = f"{submitter}_2024-01-01_aaaaaaa1"
    db.add_submission(sid)
    db.modify_submission(sid, "local_case_id", "patient-A")
    assert _run_backfill(db, _archived(metadata), db.get_submission(sid)) == _BackfillOutcome(
        _BackfillResult.UPDATED, True
    )

    db.delete_case(spare.id)

    assert _run_backfill(db, _archived(metadata), db.get_submission(sid)) == _BackfillOutcome(_BackfillResult.UPDATED)
    assert db.get_submission(sid).case_id is not None


def _flip_a_donors_mv_consent(db: SubmissionDb, submission_id: str) -> tuple[str, bool]:
    """Make one stored donor differ from metadata.json, and return its pseudonym and the stored value."""
    donor = db.get_donors(submission_id)[0]
    donor.mv_consented = not donor.mv_consented
    db.update_donor(donor)
    return donor.pseudonym, donor.mv_consented


def test_backfill_skips_the_submission_when_a_donor_changed_without_force(
    db: SubmissionDb, metadata: GrzSubmissionMetadata, submission_id: str
) -> None:
    """A stored donor that differs from metadata.json skips the submission, like any other overwrite."""
    _populate_full_row(db, submission_id, metadata)
    pseudonym, stored = _flip_a_donors_mv_consent(db, submission_id)
    db.modify_submission(submission_id, "submission_size", None)  # something additive that must wait

    result = _backfill_submission(
        current_submission=db.get_submission(submission_id),
        raw_json=_metadata_json(metadata),
        db_service=db,
        dry_run=False,
        force=False,
        ignore_fields=set(),
    )

    assert result.status == _BackfillResult.WOULD_OVERWRITE
    assert db.get_submission(submission_id).submission_size is None
    assert db.get_donors(submission_id, pseudonym)[0].mv_consented == stored, "the donor must not be overwritten"


def test_backfill_allow_overwrite_donors_writes_a_changed_donor(
    db: SubmissionDb, metadata: GrzSubmissionMetadata, submission_id: str
) -> None:
    """``donors`` is how an operator permits donor overwrites without permitting every other one."""
    _populate_full_row(db, submission_id, metadata)
    pseudonym, stored = _flip_a_donors_mv_consent(db, submission_id)

    result = _backfill_submission(
        current_submission=db.get_submission(submission_id),
        raw_json=_metadata_json(metadata),
        db_service=db,
        dry_run=False,
        force=False,
        ignore_fields=set(),
        allow_overwrite=frozenset({DONORS_KEY}),
    )

    assert result.status == _BackfillResult.UPDATED
    assert db.get_donors(submission_id, pseudonym)[0].mv_consented == (not stored)


def test_backfill_force_overwrites_a_changed_donor(
    db: SubmissionDb, metadata: GrzSubmissionMetadata, submission_id: str
) -> None:
    _populate_full_row(db, submission_id, metadata)
    pseudonym, stored = _flip_a_donors_mv_consent(db, submission_id)

    result = _backfill_submission(
        current_submission=db.get_submission(submission_id),
        raw_json=_metadata_json(metadata),
        db_service=db,
        dry_run=False,
        force=True,
        ignore_fields=set(),
    )

    assert result.status == _BackfillResult.UPDATED
    assert db.get_donors(submission_id, pseudonym)[0].mv_consented == (not stored)


def test_backfill_holds_back_the_whole_submission_when_the_case_link_changed(
    db: SubmissionDb, metadata: GrzSubmissionMetadata
) -> None:
    """A relinked submission keeps its case, and nothing else is written either.

    Replacing the link would undo a deliberate relink, and --allow-overwrite cannot name it, so
    only --force gets past this one.
    """
    initial_metadata = _as_initial(metadata)
    sid = initial_metadata.submission_id
    _populate_full_row(db, sid, initial_metadata)
    other = db.create_case(initial_metadata.submission.submitter_id, "some-other-case")
    db.set_submission_case(sid, other.id)
    db.modify_submission(sid, "submission_size", None)  # something additive that must still wait

    result = _backfill_submission(
        current_submission=db.get_submission(sid),
        raw_json=_metadata_json(initial_metadata),
        db_service=db,
        dry_run=False,
        force=False,
        ignore_fields=set(),
    )

    assert result.status == _BackfillResult.WOULD_OVERWRITE
    persisted = db.get_submission(sid)
    assert persisted.case_id == other.id
    assert persisted.submission_size is None


def _invoke_backfill_command(config_path: Path, submission_id: str) -> click.testing.Result:
    runner = click.testing.CliRunner()
    return runner.invoke(
        grzctl.cli.build_cli(),
        ["--config", str(config_path), "db", "backfill", "--submission-id", submission_id],
    )


def test_fetch_from_archives_finds_nothing_when_no_archive_holds_the_metadata(
    s3_client_mock: Any, submission_id: str
) -> None:
    for bucket in ARCHIVE_BUCKETS:
        s3_client_mock.create_bucket(Bucket=bucket)
    archive_targets = [(bucket, bucket, s3_client_mock) for bucket in ARCHIVE_BUCKETS]

    assert _fetch_metadata_json_from_archives(submission_id, archive_targets) == {}


def test_backfill_writes_the_metadata_from_the_one_archive_that_holds_it(
    db: SubmissionDb,
    s3_client_mock: Any,
    migrated_database_config_path: Path,
    metadata: GrzSubmissionMetadata,
    submission_id: str,
) -> None:
    db.add_submission(submission_id)
    for bucket in ARCHIVE_BUCKETS:
        s3_client_mock.create_bucket(Bucket=bucket)
    _put_metadata(s3_client_mock, "non_consented", submission_id, metadata)

    result = _invoke_backfill_command(migrated_database_config_path, submission_id)

    assert result.exit_code == 0, result.stderr
    assert db.get_submission(submission_id).submission_metadata == metadata.to_redacted_dict()


def test_backfill_writes_nothing_when_both_archives_hold_the_metadata(
    db: SubmissionDb,
    s3_client_mock: Any,
    migrated_database_config_path: Path,
    metadata: GrzSubmissionMetadata,
    submission_id: str,
) -> None:
    """A metadata.json in both archives is an error, and neither copy reaches the database.

    A copy committed before the conflict is noticed cannot be repaired by a rerun, because the
    rerun finds both copies again.
    """
    db.add_submission(submission_id)
    for bucket in ARCHIVE_BUCKETS:
        s3_client_mock.create_bucket(Bucket=bucket)
        _put_metadata(s3_client_mock, bucket, submission_id, metadata)

    result = _invoke_backfill_command(migrated_database_config_path, submission_id)

    assert result.exit_code == 1, "any error fails the run"
    assert "found in both" in result.stderr
    assert db.get_submission(submission_id).submission_metadata is None


def test_backfill_writes_nothing_when_an_archive_cannot_be_read(
    db: SubmissionDb,
    s3_client_mock: Any,
    migrated_database_config_path: Path,
    metadata: GrzSubmissionMetadata,
    submission_id: str,
) -> None:
    """An archive that cannot be read might hold a second copy, so the copy found in the other one is not written."""
    db.add_submission(submission_id)
    s3_client_mock.create_bucket(Bucket="non_consented")  # no consented bucket, so reading it fails
    _put_metadata(s3_client_mock, "non_consented", submission_id, metadata)

    result = _invoke_backfill_command(migrated_database_config_path, submission_id)

    assert result.exit_code == 1, "any error fails the run"
    assert "S3 error in consented archive" in result.stderr
    assert db.get_submission(submission_id).submission_metadata is None
