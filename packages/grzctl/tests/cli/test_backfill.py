# Backfill unit tests, see .vbw-planning/phases/02-implement-grzctl-db-backfill-command/02-02-PLAN.md
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
import sqlalchemy
from grz_db.models.submission import Submission, SubmissionBase, SubmissionDb
from grz_pydantic_models.submission.metadata import GrzSubmissionMetadata
from grz_pydantic_models_testing.example_metadata import grzctl as grzctl_metadata
from grzctl.commands.db.cli import _BACKFILL_IGNORE_FIELDS, _backfill_submission, _BackfillResult
from moto import mock_aws

BUCKET = "test-backfill-bucket"
REGION = "us-east-1"
IGNORE_FIELDS = _BACKFILL_IGNORE_FIELDS
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
    """Write *metadata* unredacted, unlike archival, which redacts first.

    For the redacted shape backfill actually reads from the archive, use ``_archived`` instead.
    """
    s3_client.put_object(
        Bucket=BUCKET,
        Key=f"{submission_id}/metadata/metadata.json",
        Body=json.dumps(metadata.get_raw_dict()).encode("utf-8"),
    )


def _populate_full_row(db: SubmissionDb, submission_id: str, metadata: GrzSubmissionMetadata) -> Submission:
    """Persist a fully-populated row by running the same diff/commit path the production code uses."""
    db.add_submission(submission_id)
    db.commit_changes(submission_id, db.diff(submission_id, metadata, submission_uploaded_date=None))
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


def test_backfill_force_reconciles_against_an_unredacted_copy(
    db: SubmissionDb, s3_client_mock: Any, metadata: GrzSubmissionMetadata, submission_id: str
) -> None:
    """An unredacted copy is authoritative, so --force reconciles the database against it.

    Only the upload date stays hard-ignored, since metadata.json does not carry one.
    """
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
    assert persisted.tan_g == metadata.submission.tan_g
    assert persisted.pseudonym == metadata.submission.local_case_id
    assert persisted.submission_size == metadata.get_submission_size()
    assert persisted.submission_metadata is not None
    assert persisted.submission_metadata == metadata.to_redacted_dict()


def test_backfill_leaves_a_missing_upload_date_alone(
    db: SubmissionDb, s3_client_mock: Any, metadata: GrzSubmissionMetadata, submission_id: str
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
    assert db.get_submission(submission_id).submission_uploaded_date is None


def test_backfill_does_not_write_an_upload_date_from_a_stale_snapshot(
    db: SubmissionDb, s3_client_mock: Any, metadata: GrzSubmissionMetadata, submission_id: str
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
    _put_metadata(s3_client_mock, submission_id, metadata)

    result = _backfill_submission(
        current_submission=stale,
        s3_client=s3_client_mock,
        bucket=BUCKET,
        db_service=db,
        dry_run=False,
        force=True,
        ignore_fields=IGNORE_FIELDS,
    )

    assert result == _BackfillResult.UPDATED
    assert db.get_submission(submission_id).submission_uploaded_date is None


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


def test_backfill_never_overwrites_stored_values_with_placeholders(
    db: SubmissionDb, s3_client_mock: Any, metadata: GrzSubmissionMetadata, submission_id: str
) -> None:
    """A redacted archive copy is restored from the row first, so the stored values survive."""
    db.add_submission(submission_id)
    db.modify_submission(submission_id, "tan_g", DIFFERENT_TAN_G)
    db.modify_submission(submission_id, "pseudonym", DIFFERENT_PSEUDONYM)
    current = db.get_submission(submission_id)
    _put_metadata(s3_client_mock, submission_id, _archived(metadata))

    assert _run_backfill(db, s3_client_mock, current) == _BackfillResult.UPDATED

    persisted = db.get_submission(submission_id)
    assert persisted.tan_g == DIFFERENT_TAN_G
    assert persisted.pseudonym == DIFFERENT_PSEUDONYM
    # and the restored pseudonym is what keyed the case, not the placeholder
    assert [case.local_case_id for case, _count in db.list_cases()] == [DIFFERENT_PSEUDONYM]


def _as_initial(metadata: GrzSubmissionMetadata) -> GrzSubmissionMetadata:
    raw = json.loads(metadata.model_dump_json(by_alias=True))
    raw["submission"]["submissionType"] = "initial"
    return GrzSubmissionMetadata.model_validate(raw)


def test_backfill_submission_links_case_by_default(
    db: SubmissionDb, s3_client_mock: Any, metadata: GrzSubmissionMetadata
) -> None:
    """Backfill links submissions to cases by default (skip via --ignore-field case_id)."""
    initial_metadata = _as_initial(metadata)
    sid = initial_metadata.submission_id
    db.add_submission(sid)
    db.commit_changes(sid, db.diff(sid, initial_metadata, submission_uploaded_date=None, ignore_fields={"case_id"}))
    current = db.get_submission(sid)
    assert current.case_id is None
    _put_metadata(s3_client_mock, sid, initial_metadata)

    result = _backfill_submission(
        current_submission=current,
        s3_client=s3_client_mock,
        bucket=BUCKET,
        db_service=db,
        dry_run=False,
        force=False,
        ignore_fields=IGNORE_FIELDS | {"case_id"},
    )
    assert result == _BackfillResult.UP_TO_DATE
    assert db.get_submission(sid).case_id is None

    result = _backfill_submission(
        current_submission=db.get_submission(sid),
        s3_client=s3_client_mock,
        bucket=BUCKET,
        db_service=db,
        dry_run=False,
        force=False,
        ignore_fields=IGNORE_FIELDS,
    )
    assert result == _BackfillResult.UPDATED
    linked = db.get_submission(sid)
    assert linked.case_id is not None
    # this copy is unredacted, so it is authoritative and fills the NULL pseudonym too
    assert linked.pseudonym == initial_metadata.submission.local_case_id


def _archived(metadata: GrzSubmissionMetadata) -> GrzSubmissionMetadata:
    """The copy archival uploads, in its older spelling: tanG zeroed and localCaseId emptied.

    Mirrors ``S3BotoUploadWorker.archive``, so the case key in an archive bucket is a
    placeholder shared by every submission of that submitter.
    """
    raw = metadata.to_redacted_dict()
    raw["submission"]["submissionType"] = "initial"
    # the older archival spelling, which is still what most objects in the archive carry
    raw["submission"]["localCaseId"] = ""
    return GrzSubmissionMetadata.model_validate(raw)


def _run_backfill(db: SubmissionDb, s3_client: Any, current: Submission) -> _BackfillResult:
    return _backfill_submission(
        current_submission=current,
        s3_client=s3_client,
        bucket=BUCKET,
        db_service=db,
        dry_run=False,
        force=False,
        ignore_fields=IGNORE_FIELDS,
    )


def test_backfill_keys_cases_on_the_stored_pseudonym(
    db: SubmissionDb, s3_client_mock: Any, metadata: GrzSubmissionMetadata
) -> None:
    """Two patients whose archived metadata both read localCaseId "" must not share a case."""
    archived = _archived(metadata)
    submitter = metadata.submission.submitter_id
    rows = []
    for sid, tan_g, pseudonym in (
        (f"{submitter}_2024-01-01_aaaaaaa1", "a" * 64, "patient-A"),
        (f"{submitter}_2024-01-02_aaaaaaa2", "b" * 64, "patient-B"),
    ):
        db.add_submission(sid)
        db.modify_submission(sid, "tan_g", tan_g)
        db.modify_submission(sid, "pseudonym", pseudonym)
        db.modify_submission(sid, "submission_type", "initial")
        _put_metadata(s3_client_mock, sid, archived)
        rows.append(db.get_submission(sid))

    for row in rows:
        assert _run_backfill(db, s3_client_mock, row) == _BackfillResult.UPDATED

    assert {(case.submitter_id, case.local_case_id) for case, _count in db.list_cases()} == {
        (submitter, "patient-A"),
        (submitter, "patient-B"),
    }
    linked = {row.id: db.get_submission(row.id).case_id for row in rows}
    assert None not in linked.values()
    assert len(set(linked.values())) == 2


def test_backfill_without_a_stored_pseudonym_skips_the_case_link(
    db: SubmissionDb, s3_client_mock: Any, metadata: GrzSubmissionMetadata, submission_id: str
) -> None:
    """Nothing to restore from, so the placeholders are ignored rather than written or keyed on."""
    current = db.add_submission(submission_id)
    _put_metadata(s3_client_mock, submission_id, _archived(metadata))

    assert _run_backfill(db, s3_client_mock, current) == _BackfillResult.UPDATED

    persisted = db.get_submission(submission_id)
    assert persisted.submission_size == metadata.get_submission_size()
    assert persisted.case_id is None
    assert persisted.pseudonym is None
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
    db: SubmissionDb, s3_client_mock: Any, metadata: GrzSubmissionMetadata
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
    db.modify_submission(sid, "pseudonym", "patient-A")
    current = db.get_submission(sid)
    _put_metadata(s3_client_mock, sid, _archived(metadata))

    assert _run_backfill(db, s3_client_mock, current) == _BackfillResult.LINK_UNRESOLVED

    persisted = db.get_submission(sid)
    assert persisted.case_id is None
    assert persisted.submission_size == metadata.get_submission_size()
    assert persisted.submission_metadata is not None
    assert db.get_donors(sid)


def test_backfill_links_once_the_ambiguity_is_gone(
    db: SubmissionDb, s3_client_mock: Any, metadata: GrzSubmissionMetadata
) -> None:
    """Re-running after an operator deletes the spurious duplicate case completes the link."""
    submitter = metadata.submission.submitter_id
    kept = db.create_case(submitter, "patient-A")
    _second_case_for_the_same_key(db, submitter, "patient-A")
    spare = next(case for case, _n in db.list_cases() if case.id != kept.id)

    sid = f"{submitter}_2024-01-01_aaaaaaa1"
    db.add_submission(sid)
    db.modify_submission(sid, "pseudonym", "patient-A")
    _put_metadata(s3_client_mock, sid, _archived(metadata))
    assert _run_backfill(db, s3_client_mock, db.get_submission(sid)) == _BackfillResult.LINK_UNRESOLVED

    db.delete_case(spare.id)

    assert _run_backfill(db, s3_client_mock, db.get_submission(sid)) == _BackfillResult.UPDATED
    assert db.get_submission(sid).case_id is not None


def test_backfill_ignores_only_fields_that_exist() -> None:
    """Every name in ``_BACKFILL_IGNORE_FIELDS`` must be an actual Submission column.

    A name that matches no column silently ignores nothing (set difference against a
    non-member is a no-op), so a typo here would only surface once it let backfill overwrite
    a field it was meant to skip; this asserts it directly instead.
    """
    assert _BACKFILL_IGNORE_FIELDS <= SubmissionBase.model_fields.keys()
