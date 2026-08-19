import datetime
import json
from collections.abc import Callable

import pytest
from grz_db.errors import DuplicateTanGError
from grz_db.models.submission import OutdatedDatabaseSchemaError, Submission, SubmissionDb
from grz_pydantic_models.submission.metadata import (
    REDACTED_LOCAL_CASE_ID,
    REDACTED_TAN,
    GrzSubmissionMetadata,
)

TWO_TB = 2 * 1024**4  # 2,199,023,255,552 bytes
SUBMISSION_ID = "123456789_2024-01-01_abcdef01"
SUBMISSION_ID_2 = "123456789_2024-01-02_abcdef02"
SUBMISSION_ID_3 = "123456789_2024-01-03_abcdef03"
TAN_G_1 = "a" * 64
TAN_G_2 = "b" * 64


@pytest.fixture(scope="function")
def submission(db: SubmissionDb):
    """A freshly added submission, shared across tests in this module."""
    return db.add_submission(SUBMISSION_ID)


def test_submission_size_can_store_2tb(db: SubmissionDb, submission) -> None:
    """submission_size uses BigInteger and must be able to store values >= 2 TB."""
    submission.submission_size = TWO_TB
    db.update_submission(submission)

    result = db.get_submission(SUBMISSION_ID)
    assert result is not None
    assert result.submission_size == TWO_TB


def test_submission_metadata_json_roundtrip(db: SubmissionDb, submission) -> None:
    """submission_metadata (JSON/JSONB) must survive a write-read roundtrip with nested structures."""
    metadata = {
        "string_field": "hello",
        "int_field": 42,
        "float_field": 3.14,
        "bool_field": True,
        "null_field": None,
        "list_field": [1, "two", 3.0],
        "nested": {"a": 1, "b": [2, 3]},
    }

    submission.submission_metadata = metadata
    db.update_submission(submission)

    result = db.get_submission(SUBMISSION_ID)
    assert result is not None
    assert result.submission_metadata == metadata


def test_from_metadata_sets_fields_from_metadata(metadata: GrzSubmissionMetadata) -> None:
    """Submission.from_metadata must map every metadata field correctly and leave system fields unset."""
    explicit_date = datetime.date(2025, 3, 1)
    submission = Submission.from_metadata(SUBMISSION_ID, metadata, explicit_date)

    # --- fields that must be populated from the metadata object ---
    assert submission.id == SUBMISSION_ID
    assert submission.tan_g == metadata.submission.tan_g
    assert submission.submission_type == metadata.submission.submission_type
    assert submission.submitter_id == metadata.submission.submitter_id
    assert submission.coverage_type == metadata.submission.coverage_type
    assert submission.disease_type == metadata.submission.disease_type
    assert submission.genomic_study_type == metadata.submission.genomic_study_type
    assert submission.genomic_study_subtype == metadata.submission.genomic_study_subtype
    assert submission.pseudonym == metadata.submission.local_case_id
    assert submission.data_node_id == metadata.submission.genomic_data_center_id
    assert submission.submission_uploaded_date == explicit_date  # explicit date takes precedence
    assert submission.submission_size == metadata.get_submission_size()
    assert submission.submission_metadata == metadata.to_redacted_dict()

    # --- system-managed fields must not be in model_fields_set ---
    system_fields = {"basic_qc_passed", "detailed_qc_passed", "selected_for_qc"}
    assert system_fields.isdisjoint(submission.model_fields_set), (
        f"System fields unexpectedly set: {system_fields & submission.model_fields_set}"
    )


def test_get_submissions_returns_matching(db: SubmissionDb) -> None:
    """get_submissions returns one entry per requested ID in the same order."""
    db.add_submission(SUBMISSION_ID)
    db.add_submission(SUBMISSION_ID_2)
    db.add_submission(SUBMISSION_ID_3)

    result = db.get_submissions([SUBMISSION_ID, SUBMISSION_ID_3])

    assert len(result) == 2
    assert result[0] is not None and result[0].id == SUBMISSION_ID
    assert result[1] is not None and result[1].id == SUBMISSION_ID_3


def test_get_submissions_empty_input(db: SubmissionDb) -> None:
    """get_submissions with an empty list returns an empty list."""
    assert db.get_submissions([]) == []


def test_get_submissions_unknown_ids_map_to_none(db: SubmissionDb) -> None:
    """get_submissions maps unknown IDs to None, preserving position."""
    db.add_submission(SUBMISSION_ID)
    missing = "000000000_2000-01-01_deadbeef"

    result = db.get_submissions([SUBMISSION_ID, missing])

    assert len(result) == 2
    assert result[0] is not None and result[0].id == SUBMISSION_ID
    assert result[1] is None


def test_get_submissions_includes_states(db: SubmissionDb) -> None:
    """States relationship is eagerly loaded so it can be accessed outside the session."""
    from grz_db.models.submission import SubmissionStateEnum

    db.add_submission(SUBMISSION_ID)
    db.update_submission_state(SUBMISSION_ID, SubmissionStateEnum.UPLOADED)

    result = db.get_submissions([SUBMISSION_ID])

    assert result[0] is not None
    assert len(result[0].states) >= 1


def _set_tan_g_via_modify(db: SubmissionDb, sub: Submission, tan_g: str) -> None:
    db.modify_submission(sub.id, "tan_g", tan_g)


def _set_tan_g_via_update(db: SubmissionDb, sub: Submission, tan_g: str) -> None:
    sub.tan_g = tan_g
    db.update_submission(sub)


@pytest.mark.parametrize(
    "set_tan_g",
    [_set_tan_g_via_modify, _set_tan_g_via_update],
    ids=["modify_submission", "update_submission"],
)
def test_raises_on_duplicate_tan_g(
    db: SubmissionDb,
    set_tan_g: Callable[[SubmissionDb, Submission, str], None],
) -> None:
    """Raise DuplicateTanGError when tan_g collides with an existing row."""
    sub1 = db.add_submission(SUBMISSION_ID)
    set_tan_g(db, sub1, TAN_G_1)

    sub2 = db.add_submission(SUBMISSION_ID_2)
    with pytest.raises(DuplicateTanGError):
        set_tan_g(db, sub2, TAN_G_1)


def _redacted(metadata: GrzSubmissionMetadata, local_case_id: str) -> GrzSubmissionMetadata:
    """A copy of *metadata* redacted the way archival redacts it.

    Archival writes the ``REDACTED_LOCAL_CASE_ID`` sentinel; older objects in the
    archive carry an empty ``localCaseId``, so both spellings are in the archive and
    both have to be recognised.
    """
    raw = json.loads(metadata.model_dump_json(by_alias=True))
    raw["submission"]["tanG"] = REDACTED_TAN
    raw["submission"]["localCaseId"] = local_case_id
    return GrzSubmissionMetadata.model_validate(raw)


@pytest.mark.parametrize("placeholder", ["", REDACTED_LOCAL_CASE_ID], ids=["empty", "sentinel"])
def test_restore_redacted_fields_uses_the_stored_row(
    db: SubmissionDb, metadata: GrzSubmissionMetadata, placeholder: str
) -> None:
    """Both redaction spellings are restored from the columns populate wrote."""
    submission = db.add_submission(SUBMISSION_ID)
    db.modify_submission(SUBMISSION_ID, "tan_g", TAN_G_1)
    db.modify_submission(SUBMISSION_ID, "pseudonym", "the-real-case")
    submission = db.get_submission(SUBMISSION_ID)

    archived = _redacted(metadata, placeholder)
    unrestored = submission.restore_redacted_fields(archived)

    assert unrestored == frozenset()
    assert archived.submission.tan_g == TAN_G_1
    assert archived.submission.local_case_id == "the-real-case"


@pytest.mark.parametrize("placeholder", ["", REDACTED_LOCAL_CASE_ID], ids=["empty", "sentinel"])
def test_restore_redacted_fields_reports_what_it_cannot_restore(
    db: SubmissionDb, metadata: GrzSubmissionMetadata, placeholder: str
) -> None:
    """A row with nothing stored leaves the placeholders in place and names the columns."""
    submission = db.add_submission(SUBMISSION_ID)

    archived = _redacted(metadata, placeholder)
    unrestored = submission.restore_redacted_fields(archived)

    assert unrestored == frozenset({"tan_g", "pseudonym"})
    assert archived.submission.tan_g == REDACTED_TAN
    assert archived.submission.local_case_id == placeholder


def test_restore_redacted_fields_leaves_unredacted_values_alone(
    db: SubmissionDb, metadata: GrzSubmissionMetadata
) -> None:
    """An unredacted copy is authoritative, so it is diffed rather than overwritten."""
    db.add_submission(SUBMISSION_ID)
    db.modify_submission(SUBMISSION_ID, "tan_g", TAN_G_1)
    db.modify_submission(SUBMISSION_ID, "pseudonym", "stored-case")
    submission = db.get_submission(SUBMISSION_ID)

    unrestored = submission.restore_redacted_fields(metadata)

    assert unrestored == frozenset()
    assert metadata.submission.tan_g != TAN_G_1
    assert metadata.submission.local_case_id != "stored-case"


def test_restore_redacted_fields_restores_each_column_independently(
    db: SubmissionDb, metadata: GrzSubmissionMetadata
) -> None:
    """A stored pseudonym still helps even when the tanG cannot be recovered."""
    db.add_submission(SUBMISSION_ID)
    db.modify_submission(SUBMISSION_ID, "pseudonym", "the-real-case")
    submission = db.get_submission(SUBMISSION_ID)

    archived = _redacted(metadata, "")
    unrestored = submission.restore_redacted_fields(archived)

    assert unrestored == frozenset({"tan_g"})
    assert archived.submission.local_case_id == "the-real-case"


def test_the_schema_is_checked_once_per_database(db: SubmissionDb, monkeypatch: pytest.MonkeyPatch) -> None:
    """The schema check scans the migration directory and opens a connection of its own, and
    each write opens its own session, so caching the answer on the instance is what keeps
    that cost from being paid on every write. The answer cannot change while the process is
    running, so the cache cannot go stale.
    """
    calls = 0
    original = SubmissionDb._at_latest_schema

    def counted(self: SubmissionDb) -> bool:
        nonlocal calls
        calls += 1
        return original(self)

    monkeypatch.setattr(SubmissionDb, "_at_latest_schema", counted)
    db._schema_confirmed = False  # a freshly constructed instance has not asked yet

    db.add_submission(SUBMISSION_ID)
    db.modify_submission(SUBMISSION_ID, "basic_qc_passed", True)
    db.get_submission(SUBMISSION_ID)

    assert calls == 1


def test_a_database_that_is_behind_is_asked_again(db: SubmissionDb, monkeypatch: pytest.MonkeyPatch) -> None:
    """Only a returned answer is remembered, never a raised one.

    ``db init`` and ``db upgrade`` build this object precisely when the schema is behind, so a
    remembered failure would outlive the upgrade that fixed it.
    """
    behind = True
    monkeypatch.setattr(SubmissionDb, "_at_latest_schema", lambda self: not behind)
    db._schema_confirmed = False

    with pytest.raises(OutdatedDatabaseSchemaError):
        db.get_submission(SUBMISSION_ID)

    behind = False
    assert db.get_submission(SUBMISSION_ID) is None, "the upgrade must be picked up"
