import datetime

import pytest
from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives.asymmetric import ed25519
from grz_db.models.author import Author
from grz_db.models.submission import Submission, SubmissionDb
from grz_pydantic_models.submission.metadata import GrzSubmissionMetadata

TWO_TB = 2 * 1024**4  # 2,199,023,255,552 bytes
SUBMISSION_ID = "123456789_2024-01-01_abcdef01"
SUBMISSION_ID_2 = "123456789_2024-01-02_abcdef02"
SUBMISSION_ID_3 = "123456789_2024-01-03_abcdef03"


@pytest.fixture(scope="function")
def test_author() -> Author:
    key = ed25519.Ed25519PrivateKey.generate()
    private_key_bytes = key.private_bytes(
        encoding=serialization.Encoding.PEM,
        format=serialization.PrivateFormat.OpenSSH,
        encryption_algorithm=serialization.NoEncryption(),
    )
    return Author(name="alice", private_key_bytes=private_key_bytes, private_key_passphrase="")


@pytest.fixture(scope="function")
def db(db_test_connection: str, test_author: Author) -> SubmissionDb:
    submission_db = SubmissionDb(db_url=db_test_connection, author=test_author, debug=False)
    submission_db.initialize_schema()
    return submission_db


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
    assert submission.submission_date == explicit_date  # explicit date takes precedence
    assert submission.submission_size == metadata.get_submission_size()
    assert submission.submission_metadata == metadata.to_redacted_dict()

    # --- explicit date is preferred over the metadata date ---
    submission_fallback = Submission.from_metadata(SUBMISSION_ID, metadata, None)
    assert submission_fallback.submission_date == metadata.submission.submission_date

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


def test_get_submissions_all_unknown(db: SubmissionDb) -> None:
    """get_submissions returns all-None when none of the IDs exist."""
    ids = ["000000000_2000-01-01_deadbeef", "000000000_2000-01-01_cafebabe"]
    result = db.get_submissions(ids)

    assert result == [None, None]


def test_get_submissions_includes_states(db: SubmissionDb) -> None:
    """States relationship is eagerly loaded so it can be accessed outside the session."""
    from grz_db.models.submission import SubmissionStateEnum

    db.add_submission(SUBMISSION_ID)
    db.update_submission_state(SUBMISSION_ID, SubmissionStateEnum.UPLOADED)

    result = db.get_submissions([SUBMISSION_ID])

    assert result[0] is not None
    assert len(result[0].states) >= 1
