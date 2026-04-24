import pytest
from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives.asymmetric import ed25519
from grz_db.models.author import Author
from grz_db.models.submission import SubmissionDb

TWO_TB = 2 * 1024**4  # 2,199,023,255,552 bytes
SUBMISSION_ID = "123456789_2024-01-01_abcdef01"


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
    db.modify_submission(submission_id=SUBMISSION_ID, key="submission_size", value=TWO_TB)

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
    db.modify_submission(submission_id=SUBMISSION_ID, key="submission_metadata", value=metadata)

    result = db.get_submission(SUBMISSION_ID)
    assert result is not None
    assert result.submission_metadata == metadata