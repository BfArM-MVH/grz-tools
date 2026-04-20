import pytest
from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives.asymmetric import ed25519
from grz_db.models.author import Author
from grz_db.models.submission import FailureReasonEnum, SubmissionDb, SubmissionStateEnum

SUBMITTER_ID = "123456789"


@pytest.fixture(scope="function")
def test_author() -> Author:
    key = ed25519.Ed25519PrivateKey.generate()
    private_key_bytes = key.private_bytes(
        encoding=serialization.Encoding.PEM,
        format=serialization.PrivateFormat.OpenSSH,
        encryption_algorithm=serialization.NoEncryption(),
    )
    return Author(
        name="alice",
        private_key_bytes=private_key_bytes,
        private_key_passphrase="",
    )


@pytest.fixture(scope="function")
def db(tmp_path, test_author: Author) -> SubmissionDb:
    db_path = tmp_path / "submissions.sqlite"
    db_url = f"sqlite:///{db_path.resolve()}"
    submission_db = SubmissionDb(db_url=db_url, author=test_author, debug=False)
    submission_db.initialize_schema()
    return submission_db


@pytest.fixture(scope="function")
def submission_id(db: SubmissionDb) -> str:
    sid = f"{SUBMITTER_ID}_2025-01-01_00000000"
    db.add_submission(sid)
    return sid


class TestFailureReasonEnum:
    def test_list_returns_all_values(self):
        values = FailureReasonEnum.list()
        assert set(values) == {e.value for e in FailureReasonEnum}

    def test_case_insensitive_lookup(self):
        assert FailureReasonEnum("DECRYPTIONERROR") == FailureReasonEnum.DECRYPTION_ERROR
        assert FailureReasonEnum("DecryptionError") == FailureReasonEnum.DECRYPTION_ERROR

    def test_invalid_value_raises(self):
        with pytest.raises(ValueError):
            FailureReasonEnum("notavalidreason")


class TestFailureReasonPersistence:
    @pytest.mark.parametrize("failure_reason", list(FailureReasonEnum))
    def test_roundtrip(self, db: SubmissionDb, submission_id: str, failure_reason: FailureReasonEnum):
        """Each FailureReasonEnum value should persist and be retrieved correctly."""
        db.update_submission_state(
            submission_id,
            SubmissionStateEnum.ERROR,
            failure_reason=failure_reason,
        )
        submission = db.get_submission(submission_id)
        latest_state = submission.get_latest_state()
        assert latest_state.state == SubmissionStateEnum.ERROR
        assert latest_state.failure_reason == failure_reason

    def test_failure_reason_none_by_default(self, db: SubmissionDb, submission_id: str):
        """failure_reason should be None when not provided."""
        db.update_submission_state(submission_id, SubmissionStateEnum.UPLOADED)
        submission = db.get_submission(submission_id)
        latest_state = submission.get_latest_state()
        assert latest_state.failure_reason is None

    def test_failure_reason_only_on_latest_state(self, db: SubmissionDb, submission_id: str):
        """failure_reason should only be set on the Error state, not bleed into subsequent states."""
        db.update_submission_state(
            submission_id,
            SubmissionStateEnum.ERROR,
            failure_reason=FailureReasonEnum.DECRYPTION_ERROR,
        )
        db.update_submission_state(submission_id, SubmissionStateEnum.UPLOADED)
        submission = db.get_submission(submission_id)
        states = sorted(submission.states, key=lambda s: s.timestamp)
        assert states[-1].state == SubmissionStateEnum.UPLOADED
        assert states[-1].failure_reason is None
        assert states[-2].failure_reason == FailureReasonEnum.DECRYPTION_ERROR
