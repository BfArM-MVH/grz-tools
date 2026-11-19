import pytest
from grz_db.models.submission import SubmissionDb, SubmissionStateEnum

SUBMITTER_ID = "123456789"


@pytest.fixture(scope="function")
def submission_id(db: SubmissionDb) -> str:
    sid = f"{SUBMITTER_ID}_2025-01-01_00000000"
    db.add_submission(sid)
    return sid


class TestDeletionStates:
    @pytest.mark.parametrize("state", [SubmissionStateEnum.DELETING, SubmissionStateEnum.DELETED])
    def test_roundtrip(self, db: SubmissionDb, submission_id: str, state: SubmissionStateEnum):
        """The deletion states persist, which on PostgreSQL needs them in the submissionstateenum type."""
        db.update_submission_state(submission_id, state)
        assert db.get_submission(submission_id).get_latest_state().state == state
