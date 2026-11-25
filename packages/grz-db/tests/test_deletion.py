import datetime

import pytest
from grz_db.models.submission import Submission, SubmissionDb, SubmissionStateEnum
from sqlmodel import select

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


def _queried(db: SubmissionDb, submission_id: str) -> tuple[bool, datetime.datetime | None, bool]:
    """``is_deleted``, ``deleted_at`` and ``deletion_started`` as the database computes them."""
    with db.transaction() as session:
        return session.exec(
            select(Submission.is_deleted, Submission.deleted_at, Submission.deletion_started).where(
                Submission.id == submission_id
            )
        ).one()


def _first_deleted_timestamp(submission: Submission) -> datetime.datetime:
    return min(state.timestamp for state in submission.states if state.state == SubmissionStateEnum.DELETED)


class TestDeletionProperties:
    @pytest.mark.parametrize(
        ("history", "is_deleted", "deletion_started"),
        [
            pytest.param([], False, False, id="no-states"),
            pytest.param([SubmissionStateEnum.DOWNLOADED], False, False, id="downloaded"),
            pytest.param([SubmissionStateEnum.DELETING], False, True, id="interrupted"),
            pytest.param([SubmissionStateEnum.DELETING, SubmissionStateEnum.DELETED], True, True, id="deleted"),
            pytest.param(
                [SubmissionStateEnum.DELETED, SubmissionStateEnum.DOWNLOADED], True, True, id="state-after-deleted"
            ),
            pytest.param([SubmissionStateEnum.DELETED, SubmissionStateEnum.DELETED], True, True, id="deleted-twice"),
        ],
    )
    def test_object_and_query_agree(
        self,
        db: SubmissionDb,
        submission_id: str,
        history: list[SubmissionStateEnum],
        is_deleted: bool,
        deletion_started: bool,
    ):
        """Each property gives the same answer on a loaded submission and in a query."""
        for state in history:
            db.update_submission_state(submission_id, state)
        submission = db.get_submission(submission_id)
        deleted_at = _first_deleted_timestamp(submission) if is_deleted else None

        assert (submission.is_deleted, submission.deleted_at, submission.deletion_started) == (
            is_deleted,
            deleted_at,
            deletion_started,
        )
        assert _queried(db, submission_id) == (is_deleted, deleted_at, deletion_started)

    def test_deleted_at_is_the_first_deleted_state(self, db: SubmissionDb, submission_id: str):
        """A second ``Deleted`` state, as ``grzctl delete --check`` records it, leaves ``deleted_at`` unchanged."""
        db.update_submission_state(submission_id, SubmissionStateEnum.DELETED)
        first = db.get_submission(submission_id).get_latest_state().timestamp
        db.update_submission_state(submission_id, SubmissionStateEnum.DELETED)

        submission = db.get_submission(submission_id)
        assert len(submission.states) == 2
        assert submission.deleted_at == first
        assert _queried(db, submission_id)[1] == first
