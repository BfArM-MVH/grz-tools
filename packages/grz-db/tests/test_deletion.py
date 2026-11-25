import datetime
from collections.abc import Callable

import pytest
from grz_db.errors import SubmissionDeletedError
from grz_db.models.submission import (
    ChangeRequestEnum,
    DetailedQCResult,
    Submission,
    SubmissionDb,
    SubmissionStateEnum,
    SubmissionStateLog,
    SubmissionType,
)
from grz_pydantic_models.submission.metadata import (
    GrzSubmissionMetadata,
    LibraryType,
    SequenceSubtype,
    SequenceType,
)
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

    def test_state_after_deleted_keeps_it_deleted(self, db: SubmissionDb, submission_id: str):
        """A state recorded after ``Deleted`` leaves the submission deleted.

        grz-db refuses such a state, so it is written past that guard, as a grzctl from before
        the guard could have done.
        """
        db.update_submission_state(submission_id, SubmissionStateEnum.DELETED)
        first = db.get_submission(submission_id).get_latest_state().timestamp
        with db.transaction() as session:
            session.add(
                SubmissionStateLog(
                    submission_id=submission_id,
                    state=SubmissionStateEnum.DOWNLOADED,
                    author_name="alice",
                    signature="unsigned",
                )
            )

        submission = db.get_submission(submission_id)
        assert submission.get_latest_state().state == SubmissionStateEnum.DOWNLOADED
        assert (submission.is_deleted, submission.deleted_at, submission.deletion_started) == (True, first, True)
        assert _queried(db, submission_id) == (True, first, True)

    def test_deleted_at_is_the_first_deleted_state(self, db: SubmissionDb, submission_id: str):
        """A second ``Deleted`` state, as ``grzctl delete --check`` records it, leaves ``deleted_at`` unchanged."""
        db.update_submission_state(submission_id, SubmissionStateEnum.DELETED)
        first = db.get_submission(submission_id).get_latest_state().timestamp
        db.update_submission_state(submission_id, SubmissionStateEnum.DELETED)

        submission = db.get_submission(submission_id)
        assert len(submission.states) == 2
        assert submission.deleted_at == first
        assert _queried(db, submission_id)[1] == first


def _qc_result(submission_id: str) -> DetailedQCResult:
    return DetailedQCResult(
        submission_id=submission_id,
        lab_datum_id="0",
        pseudonym="index",
        sequence_type=SequenceType.dna,
        sequence_subtype=SequenceSubtype.germline,
        library_type=LibraryType.wes,
        percent_bases_above_quality_threshold_minimum_quality=30.0,
        percent_bases_above_quality_threshold_percent=90.0,
        percent_bases_above_quality_threshold_passed_qc=True,
        percent_bases_above_quality_threshold_percent_deviation=0.0,
        mean_depth_of_coverage=100.0,
        mean_depth_of_coverage_passed_qc=True,
        mean_depth_of_coverage_percent_deviation=0.0,
        targeted_regions_min_coverage=20.0,
        targeted_regions_above_min_coverage=0.9,
        targeted_regions_above_min_coverage_passed_qc=True,
        targeted_regions_above_min_coverage_percent_deviation=0.0,
    )


def _update_submission(db: SubmissionDb, submission_id: str) -> None:
    submission = db.get_submission(submission_id)
    submission.local_case_id = "changed"
    db.update_submission(submission)


GUARDED_WRITES: dict[str, Callable[[SubmissionDb, str, GrzSubmissionMetadata, int], object]] = {
    "update_submission_state": lambda db, sid, _metadata, _case: db.update_submission_state(
        sid, SubmissionStateEnum.DOWNLOADING
    ),
    "modify_submission": lambda db, sid, _metadata, _case: db.modify_submission(sid, "local_case_id", "changed"),
    "set_selected_for_qc": lambda db, sid, _metadata, _case: db.set_selected_for_qc(sid, True),
    "update_submission": lambda db, sid, _metadata, _case: _update_submission(db, sid),
    "set_submission_inbox": lambda db, sid, _metadata, _case: db.set_submission_inbox(sid, "inbox-a"),
    "set_submission_case": lambda db, sid, _metadata, case_id: db.set_submission_case(sid, case_id),
    "populate": lambda db, sid, metadata, _case: db.populate(sid, metadata, submission_date=None),
    "commit_changes": lambda db, sid, metadata, _case: db.commit_changes(
        sid, db.diff(sid, metadata, submission_uploaded_date=None)
    ),
    "add_detailed_qc_result": lambda db, sid, _metadata, _case: db.add_detailed_qc_result(_qc_result(sid)),
    "add_detailed_qc_results": lambda db, sid, _metadata, _case: db.add_detailed_qc_results([_qc_result(sid)]),
}


@pytest.fixture(scope="function")
def metadata_submission_id(db: SubmissionDb, metadata: GrzSubmissionMetadata) -> str:
    """A submission, not yet populated, with the ID of the ``metadata`` example."""
    sid = metadata.submission_id
    db.add_submission(sid)
    return sid


class TestWritesRefusedOnceDeletionStarted:
    @pytest.mark.parametrize(
        "history",
        [
            pytest.param([SubmissionStateEnum.DELETING], id="interrupted"),
            pytest.param([SubmissionStateEnum.DELETING, SubmissionStateEnum.DELETED], id="deleted"),
        ],
    )
    @pytest.mark.parametrize("write", GUARDED_WRITES)
    def test_refused_and_nothing_written(
        self,
        db: SubmissionDb,
        metadata: GrzSubmissionMetadata,
        metadata_submission_id: str,
        history: list[SubmissionStateEnum],
        write: str,
    ):
        """Each write path raises SubmissionDeletedError and leaves the submission as it was."""
        case_id = db.create_case(SUBMITTER_ID, "case-1").id
        for state in history:
            db.update_submission_state(metadata_submission_id, state)
        before = db.get_submission(metadata_submission_id)

        with pytest.raises(SubmissionDeletedError):
            GUARDED_WRITES[write](db, metadata_submission_id, metadata, case_id)

        after = db.get_submission(metadata_submission_id)
        assert after.model_dump() == before.model_dump()
        assert len(after.states) == len(history)
        assert db.get_detailed_qc_results(metadata_submission_id) == ()


class TestWritesAllowedOnceDeletionStarted:
    def test_deletion_states_can_still_be_written(self, db: SubmissionDb, submission_id: str):
        db.update_submission_state(submission_id, SubmissionStateEnum.DELETING)
        db.update_submission_state(submission_id, SubmissionStateEnum.DELETED)
        db.update_submission_state(submission_id, SubmissionStateEnum.DELETED)
        assert len(db.get_submission(submission_id).states) == 3

    def test_change_request_can_still_be_recorded(self, db: SubmissionDb, submission_id: str):
        db.update_submission_state(submission_id, SubmissionStateEnum.DELETED)
        db.add_change_request(
            submission_id,
            ChangeRequestEnum.DELETE,
            requester_name="Erika Mustermann",
            requester_email="erika.mustermann@example.org",
            requested_at=datetime.date(2026, 11, 1),
            request_email_content="Bitte löschen.",
        )
        assert len(db.list_change_requests()) == 1

    def test_case_can_still_be_unlinked(self, db: SubmissionDb, submission_id: str):
        db.modify_submission(submission_id, "submission_type", SubmissionType.initial.value)
        case_id = db.create_case(SUBMITTER_ID, "case-1").id
        db.set_submission_case(submission_id, case_id)
        db.update_submission_state(submission_id, SubmissionStateEnum.DELETED)

        db.clear_submission_case(submission_id)
        assert db.get_submission(submission_id).case_id is None
