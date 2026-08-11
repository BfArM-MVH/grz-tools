import datetime
import math

import pytest
from grz_db.errors import (
    SubmissionBasicQCNotPassedError,
    SubmissionDateIsNoneError,
    SubmissionTypeIsNoneError,
)
from grz_db.models.submission import SubmissionDb, SubmissionStateEnum, SubmissionStateLog, SubmissionType
from sqlmodel import Session

SUBMITTER_ID = "123456789"
DEFAULT_HISTORY = ["uploaded", "downloading", "downloaded", "decrypting", "decrypted", "validating", "validated"]


def _update_submission_state(
    db: SubmissionDb, submission_id: str, state: SubmissionStateEnum, timestamp: datetime.datetime
):
    """
    Helper to manually insert a submission state log, to be able to control timestamps.
    """
    with Session(db.engine) as session:
        log = SubmissionStateLog(
            submission_id=submission_id,
            state=state,
            timestamp=timestamp,
            author_name="alice",
            signature="dummy",
        )
        session.add(log)
        session.commit()


def _add_submission_with_history(
    db: SubmissionDb,
    submission_id: str,
    submitter_id: str,
    submission_date: datetime.date,
    states: list[str],
    base_timestamp: datetime.datetime,
    is_qced: bool = False,
):
    """
    Helper to manually insert a submission and its state history.
    """
    db.add_submission(submission_id)

    db.modify_submission(submission_id, "submission_uploaded_date", str(submission_date.isoformat()))
    db.modify_submission(submission_id, "submission_type", SubmissionType.initial)
    db.modify_submission(submission_id, "submitter_id", submitter_id)
    db.modify_submission(submission_id, "basic_qc_passed", "true")

    if is_qced:
        db.modify_submission(submission_id, "detailed_qc_passed", "true")

    current_timestamp = base_timestamp
    for state_str in states:
        state_enum = SubmissionStateEnum(state_str.capitalize())
        _update_submission_state(db, submission_id, state_enum, current_timestamp)
        current_timestamp += datetime.timedelta(seconds=1)


def _add_submission_block(
    db: SubmissionDb, base_date: datetime.date, start_time: datetime.datetime, count: int
) -> list:
    """Add *count* QC candidates ten minutes apart and return them in submission order."""
    submissions = []
    for i in range(count):
        submission_id = f"{SUBMITTER_ID}_{base_date}_{i:0>8}"
        _add_submission_with_history(
            db,
            submission_id,
            SUBMITTER_ID,
            base_date,
            DEFAULT_HISTORY,
            base_timestamp=start_time + datetime.timedelta(minutes=i * 10),
            is_qced=False,
        )
        submission = db.get_submission(submission_id)
        assert submission is not None
        submissions.append(submission)
    return submissions


class TestQcStrategy:
    """Tests the QC selection strategy method db.should_qc."""

    def test_should_qc_returns_existing_selected_for_qc_true(self, db: SubmissionDb):
        test_date = datetime.date(2025, 12, 1)
        base_timestamp = datetime.datetime.combine(test_date, datetime.time(10, 0), tzinfo=datetime.UTC)
        submission_id = f"{SUBMITTER_ID}_{test_date}_00000000"

        _add_submission_with_history(
            db,
            submission_id,
            SUBMITTER_ID,
            test_date,
            DEFAULT_HISTORY,
            base_timestamp=base_timestamp,
            is_qced=False,
        )
        db.modify_submission(submission_id, "selected_for_qc", "true")

        should_run = db.should_qc(
            submission_id=submission_id,
            target_percentage=2.0,
            salt="any_salt",
        )

        assert should_run is True

    def test_should_qc_returns_existing_selected_for_qc_false(self, db: SubmissionDb):
        test_date = datetime.date(2025, 12, 1)
        base_timestamp = datetime.datetime.combine(test_date, datetime.time(10, 0), tzinfo=datetime.UTC)
        submission_id = f"{SUBMITTER_ID}_{test_date}_00000000"

        _add_submission_with_history(
            db,
            submission_id,
            SUBMITTER_ID,
            test_date,
            DEFAULT_HISTORY,
            base_timestamp=base_timestamp,
            is_qced=False,
        )
        db.modify_submission(submission_id, "selected_for_qc", "false")

        should_run = db.should_qc(
            submission_id=submission_id,
            target_percentage=2.0,
            salt="any_salt",
        )

        assert should_run is False

    def test_should_qc_persists_selected_for_qc_result(self, db: SubmissionDb):
        test_date = datetime.date(2025, 12, 1)
        base_timestamp = datetime.datetime.combine(test_date, datetime.time(10, 0), tzinfo=datetime.UTC)
        submission_id = f"{SUBMITTER_ID}_{test_date}_00000000"

        _add_submission_with_history(
            db,
            submission_id,
            SUBMITTER_ID,
            test_date,
            DEFAULT_HISTORY,
            base_timestamp=base_timestamp,
            is_qced=False,
        )

        should_run = db.should_qc(
            submission_id=submission_id,
            target_percentage=2.0,
            salt="any_salt",
        )

        assert should_run is True

        submission = db.get_submission(submission_id)
        assert submission is not None
        assert submission.selected_for_qc is True

    def test_should_qc_persists_a_negative_decision(self, db: SubmissionDb):
        """A zero target leaves no branch that can select, so the stored flag must go to False.

        The month rule ignores the target and fires whenever nothing has been selected yet, so an
        already-selected submission is needed for the decision to reach the target-driven branches.
        """
        base_date = datetime.date(2025, 12, 1)
        start_time = datetime.datetime.combine(base_date, datetime.time(9, 0), tzinfo=datetime.UTC)

        selected_submission_id = f"{SUBMITTER_ID}_{base_date}_00000000"
        _add_submission_with_history(
            db, selected_submission_id, SUBMITTER_ID, base_date, DEFAULT_HISTORY, base_timestamp=start_time
        )
        db.modify_submission(selected_submission_id, "selected_for_qc", "true")

        candidate_submission_id = f"{SUBMITTER_ID}_{base_date}_00000001"
        _add_submission_with_history(
            db,
            candidate_submission_id,
            SUBMITTER_ID,
            base_date,
            DEFAULT_HISTORY,
            base_timestamp=start_time + datetime.timedelta(minutes=10),
        )

        assert db.should_qc(candidate_submission_id, 0.0, "any-salt") is False
        submission = db.get_submission(candidate_submission_id)
        assert submission is not None
        assert submission.selected_for_qc is False

    def test_a_submission_with_qc_history_counts_toward_the_target(self, db: SubmissionDb):
        """A submission already QCed counts even though its selected_for_qc flag was never set.

        Asserted on the counting rule rather than on should_qc's verdict: routing through the full
        decision would leave the outcome resting on whether the random branch happened to fire.
        """
        base_date = datetime.date(2025, 12, 1)
        start_time = datetime.datetime.combine(base_date, datetime.time(9, 0), tzinfo=datetime.UTC)

        qced_submission_id = f"{SUBMITTER_ID}_{base_date}_00000000"
        _add_submission_with_history(
            db,
            qced_submission_id,
            SUBMITTER_ID,
            base_date,
            [*DEFAULT_HISTORY, "qcing", "qced", "cleaning", "cleaned"],
            base_timestamp=start_time,
            is_qced=True,
        )

        untouched_submission_id = f"{SUBMITTER_ID}_{base_date}_00000001"
        _add_submission_with_history(
            db,
            untouched_submission_id,
            SUBMITTER_ID,
            base_date,
            DEFAULT_HISTORY,
            base_timestamp=start_time + datetime.timedelta(minutes=10),
        )

        qced = db.get_submission(qced_submission_id)
        untouched = db.get_submission(untouched_submission_id)
        assert qced is not None and untouched is not None
        assert qced.selected_for_qc is not True, "the flag is unset; only the state history should count it"
        assert db._submission_counts_as_selected_for_qc(qced) is True
        assert db._submission_counts_as_selected_for_qc(untouched) is False

    def test_should_qc_still_selects_once_per_block_when_the_target_is_already_met(self, db: SubmissionDb):
        """With the target already exceeded, only the random branch can select, so it must still fire.

        The month and quarter rules are both satisfied here, so this is the one place that proves
        ``should_qc`` actually consults the random selection rather than merely returning early.
        """
        target_percentage = 20.0
        block_size = math.floor(1 / (target_percentage / 100.0))
        base_date = datetime.date(2025, 12, 1)
        start_time = datetime.datetime.combine(base_date, datetime.time(9, 0), tzinfo=datetime.UTC)

        submissions = _add_submission_block(db, base_date, start_time, count=block_size * 2)
        # push the selected share far above the target so neither the month nor the quarter rule fires
        for submission in submissions[:block_size]:
            db.modify_submission(submission.id, "selected_for_qc", "true")

        selected = [
            submission.id
            for submission in submissions[block_size:]
            if db.should_qc(submission.id, target_percentage, "test-salt")
        ]
        assert len(selected) == 1, f"expected the random branch to pick exactly one, picked {selected}"

    def test_ratio_catchup_keeps_the_selected_share_at_or_above_target(self, db: SubmissionDb):
        """The running selected share must never drop below the target; that is what catch-up means.

        Asserted as a property of the outcome rather than by recomputing the selection formula:
        a test that re-derives the seed would agree with a broken implementation.
        """
        target_percentage = 20.0
        target_proportion = target_percentage / 100.0
        salt = "ratio-test"
        base_date = datetime.date(2025, 12, 1)
        start_time = datetime.datetime.combine(base_date, datetime.time(9, 0), tzinfo=datetime.UTC)

        selected_count = 0
        for i in range(12):
            submission_id = f"{SUBMITTER_ID}_{base_date}_{i:0>8}"
            submission_timestamp = start_time + datetime.timedelta(minutes=i * 10)
            _add_submission_with_history(
                db,
                submission_id,
                SUBMITTER_ID,
                base_date,
                DEFAULT_HISTORY,
                base_timestamp=submission_timestamp,
                is_qced=False,
            )

            if db.should_qc(submission_id, target_percentage, salt):
                _update_submission_state(
                    db, submission_id, SubmissionStateEnum.QCING, submission_timestamp + datetime.timedelta(minutes=1)
                )
                _update_submission_state(
                    db, submission_id, SubmissionStateEnum.QCED, submission_timestamp + datetime.timedelta(minutes=2)
                )
                db.modify_submission(submission_id, "detailed_qc_passed", "true")
                selected_count += 1

            assert selected_count / (i + 1) >= target_proportion, (
                f"after {i + 1} submissions only {selected_count} were selected, "
                f"below the {target_proportion:.0%} target"
            )

    @pytest.mark.parametrize("target_percentage", [2.0, 20.0])
    def test_exactly_one_submission_per_block_is_randomly_selected(self, db: SubmissionDb, target_percentage: float):
        """Each block of ``floor(1 / target)`` submissions must contribute exactly one selection.

        This is the guarantee the block scheme exists to provide, and it holds whatever seed the
        implementation derives. Recomputing that seed here would instead make the test agree with
        a broken implementation, since both sides would share the same mistake.
        """
        target_proportion = target_percentage / 100.0
        salt = "test-salt"
        base_date = datetime.date(2025, 7, 1)
        start_time = datetime.datetime.combine(base_date, datetime.time(8, 0), tzinfo=datetime.UTC)
        block_size = math.floor(1 / target_proportion)

        submissions = _add_submission_block(db, base_date, start_time, count=block_size * 2)

        for block in range(2):
            offset = block * block_size
            selected = [
                index
                for index in range(offset, offset + block_size)
                if db._is_randomly_selected_for_qc(submissions[index], submissions, target_proportion, salt)
            ]
            assert len(selected) == 1, f"block {block} selected {selected}, expected exactly one submission"

    def test_random_selection_is_reproducible(self, db: SubmissionDb):
        """The selection is documented as random but reproducible, so the same inputs must agree."""
        target_proportion = 0.02
        salt = "test-salt"
        base_date = datetime.date(2025, 7, 1)
        start_time = datetime.datetime.combine(base_date, datetime.time(8, 0), tzinfo=datetime.UTC)
        block_size = math.floor(1 / target_proportion)

        submissions = _add_submission_block(db, base_date, start_time, count=block_size)

        first = [db._is_randomly_selected_for_qc(s, submissions, target_proportion, salt) for s in submissions]
        second = [db._is_randomly_selected_for_qc(s, submissions, target_proportion, salt) for s in submissions]
        assert first == second

    def test_no_submission_is_randomly_selected_at_zero_target(self, db: SubmissionDb):
        """A zero target short-circuits the random branch, so nothing may be selected by it."""
        base_date = datetime.date(2025, 7, 1)
        start_time = datetime.datetime.combine(base_date, datetime.time(8, 0), tzinfo=datetime.UTC)

        submissions = _add_submission_block(db, base_date, start_time, count=5)

        assert not any(db._is_randomly_selected_for_qc(s, submissions, 0.0, "test-salt") for s in submissions)

    def test_should_qc_raises_on_missing_submission_date(self, db: SubmissionDb):
        """Verify SubmissionDateIsNoneError when submission_date is None."""
        submission_id = f"{SUBMITTER_ID}_2025-12-01_00000000"
        db.add_submission(submission_id)
        db.modify_submission(submission_id, "submission_type", SubmissionType.initial.value)
        db.modify_submission(submission_id, "basic_qc_passed", "true")

        with pytest.raises(SubmissionDateIsNoneError):
            db.should_qc(submission_id, 2.0, "salt")

    def test_should_qc_raises_on_missing_submission_type(self, db: SubmissionDb):
        """Verify SubmissionTypeIsNoneError when submission_type is None."""
        submission_id = f"{SUBMITTER_ID}_2025-12-01_00000000"
        db.add_submission(submission_id)
        db.modify_submission(submission_id, "submission_uploaded_date", "2025-12-01")
        db.modify_submission(submission_id, "basic_qc_passed", "true")

        with pytest.raises(SubmissionTypeIsNoneError):
            db.should_qc(submission_id, 2.0, "salt")

    def test_should_qc_raises_on_failed_basic_qc(self, db: SubmissionDb):
        """Verify SubmissionBasicQCNotPassedError when basic_qc_passed is False."""
        test_date = datetime.date(2025, 12, 1)
        submission_id = f"{SUBMITTER_ID}_{test_date}_00000000"
        _add_submission_with_history(
            db,
            submission_id,
            SUBMITTER_ID,
            test_date,
            DEFAULT_HISTORY,
            base_timestamp=datetime.datetime.combine(test_date, datetime.time(10, 0), tzinfo=datetime.UTC),
        )
        db.modify_submission(submission_id, "basic_qc_passed", "false")

        with pytest.raises(SubmissionBasicQCNotPassedError):
            db.should_qc(submission_id, 2.0, "salt")
