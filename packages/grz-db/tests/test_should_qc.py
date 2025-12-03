import datetime
import math
import random
from pathlib import Path

import pytest
from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives.asymmetric import ed25519
from grz_db.models.author import Author
from grz_db.models.submission import SubmissionDb, SubmissionStateEnum, SubmissionStateLog, SubmissionType
from sqlmodel import Session

SUBMITTER_ID = "123456789"
DEFAULT_HISTORY = ["uploaded", "downloading", "downloaded", "decrypting", "decrypted", "validating", "validated"]


@pytest.fixture(scope="function")
def test_author() -> Author:
    """Creates a test author with a valid temporary private key."""
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
def db(tmp_path: Path, test_author: Author) -> SubmissionDb:
    """Create a clean database initialized with schema for each test function."""
    db_path = tmp_path / "submissions.sqlite"
    db_url = f"sqlite:///{db_path.resolve()}"

    submission_db = SubmissionDb(db_url=db_url, author=test_author, debug=False)
    submission_db.initialize_schema()

    return submission_db


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

    db.modify_submission(submission_id, "submission_date", str(submission_date.isoformat()))
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


class TestQcStrategy:
    """Tests the QC selection strategy method db.should_qc."""

    def test_first_of_month_always_runs(self, db: SubmissionDb):
        """
        Test that the first (validated, initial) submission of a month is always QCed.
        """
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
        assert should_run is True, "The first submission of the month must be selected for QC."

    def test_quarterly_ratio_catchup(self, db: SubmissionDb):
        """
        Tests that the quarterly target is met.
        Target: 20%.
        """
        target_percentage = 20.0
        salt = "ratio-test"
        base_date = datetime.date(2025, 12, 1)
        start_time = datetime.datetime.combine(base_date, datetime.time(9, 0), tzinfo=datetime.UTC)

        submission_id = f"{SUBMITTER_ID}_{base_date}_00000000"
        _add_submission_with_history(
            db,
            submission_id,
            SUBMITTER_ID,
            base_date,
            [*DEFAULT_HISTORY, "qcing", "qced"],
            base_timestamp=start_time,
            is_qced=True,
        )

        for i in range(1, 10):
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

            should_run = db.should_qc(submission_id, target_percentage, salt)

            # Ratio hits 20% exactly at indices 4 (1/5) and 9 (2/10)
            if i in (4, 9):
                assert should_run is True, f"Index {i} should have been selected for QC (Ratio catchup)"

                _update_submission_state(
                    db, submission_id, SubmissionStateEnum.QCING, submission_timestamp + datetime.timedelta(minutes=1)
                )
                _update_submission_state(
                    db, submission_id, SubmissionStateEnum.QCED, submission_timestamp + datetime.timedelta(minutes=2)
                )

                db.modify_submission(submission_id, "detailed_qc_passed", "true")
            else:
                assert should_run is False, f"Index {i} should NOT have been selected for QC"

    @pytest.mark.parametrize("target_percentage", [1.0, 2.0, 4.0, 5.0, 10.0, 20.0, 100.0])
    def test_random_selection(self, db: SubmissionDb, target_percentage: float):
        """
        Test the random selection logic.
        We simulate the logic (Month -> Ratio -> Random) to verify expectation.
        """
        salt = "test-salt"
        base_date = datetime.date(2025, 7, 1)
        start_time = datetime.datetime.combine(base_date, datetime.time(8, 0), tzinfo=datetime.UTC)

        block_size = math.floor(1 / (target_percentage / 100.0))
        limit = block_size

        qced_count = 0
        total_count = 0

        for i in range(limit):
            submission_id = f"{SUBMITTER_ID}_{base_date}_{i:0>8}"
            submission_timestamp = start_time + datetime.timedelta(minutes=i * 10)
            total_count += 1

            is_first_of_month_trigger = qced_count == 0

            current_ratio = qced_count / total_count
            is_ratio_trigger = current_ratio <= (target_percentage / 100.0)

            block_index = total_count // block_size
            seed = f"{SUBMITTER_ID}-{base_date.year}-3-{block_index}-{salt}"
            rng = random.Random(seed)
            target_index_in_block = rng.randint(0, block_size - 1)

            is_random_trigger = i == target_index_in_block

            expect_qc = False
            if is_first_of_month_trigger:
                expect_qc = True
            elif is_ratio_trigger:
                expect_qc = True
            elif is_random_trigger:
                expect_qc = True

            _add_submission_with_history(
                db,
                submission_id,
                SUBMITTER_ID,
                base_date,
                DEFAULT_HISTORY,
                base_timestamp=submission_timestamp,
                is_qced=False,
            )

            should_run = db.should_qc(submission_id=submission_id, target_percentage=target_percentage, salt=salt)

            assert should_run is expect_qc, (
                f"Index {i}: Expect={expect_qc}, Got={should_run}. "
                f"(Ratio: {current_ratio:.3f} vs {target_percentage / 100.0}, "
                f"Stats: QCed={qced_count}, Total={total_count})"
            )

            if should_run:
                _update_submission_state(
                    db, submission_id, SubmissionStateEnum.QCING, submission_timestamp + datetime.timedelta(minutes=1)
                )
                _update_submission_state(
                    db, submission_id, SubmissionStateEnum.QCED, submission_timestamp + datetime.timedelta(minutes=2)
                )
                db.modify_submission(submission_id, "detailed_qc_passed", "true")
                qced_count += 1
