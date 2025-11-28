import datetime
import hashlib
from collections.abc import Generator
from pathlib import Path

import pytest
from grz_db.models.author import Author
from grz_db.models.submission import (
    Submission,
    SubmissionDb,
    SubmissionStateEnum,
    SubmissionStateLog,
)


@pytest.fixture(name="db")
def db_fixture(author: Author, tmp_path: Path) -> Generator[SubmissionDb, None, None]:
    db_file = tmp_path / "test.sqlite"
    db_url = f"sqlite:///{db_file}"
    db = SubmissionDb(db_url=db_url, author=author)
    db.initialize_schema()

    yield db


@pytest.fixture(name="author")
def author_fixture() -> Author:
    return Author(name="test_author", private_key_bytes=b"", private_key_passphrase="")


def add_test_submission(
    db: SubmissionDb,
    author: Author,
    submitter_id: str,
    uploaded_timestamp: datetime.datetime,
    is_validated: bool,
    is_qcing: bool,
):
    """Helper function to populate the in-memory database with a test submission."""
    date_str = uploaded_timestamp.date().isoformat()
    unique_hash = hashlib.sha256(str(uploaded_timestamp.timestamp()).encode()).hexdigest()[:8]
    sub_id = f"{submitter_id}_{date_str}_{unique_hash}"
    with db._get_session() as session:
        submission = Submission(id=sub_id, submitter_id=submitter_id, submission_date=uploaded_timestamp.date())
        session.add(submission)

        uploaded_state = SubmissionStateLog(
            submission_id=sub_id,
            state=SubmissionStateEnum.UPLOADED,
            timestamp=uploaded_timestamp,
            author_name=author.name,
            signature="dummy",
        )
        session.add(uploaded_state)

        if is_validated:
            validated_state = SubmissionStateLog(
                submission_id=sub_id,
                state=SubmissionStateEnum.VALIDATED,
                timestamp=uploaded_timestamp + datetime.timedelta(minutes=1),
                author_name=author.name,
                signature="dummy",
            )
            session.add(validated_state)

        if is_qcing:
            if not is_validated:
                raise ValueError("A submission cannot be in QCing state without being Validated.")
            qcing_state = SubmissionStateLog(
                submission_id=sub_id,
                state=SubmissionStateEnum.QCING,
                timestamp=uploaded_timestamp + datetime.timedelta(minutes=2),
                author_name=author.name,
                signature="dummy",
            )
            session.add(qcing_state)

        session.commit()


class TestGetMonthlyQCStats:
    LE_A_ID = "111111111"
    LE_B_ID = "999999999"

    def test_interspersed_qc(self, db: SubmissionDb, author: Author):
        """
        Tests a standard month with 5 validated submissions, where the 3rd one was picked for QC.
        """
        year, month = 2025, 7

        add_test_submission(
            db, author, self.LE_A_ID, datetime.datetime(year, month, 5, 10, 0), is_validated=True, is_qcing=False
        )
        add_test_submission(
            db, author, self.LE_A_ID, datetime.datetime(year, month, 6, 11, 0), is_validated=True, is_qcing=False
        )
        add_test_submission(
            db,
            author,
            self.LE_A_ID,
            datetime.datetime(year, month, month, 12, 0),
            is_validated=True,
            is_qcing=True,
        )
        add_test_submission(
            db, author, self.LE_A_ID, datetime.datetime(year, month, 8, 13, 0), is_validated=True, is_qcing=False
        )
        add_test_submission(
            db, author, self.LE_A_ID, datetime.datetime(year, month, 9, 14, 0), is_validated=True, is_qcing=False
        )

        validated, qcing, since_last_qcing = db.get_monthly_qc_stats(self.LE_A_ID, year, month)

        assert validated == 5
        assert qcing == 1
        assert since_last_qcing == 2

    def test_no_qc_yet_in_month(self, db: SubmissionDb, author: Author):
        """
        Tests a month where validated submissions exist, but none have been selected for QC.
        """
        year, month = 2025, 7

        add_test_submission(
            db, author, self.LE_A_ID, datetime.datetime(year, month, 5, 10, 0), is_validated=True, is_qcing=False
        )
        add_test_submission(
            db, author, self.LE_A_ID, datetime.datetime(year, month, 6, 11, 0), is_validated=True, is_qcing=False
        )

        validated, qcing, since_last_qcing = db.get_monthly_qc_stats(self.LE_A_ID, year, month)

        assert validated == 2
        assert qcing == 0
        assert since_last_qcing == 2

    def test_empty_month_returns_zero(self, db: SubmissionDb, author: Author):
        validated, qcing, since_last_qcing = db.get_monthly_qc_stats(self.LE_A_ID, 2025, 8)
        assert (validated, qcing, since_last_qcing) == (0, 0, 0)

    def test_last_submission_was_qcing(self, db: SubmissionDb, author: Author):
        year, month = 2025, 7

        add_test_submission(
            db, author, self.LE_A_ID, datetime.datetime(year, month, 5, 10, 0), is_validated=True, is_qcing=False
        )
        add_test_submission(
            db, author, self.LE_A_ID, datetime.datetime(year, month, 6, 11, 0), is_validated=True, is_qcing=True
        )

        validated, qcing, since_last_qcing = db.get_monthly_qc_stats(self.LE_A_ID, year, month)

        assert validated == 2
        assert qcing == 1
        assert since_last_qcing == 0

    def test_multi_le_and_month_isolation(self, db: SubmissionDb, author: Author):
        year = 2025

        # LE_A in month 7
        le = self.LE_A_ID
        month = 7
        add_test_submission(db, author, le, datetime.datetime(year, month, 5, 10, 0), is_validated=True, is_qcing=False)
        add_test_submission(db, author, le, datetime.datetime(year, month, 6, 11, 0), is_validated=True, is_qcing=True)
        add_test_submission(db, author, le, datetime.datetime(year, month, 7, 12, 0), is_validated=True, is_qcing=False)

        # LE_B in the same month (7)
        add_test_submission(
            db, author, self.LE_B_ID, datetime.datetime(year, 7, 8, 13, 0), is_validated=True, is_qcing=True
        )

        # LE_A in a different month (8)
        add_test_submission(
            db, author, self.LE_A_ID, datetime.datetime(year, 8, 1, 10, 0), is_validated=True, is_qcing=True
        )

        # check only LE_A in month 7
        validated, qcing, since_last_qcing = db.get_monthly_qc_stats(self.LE_A_ID, year, 7)

        assert validated == 3
        assert qcing == 1
        assert since_last_qcing == 1

    def test_non_validated_submissions_are_ignored(self, db: SubmissionDb, author: Author):
        year, month = 2025, 7

        add_test_submission(
            db, author, self.LE_A_ID, datetime.datetime(year, month, 5, 10, 0), is_validated=True, is_qcing=False
        )
        add_test_submission(
            db, author, self.LE_A_ID, datetime.datetime(year, month, 6, 11, 0), is_validated=False, is_qcing=False
        )
        add_test_submission(
            db, author, self.LE_A_ID, datetime.datetime(year, month, 7, 12, 0), is_validated=True, is_qcing=False
        )

        validated, qcing, since_last_qcing = db.get_monthly_qc_stats(self.LE_A_ID, year, month)

        assert validated == 2
        assert qcing == 0
        assert since_last_qcing
