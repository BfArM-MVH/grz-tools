"""Coverage for the upload timing migration (a7c3e9d1f5b2).

Seeds a submission at the revision just before the migration, upgrades, and asserts that the three
nullable columns appear, that the pre-migration row stays untouched (NULL), and that a recorded
timing can be read back again.
"""

import datetime

import sqlalchemy
from grz_db.models.submission import SubmissionDb

PRE_UPLOAD_TIMING_MIGRATION_REVISION = "fe9309db9fca"
UPLOAD_TIMING_MIGRATION = "a7c3e9d1f5b2"


def test_upload_timing_migration_adds_the_columns(db_test_connection: str):
    """The migration (a7c3e9d1f5b2) adds the nullable upload timing columns.

    Submissions recorded before the migration stay untouched (NULL): their inbox objects are
    gone once ``grzctl clean`` ran, so there is nothing to backfill them from.
    """
    submission_id = "111111111_2025-01-01_00000001"
    db = SubmissionDb(db_url=db_test_connection, author=None)
    db.upgrade_schema(revision=PRE_UPLOAD_TIMING_MIGRATION_REVISION)

    engine = sqlalchemy.create_engine(db_test_connection)
    submissions = sqlalchemy.Table("submissions", sqlalchemy.MetaData(), autoload_with=engine)
    with engine.begin() as conn:
        conn.execute(submissions.insert(), dict(id=submission_id, submission_type="initial"))

    db.upgrade_schema(revision=UPLOAD_TIMING_MIGRATION)

    upgraded_submissions = sqlalchemy.Table("submissions", sqlalchemy.MetaData(), autoload_with=engine)
    columns = {column.name for column in upgraded_submissions.columns}
    assert {"upload_started_at", "upload_finished_at", "upload_resumed"} <= columns

    # the ORM refuses a database behind the head, so bring it there for the reads below
    db.upgrade_schema()
    submission = db.get_submission(submission_id)
    assert submission is not None
    assert (submission.upload_started_at, submission.upload_finished_at, submission.upload_resumed) == (None,) * 3

    started = datetime.datetime(2025, 1, 1, 10, 0, tzinfo=datetime.UTC)
    finished = started + datetime.timedelta(minutes=30)
    db.set_submission_upload_timing(submission_id, started, finished, resumed=False)
    assert (recorded := db.get_submission(submission_id)) is not None
    assert recorded.upload_finished_at - recorded.upload_started_at == datetime.timedelta(minutes=30)
    assert recorded.upload_resumed is False
