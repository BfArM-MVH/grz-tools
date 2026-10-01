"""Coverage for the inbox migration (fe9309db9fca).

Seeds a submission at the revision just before the inbox migration, upgrades, and asserts that the
nullable ``inbox`` column appears, that the pre-migration row stays untouched (NULL), and that the
recorded inbox can be read back again.
"""

import sqlalchemy
from grz_db.models.submission import SubmissionDb

PRE_INBOX_MIGRATION_REVISION = "b7e4c1f9a2d3"
INBOX_MIGRATION = "fe9309db9fca"


def test_inbox_migration_adds_the_inbox_column(db_test_connection: str):
    """The migration (fe9309db9fca) adds the nullable ``inbox`` column.

    Submissions recorded before the migration stay untouched (NULL).
    The values can be derived by backfill afterwards.
    The S3 bucket is not stored.
    It follows from the config.
    """
    submission_id = "111111111_2025-01-01_00000001"
    db = SubmissionDb(db_url=db_test_connection, author=None)
    db.upgrade_schema(revision=PRE_INBOX_MIGRATION_REVISION)

    engine = sqlalchemy.create_engine(db_test_connection)
    submissions = sqlalchemy.Table("submissions", sqlalchemy.MetaData(), autoload_with=engine)
    with engine.begin() as conn:
        conn.execute(submissions.insert(), dict(id=submission_id, submission_type="initial"))

    db.upgrade_schema(revision=INBOX_MIGRATION)

    upgraded_submissions = sqlalchemy.Table("submissions", sqlalchemy.MetaData(), autoload_with=engine)
    assert "inbox" in {column.name for column in upgraded_submissions.columns}

    submission = db.get_submission(submission_id)
    assert submission is not None
    assert submission.inbox is None

    db.set_submission_inbox(submission_id, inbox_name="inbox")
    assert (recorded := db.get_submission(submission_id)) is not None
    assert recorded.inbox == "inbox"
