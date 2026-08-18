"""Backfill coverage for the cases migration (f8c1a4b7e2d9).

Seeds submissions at the revision just before the cases migration, upgrades, and asserts that the
backfill groups by ``(submitter_id, pseudonym)``, stores the keys on the case, links submissions,
keeps distinct submitters apart, and leaves unkeyable rows unlinked.
"""

import pytest
import sqlalchemy
from grz_db.models.submission import SubmissionDb

PRE_CASES_REVISION = "c8d4a7f2b1e6"
CASES_REVISION = "f8c1a4b7e2d9"


def _tan(i: int) -> str:
    return f"{i:064x}"


def test_cases_backfill_groups_by_submitter_and_local_case(db_test_connection: str):
    db = SubmissionDb(db_url=db_test_connection, author=None)
    db.upgrade_schema(revision=PRE_CASES_REVISION)

    engine = sqlalchemy.create_engine(db_test_connection)
    submissions = sqlalchemy.Table("submissions", sqlalchemy.MetaData(), autoload_with=engine)

    shared_a1 = "111111111_2025-01-01_00000001"
    shared_a2 = "111111111_2025-01-02_00000002"
    other_submitter = "222222222_2025-01-01_00000003"
    null_pseudonym = "111111111_2025-01-03_00000004"
    null_submitter = "333333333_2025-01-01_00000005"
    redacted_1 = "111111111_2025-01-04_00000006"
    redacted_2 = "111111111_2025-01-05_00000007"
    empty_pseudonym = "111111111_2025-01-06_00000008"
    test_type = "111111111_2025-01-07_00000009"

    redacted = "REDACTED_LOCAL_CASE_ID"
    rows = [
        # two submissions, same (submitter, local case id) -> one case, both linked
        dict(id=shared_a1, tan_g=_tan(1), pseudonym="caseX", submitter_id="111111111", submission_type="initial"),
        dict(id=shared_a2, tan_g=_tan(2), pseudonym="caseX", submitter_id="111111111", submission_type="addition"),
        # different submitter, same local case id -> a distinct case
        dict(id=other_submitter, tan_g=_tan(3), pseudonym="caseX", submitter_id="222222222", submission_type="initial"),
        # no pseudonym -> unlinked
        dict(id=null_pseudonym, tan_g=_tan(4), pseudonym=None, submitter_id="111111111", submission_type="initial"),
        # pseudonym but no submitter -> cannot form the key, unlinked
        dict(id=null_submitter, tan_g=_tan(5), pseudonym="caseY", submitter_id=None, submission_type="initial"),
        # redaction placeholders are not real case keys: two redacted 'initial' rows for the same
        # submitter must not be grouped into one case, and must not falsely trigger the migration's
        # duplicate-detection query either
        dict(id=redacted_1, tan_g=_tan(6), pseudonym=redacted, submitter_id="111111111", submission_type="initial"),
        dict(id=redacted_2, tan_g=_tan(7), pseudonym=redacted, submitter_id="111111111", submission_type="initial"),
        dict(id=empty_pseudonym, tan_g=_tan(8), pseudonym="", submitter_id="111111111", submission_type="initial"),
        # test submissions are never case-tracked, even with a full (submitter, local case id) key
        dict(id=test_type, tan_g=_tan(9), pseudonym="caseX", submitter_id="111111111", submission_type="test"),
    ]
    with engine.begin() as conn:
        conn.execute(submissions.insert(), rows)

    db.upgrade_schema(revision=CASES_REVISION)

    case_id_a1 = db.get_submission(shared_a1).case_id
    case_id_a2 = db.get_submission(shared_a2).case_id
    case_id_b = db.get_submission(other_submitter).case_id
    assert case_id_a1 is not None and case_id_b is not None
    assert case_id_a1 == case_id_a2  # shared (submitter, local case id) -> same case
    assert case_id_a1 != case_id_b  # distinct submitter -> distinct case

    # the case stores the resolution keys, psn is not assigned yet
    case_a = db.get_case(case_id_a1)
    assert case_a.submitter_id == "111111111" and case_a.local_case_id == "caseX"
    assert case_a.psn is None and db.get_case(case_id_b).psn is None

    assert {s.id for s in db.list_submissions_for_case(case_id_a1)} == {shared_a1, shared_a2}
    assert {s.id for s in db.list_submissions_for_case(case_id_b)} == {other_submitter}

    # unkeyable rows stay unlinked, and no spurious case is created for them
    assert db.get_submission(null_pseudonym).case_id is None
    assert db.get_submission(null_submitter).case_id is None
    assert db.get_submission(redacted_1).case_id is None
    assert db.get_submission(redacted_2).case_id is None
    assert db.get_submission(empty_pseudonym).case_id is None
    assert db.get_submission(test_type).case_id is None  # despite sharing case A's key

    assert len(db.list_cases()) == 2


def test_cases_backfill_rejects_duplicate_qc_passed_initials(db_test_connection: str):
    """Two QC-passed 'initial' submissions sharing (submitter_id, pseudonym) are genuinely inconsistent.

    A case may have at most one initial submission that passed basic QC, so the migration must
    fail fast instead of silently choosing one of the two as the case's initial submission. This
    exercises the duplicate-detection query (COUNT(*) ... GROUP BY ... HAVING COUNT(*) > 1).
    """
    db = SubmissionDb(db_url=db_test_connection, author=None)
    db.upgrade_schema(revision=PRE_CASES_REVISION)

    engine = sqlalchemy.create_engine(db_test_connection)
    submissions = sqlalchemy.Table("submissions", sqlalchemy.MetaData(), autoload_with=engine)

    initial_a = "111111111_2025-01-01_00000001"
    initial_b = "111111111_2025-01-02_00000002"
    rows = [
        # two QC-passed 'initial' submissions for the same (submitter, local case id) -> inconsistent
        dict(
            id=initial_a,
            tan_g=_tan(1),
            pseudonym="caseX",
            submitter_id="111111111",
            submission_type="initial",
            basic_qc_passed=True,
        ),
        dict(
            id=initial_b,
            tan_g=_tan(2),
            pseudonym="caseX",
            submitter_id="111111111",
            submission_type="initial",
            basic_qc_passed=True,
        ),
    ]
    with engine.begin() as conn:
        conn.execute(submissions.insert(), rows)

    with pytest.raises(RuntimeError, match="more than one QC-passed 'initial' submission shares"):
        db.upgrade_schema(revision=CASES_REVISION)

    # the failed upgrade must not have created the cases table
    inspector = sqlalchemy.inspect(engine)
    assert "cases" not in inspector.get_table_names()


def test_cases_backfill_links_competing_initials_that_did_not_all_pass_qc(db_test_connection: str):
    """Legacy rows sharing one (submitter_id, pseudonym) key, where basic QC did not pass for all
    of them, still link to a single case; the migration must recognize the QC-passed row among
    them as the case's initial submission.
    """
    db = SubmissionDb(db_url=db_test_connection, author=None)
    db.upgrade_schema(revision=PRE_CASES_REVISION)

    engine = sqlalchemy.create_engine(db_test_connection)
    submissions = sqlalchemy.Table("submissions", sqlalchemy.MetaData(), autoload_with=engine)

    failed = "111111111_2025-01-01_00000001"
    corrected = "111111111_2025-01-02_00000002"
    pending = "111111111_2025-01-03_00000003"
    rows = [
        dict(
            id=failed,
            tan_g=_tan(1),
            pseudonym="caseX",
            submitter_id="111111111",
            submission_type="initial",
            basic_qc_passed=False,
        ),
        dict(
            id=corrected,
            tan_g=_tan(2),
            pseudonym="caseX",
            submitter_id="111111111",
            submission_type="initial",
            basic_qc_passed=True,
        ),
        dict(
            id=pending,
            tan_g=_tan(3),
            pseudonym="caseX",
            submitter_id="111111111",
            submission_type="initial",
            basic_qc_passed=None,
        ),
    ]
    with engine.begin() as conn:
        conn.execute(submissions.insert(), rows)

    db.upgrade_schema(revision=CASES_REVISION)

    case_ids = {db.get_submission(sid).case_id for sid in (failed, corrected, pending)}
    assert len(case_ids) == 1 and None not in case_ids
    (case_id,) = case_ids
    assert db.get_initial_submission(case_id).id == corrected


def test_cases_migration_extends_failure_reason_enum(db_test_connection: str):
    """The migration must add 'duplicate_initial' to the PostgreSQL failurereasonenum type, the
    value recorded when a duplicate initial submission fails basic QC.
    """
    db = SubmissionDb(db_url=db_test_connection, author=None)
    db.upgrade_schema(revision=CASES_REVISION)

    engine = sqlalchemy.create_engine(db_test_connection)
    if engine.dialect.name != "postgresql":
        pytest.skip("SQLite stores the enum as plain VARCHAR; only PostgreSQL needs the type extended")
    with engine.connect() as conn:
        labels = conn.execute(
            sqlalchemy.text(
                "SELECT enumlabel FROM pg_enum JOIN pg_type ON pg_enum.enumtypid = pg_type.oid "
                "WHERE pg_type.typname = 'failurereasonenum'"
            )
        ).scalars()
        assert "duplicate_initial" in set(labels)
