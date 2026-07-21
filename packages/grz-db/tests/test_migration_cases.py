"""Backfill coverage for the cases migration (f8c1a4b7e2d9).

Seeds submissions at the revision just before the cases migration, upgrades, and asserts that the
backfill groups by ``(submitter_id, pseudonym)``, stores the keys on the case, links submissions,
keeps distinct submitters apart, and leaves unkeyable rows unlinked.
"""

from pathlib import Path

import pytest
import sqlalchemy
from grz_db.models.submission import SubmissionDb

PRE_CASES_REVISION = "66f36abbea34"
CASES_REVISION = "f8c1a4b7e2d9"


def _tan(i: int) -> str:
    return f"{i:064x}"


def test_cases_backfill_groups_by_submitter_and_local_case(tmp_path: Path):
    db_url = f"sqlite:///{(tmp_path / 'migrate.sqlite').resolve()}"
    db = SubmissionDb(db_url=db_url, author=None)
    db.upgrade_schema(revision=PRE_CASES_REVISION)

    engine = sqlalchemy.create_engine(db_url)
    submissions = sqlalchemy.Table("submissions", sqlalchemy.MetaData(), autoload_with=engine)

    shared_a1 = "111111111_2025-01-01_00000001"
    shared_a2 = "111111111_2025-01-02_00000002"
    other_submitter = "222222222_2025-01-01_00000003"
    null_pseudonym = "111111111_2025-01-03_00000004"
    null_submitter = "333333333_2025-01-01_00000005"

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

    assert len(db.list_cases()) == 2


def test_cases_backfill_rejects_duplicate_initials(tmp_path: Path):
    """Two 'initial' submissions sharing (submitter_id, pseudonym) are genuinely inconsistent.

    The migration must fail fast rather than silently pick a winner, exercising the
    duplicate-detection query (COUNT(*) ... GROUP BY ... HAVING COUNT(*) > 1).
    """
    db_url = f"sqlite:///{(tmp_path / 'migrate.sqlite').resolve()}"
    db = SubmissionDb(db_url=db_url, author=None)
    db.upgrade_schema(revision=PRE_CASES_REVISION)

    engine = sqlalchemy.create_engine(db_url)
    submissions = sqlalchemy.Table("submissions", sqlalchemy.MetaData(), autoload_with=engine)

    initial_a = "111111111_2025-01-01_00000001"
    initial_b = "111111111_2025-01-02_00000002"
    rows = [
        # two 'initial' submissions for the same (submitter, local case id) -> inconsistent
        dict(id=initial_a, tan_g=_tan(1), pseudonym="caseX", submitter_id="111111111", submission_type="initial"),
        dict(id=initial_b, tan_g=_tan(2), pseudonym="caseX", submitter_id="111111111", submission_type="initial"),
    ]
    with engine.begin() as conn:
        conn.execute(submissions.insert(), rows)

    with pytest.raises(RuntimeError, match="more than one 'initial' submission shares"):
        db.upgrade_schema(revision=CASES_REVISION)

    # the failed upgrade must not have created the cases table
    inspector = sqlalchemy.inspect(engine)
    assert "cases" not in inspector.get_table_names()
