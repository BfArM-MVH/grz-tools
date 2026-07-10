"""Tests for the case tracking layer: the ``Case`` model, the case service methods, the pluggable
case resolver, and the one-initial-per-case rule enforced by ``SubmissionDb.assign_case``.
"""

from pathlib import Path

import pytest
from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives.asymmetric import ed25519
from grz_db.errors import (
    CaseHasLinkedSubmissionsError,
    CaseNotFoundError,
    DuplicatePsnError,
    SubmissionTypeInvalidForCaseError,
)
from grz_db.models.author import Author
from grz_db.models.submission import PsnResolver, SubmissionDb, SubmissionType

SUBMITTER_A = "111111111"
SUBMITTER_B = "222222222"


@pytest.fixture
def test_author() -> Author:
    key = ed25519.Ed25519PrivateKey.generate()
    private_key_bytes = key.private_bytes(
        encoding=serialization.Encoding.PEM,
        format=serialization.PrivateFormat.OpenSSH,
        encryption_algorithm=serialization.NoEncryption(),
    )
    return Author(name="alice", private_key_bytes=private_key_bytes, private_key_passphrase="")


@pytest.fixture
def db(tmp_path: Path, test_author: Author) -> SubmissionDb:
    db_url = f"sqlite:///{(tmp_path / 'submissions.sqlite').resolve()}"
    submission_db = SubmissionDb(db_url=db_url, author=test_author, debug=False)
    submission_db.initialize_schema()
    return submission_db


def _sid(submitter_id: str, suffix: str) -> str:
    return f"{submitter_id}_2025-01-01_{suffix}"


def _add(db: SubmissionDb, submission_id: str, submission_type: SubmissionType) -> None:
    db.add_submission(submission_id)
    db.modify_submission(submission_id, "submission_type", submission_type.value)


def test_create_case_psn_uniqueness_but_key_not_constrained(db: SubmissionDb):
    # (submitter_id, local_case_id) is NOT a uniqueness constraint: duplicates are allowed.
    db.create_case(SUBMITTER_A, "caseX")
    db.create_case(SUBMITTER_A, "caseX")

    # psn IS unique when set.
    db.create_case(SUBMITTER_A, "caseY", psn="PSN-1")
    with pytest.raises(DuplicatePsnError):
        db.create_case(SUBMITTER_B, "caseZ", psn="PSN-1")


def test_assign_case_initial_creates_and_links(db: SubmissionDb):
    sid = _sid(SUBMITTER_A, "0000000a")
    _add(db, sid, SubmissionType.initial)

    case = db.assign_case(sid, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.initial)

    assert db.get_submission(sid).case_id == case.id
    assert case.submitter_id == SUBMITTER_A and case.local_case_id == "caseX" and case.psn is None
    assert db.get_initial_submission(case.id).id == sid


def test_assign_case_addition_links_existing_case(db: SubmissionDb):
    initial = _sid(SUBMITTER_A, "0000000a")
    addition = _sid(SUBMITTER_A, "0000000b")
    _add(db, initial, SubmissionType.initial)
    _add(db, addition, SubmissionType.addition)

    case = db.assign_case(
        initial, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.initial
    )
    case2 = db.assign_case(
        addition, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.addition
    )

    assert case2.id == case.id
    assert {s.id for s in db.list_submissions_for_case(case.id)} == {initial, addition}


def test_assign_case_rejects_non_initial_on_new_case(db: SubmissionDb):
    sid = _sid(SUBMITTER_A, "0000000c")
    _add(db, sid, SubmissionType.correction)
    with pytest.raises(SubmissionTypeInvalidForCaseError):
        db.assign_case(sid, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.correction)


def test_assign_case_rejects_duplicate_initial(db: SubmissionDb):
    first = _sid(SUBMITTER_A, "0000000a")
    second = _sid(SUBMITTER_A, "0000000b")
    _add(db, first, SubmissionType.initial)
    _add(db, second, SubmissionType.initial)

    db.assign_case(first, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.initial)
    with pytest.raises(SubmissionTypeInvalidForCaseError):
        db.assign_case(second, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.initial)


def test_assign_case_idempotent(db: SubmissionDb):
    sid = _sid(SUBMITTER_A, "0000000a")
    _add(db, sid, SubmissionType.initial)

    case = db.assign_case(sid, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.initial)
    case_again = db.assign_case(
        sid, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.initial
    )

    assert case_again.id == case.id
    assert [s.id for s in db.list_submissions_for_case(case.id)] == [sid]


def test_two_submitters_same_local_case_id_are_distinct_cases(db: SubmissionDb):
    sid_a = _sid(SUBMITTER_A, "0000000a")
    sid_b = _sid(SUBMITTER_B, "0000000b")
    _add(db, sid_a, SubmissionType.initial)
    _add(db, sid_b, SubmissionType.initial)

    case_a = db.assign_case(
        sid_a, submitter_id=SUBMITTER_A, local_case_id="sharedLocalId", submission_type=SubmissionType.initial
    )
    case_b = db.assign_case(
        sid_b, submitter_id=SUBMITTER_B, local_case_id="sharedLocalId", submission_type=SubmissionType.initial
    )

    assert case_a.id != case_b.id


def test_assign_case_with_psn_resolver_links_by_pseudonym(db: SubmissionDb):
    """The pluggable resolver can locate a case by psn (the future PSN-based flow)."""
    pre = db.create_case(submitter_id=SUBMITTER_A, local_case_id="caseX", psn="PSN-1")
    sid = _sid(SUBMITTER_A, "0000000a")
    _add(db, sid, SubmissionType.initial)

    linked = db.assign_case(sid, psn="PSN-1", submission_type=SubmissionType.initial, resolver=PsnResolver())

    assert linked.id == pre.id
    assert db.get_submission(sid).case_id == pre.id


def test_modify_case_psn_duplicate_and_unknown_key(db: SubmissionDb):
    case = db.create_case(SUBMITTER_A, "caseX")
    db.modify_case(case.id, "psn", "RKI-PSN-1")
    assert db.get_case(case.id).psn == "RKI-PSN-1"

    other = db.create_case(SUBMITTER_B, "caseY")
    with pytest.raises(DuplicatePsnError):
        db.modify_case(other.id, "psn", "RKI-PSN-1")

    with pytest.raises(ValueError, match="Unknown column key"):
        db.modify_case(case.id, "not_a_column", "x")
    with pytest.raises(CaseNotFoundError):
        db.modify_case(9999, "psn", "x")


def test_delete_case_empty_and_refuses_when_linked(db: SubmissionDb):
    sid = _sid(SUBMITTER_A, "0000000a")
    _add(db, sid, SubmissionType.initial)
    case = db.assign_case(sid, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.initial)

    with pytest.raises(CaseHasLinkedSubmissionsError):
        db.delete_case(case.id)

    empty = db.create_case(SUBMITTER_A, "emptyCase")
    db.delete_case(empty.id)
    assert db.get_case(empty.id) is None


def test_relink_rejects_second_initial_via_db_constraint(db: SubmissionDb):
    """set_submission_case skips the application-level check; the partial unique index is the last line of defense."""
    sid_a = _sid(SUBMITTER_A, "0000000a")
    sid_b = _sid(SUBMITTER_B, "0000000b")
    _add(db, sid_a, SubmissionType.initial)
    _add(db, sid_b, SubmissionType.initial)

    case_a = db.assign_case(
        sid_a, submitter_id=SUBMITTER_A, local_case_id="caseA", submission_type=SubmissionType.initial
    )
    db.assign_case(sid_b, submitter_id=SUBMITTER_B, local_case_id="caseB", submission_type=SubmissionType.initial)

    with pytest.raises(SubmissionTypeInvalidForCaseError):
        db.set_submission_case(sid_b, case_a.id)
