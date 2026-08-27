"""Tests for the case tracking layer: the ``Case`` model, the case service methods, the pluggable
case resolver, and the rule limiting a case to one ``initial`` submission that passed basic QC.
"""

import json
import re
import threading
from collections.abc import Callable, Iterator
from contextlib import contextmanager

import pytest
import sqlalchemy
from grz_db.errors import (
    AmbiguousCaseError,
    CaseHasLinkedSubmissionsError,
    CaseNotFoundError,
    DuplicateCaseError,
    DuplicateInitialSubmissionError,
    DuplicatePsnError,
    SubmissionNotFoundError,
    SubmissionTypeInvalidForCaseError,
)
from grz_db.models.submission import PsnResolver, SubmissionDb, SubmissionType, _Index, _rejected_by
from grz_pydantic_models.submission.metadata import REDACTED_LOCAL_CASE_ID, GrzSubmissionMetadata
from pydantic import ValidationError
from sqlalchemy.exc import IntegrityError

SUBMITTER_A = "111111111"
SUBMITTER_B = "222222222"


def _with_submission_type(metadata: GrzSubmissionMetadata, submission_type: str) -> GrzSubmissionMetadata:
    """A copy of *metadata* with its ``submissionType`` replaced.

    The shipped example is a ``test`` submission, which is never case-tracked, so tests
    exercising case linking need an ``initial`` variant.
    """
    raw = json.loads(metadata.model_dump_json(by_alias=True))
    raw["submission"]["submissionType"] = submission_type
    return GrzSubmissionMetadata.model_validate(raw)


def _sid(submitter_id: str, suffix: str) -> str:
    return f"{submitter_id}_2025-01-01_{suffix}"


def _add(db: SubmissionDb, submission_id: str, submission_type: SubmissionType) -> None:
    db.add_submission(submission_id)
    db.modify_submission(submission_id, "submission_type", submission_type.value)


def _record_basic_qc(db: SubmissionDb, submission_id: str, passed: bool) -> None:
    db.modify_submission(submission_id, "basic_qc_passed", passed)


def _lose_the_race_to(monkeypatch: pytest.MonkeyPatch, win: Callable[[str, str], None]) -> None:
    """Hand the key to *win* in the window between looking a case up and inserting one.

    Threads cannot be made to interleave on demand, so the race is forced here instead. Only
    the first pass through the window hands it over, so *win* may itself resolve a case.
    """
    original = SubmissionDb._find_case_for_link
    won = False

    def let_the_other_process_win(self, session, submission_id, **kwargs):
        nonlocal won
        submission, case = original(self, session, submission_id, **kwargs)
        if case is None and not won:
            won = True
            win(kwargs["submitter_id"], kwargs["local_case_id"])
        return submission, case

    monkeypatch.setattr(SubmissionDb, "_find_case_for_link", let_the_other_process_win)


def _second_case_for_the_same_key(db: SubmissionDb, submitter_id: str, local_case_id: str) -> None:
    """Give one key a second case, which ``ux_cases_submitter_local_case`` forbids.

    Only reachable by writing around the application, which is exactly the state the ambiguity
    guard exists to survive: the index is dropped first so the row can be inserted at all.
    """
    with db.transaction() as session:
        session.execute(sqlalchemy.text("DROP INDEX ux_cases_submitter_local_case"))
        session.execute(
            sqlalchemy.text("INSERT INTO cases (submitter_id, local_case_id) VALUES (:submitter, :local_case)"),
            {"submitter": submitter_id, "local_case": local_case_id},
        )
        session.commit()


def test_create_case_key_and_psn_are_both_unique(db: SubmissionDb):
    # a submitter's local case ID identifies one patient, so it identifies one case
    db.create_case(SUBMITTER_A, "caseX")
    with pytest.raises(DuplicateCaseError):
        db.create_case(SUBMITTER_A, "caseX")

    # ... while the key stays unconstrained where either half is absent
    db.create_case(SUBMITTER_A, None)
    db.create_case(SUBMITTER_A, None)

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

    # get_initial_submission returns None until an initial submission of the case passes basic QC
    assert db.get_initial_submission(case.id) is None
    _record_basic_qc(db, sid, True)
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


@pytest.mark.parametrize(
    "submission_type", [SubmissionType.followup, SubmissionType.addition, SubmissionType.correction]
)
def test_non_initial_opens_a_case_of_its_own(db: SubmissionDb, submission_type: SubmissionType):
    """A patient's first submission to reach this GRZ need not be the initial one.

    The submission is real, it is reported to BfArM on its own tanG, and refusing to link it
    would only keep it out of the record.
    """
    sid = _sid(SUBMITTER_A, "0000000c")
    _add(db, sid, submission_type)

    assert (
        db.resolve_case(sid, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=submission_type) is None
    )
    case = db.assign_case(sid, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=submission_type)

    assert db.get_submission(sid).case_id == case.id
    assert db.get_initial_submission(case.id) is None  # the case simply has no initial submission


def test_a_later_initial_joins_the_case_its_followup_opened(db: SubmissionDb):
    followup = _sid(SUBMITTER_A, "0000000c")
    initial = _sid(SUBMITTER_A, "0000000a")
    _add(db, followup, SubmissionType.followup)
    _add(db, initial, SubmissionType.initial)

    opened = db.assign_case(
        followup, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.followup
    )
    joined = db.assign_case(
        initial, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.initial
    )

    assert joined.id == opened.id
    _record_basic_qc(db, initial, True)
    assert db.get_initial_submission(opened.id).id == initial


def test_second_initial_still_links_once_one_passed_qc(db: SubmissionDb):
    """A rival initial submission stays on record, linked to the case it duplicates. The rule
    binds when it tries to pass basic QC, not when it is linked.
    """
    qc_passed_id, case_id = _case_with_qc_passed_initial(db)
    second = _sid(SUBMITTER_A, "0000000b")
    _add(db, second, SubmissionType.initial)

    preview = db.resolve_case(
        second, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.initial
    )
    case = db.assign_case(
        second, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.initial
    )

    assert preview is not None and preview.id == case_id and case.id == case_id
    assert {s.id for s in db.list_submissions_for_case(case_id)} == {qc_passed_id, second}
    # the case keeps its QC-passed initial submission, and the rival may never become a second one
    assert db.get_initial_submission(case_id).id == qc_passed_id
    with pytest.raises(DuplicateInitialSubmissionError):
        _record_basic_qc(db, second, True)


def test_duplicate_initial_that_already_failed_qc_links(db: SubmissionDb):
    """The state ``grzctl validate`` leaves behind: basic QC already failed, nothing else known."""
    qc_passed_id, case_id = _case_with_qc_passed_initial(db)
    second = _sid(SUBMITTER_A, "0000000b")
    _add(db, second, SubmissionType.initial)
    _record_basic_qc(db, second, False)

    case = db.assign_case(
        second, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.initial
    )

    assert case.id == case_id
    assert db.get_submission(second).case_id == case_id
    assert db.get_initial_submission(case_id).id == qc_passed_id


def test_assign_case_allows_competing_initials_while_qc_pending(db: SubmissionDb):
    """The data alone cannot tell a re-upload from a duplicate, so competing initial submissions
    may coexist while basic QC is pending.
    """
    first = _sid(SUBMITTER_A, "0000000a")
    second = _sid(SUBMITTER_A, "0000000b")
    _add(db, first, SubmissionType.initial)
    _add(db, second, SubmissionType.initial)

    case = db.assign_case(
        first, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.initial
    )
    case2 = db.assign_case(
        second, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.initial
    )

    assert case2.id == case.id
    assert {s.id for s in db.list_submissions_for_case(case.id)} == {first, second}
    assert db.get_initial_submission(case.id) is None  # basic QC is still pending


def test_corrected_initial_replaces_failed_one_without_relinking(db: SubmissionDb):
    """A corrected re-upload after failed basic QC becomes the case's QC-passed initial submission
    with no operator steps.
    """
    bad = _sid(SUBMITTER_A, "0000000a")
    corrected = _sid(SUBMITTER_A, "0000000b")
    followup = _sid(SUBMITTER_A, "0000000c")
    _add(db, bad, SubmissionType.initial)
    _add(db, corrected, SubmissionType.initial)
    _add(db, followup, SubmissionType.followup)

    case = db.assign_case(bad, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.initial)
    _record_basic_qc(db, bad, False)

    case2 = db.assign_case(
        corrected, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.initial
    )
    _record_basic_qc(db, corrected, True)

    assert case2.id == case.id
    assert db.get_initial_submission(case.id).id == corrected
    # the failed upload stays linked for the record
    assert {s.id for s in db.list_submissions_for_case(case.id)} == {bad, corrected}

    case3 = db.assign_case(
        followup, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.followup
    )
    assert case3.id == case.id


def _pass_basic_qc_via_modify(db: SubmissionDb, submission_id: str) -> None:
    db.modify_submission(submission_id, "basic_qc_passed", True)


def _pass_basic_qc_via_update(db: SubmissionDb, submission_id: str) -> None:
    submission = db.get_submission(submission_id)
    submission.basic_qc_passed = True
    db.update_submission(submission)


@pytest.mark.parametrize(
    "record_pass",
    [_pass_basic_qc_via_modify, _pass_basic_qc_via_update],
    ids=["modify_submission", "update_submission"],
)
def test_second_initial_to_pass_basic_qc_is_rejected(
    db: SubmissionDb, record_pass: Callable[[SubmissionDb, str], None]
):
    """Once one initial submission has passed basic QC, recording a pass for a competing initial
    submission of the same case raises instead of writing.
    """
    qc_passed_id, duplicate_id, case_id = _case_with_competing_initials(db)

    _record_basic_qc(db, qc_passed_id, True)
    with pytest.raises(DuplicateInitialSubmissionError, match=qc_passed_id):
        record_pass(db, duplicate_id)

    # a duplicate initial submission can still be recorded as failed; failure is a terminal basic QC outcome
    _record_basic_qc(db, duplicate_id, False)
    assert db.get_submission(duplicate_id).basic_qc_passed is False
    assert db.get_initial_submission(case_id).id == qc_passed_id


def test_followup_links_even_when_the_only_initial_failed_qc(db: SubmissionDb):
    """Whether the case has a usable initial submission is a reimbursement question for BfArM,
    which reads the Prüfbericht, not a reason to keep the followup out of the database.
    """
    initial = _sid(SUBMITTER_A, "0000000a")
    followup = _sid(SUBMITTER_A, "0000000b")
    _add(db, initial, SubmissionType.initial)
    _add(db, followup, SubmissionType.followup)

    case = db.assign_case(
        initial, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.initial
    )
    _record_basic_qc(db, initial, False)

    joined = db.assign_case(
        followup, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.followup
    )

    assert joined.id == case.id
    assert {s.id for s in db.list_submissions_for_case(case.id)} == {initial, followup}
    assert db.get_initial_submission(case.id) is None


def _case_with_qc_passed_initial(db: SubmissionDb) -> tuple[str, int]:
    """A case whose initial submission has passed basic QC; returns (id of the case's QC-passed
    initial submission, case id).
    """
    qc_passed_id = _sid(SUBMITTER_A, "0000000a")
    _add(db, qc_passed_id, SubmissionType.initial)
    case = db.assign_case(
        qc_passed_id, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.initial
    )
    _record_basic_qc(db, qc_passed_id, True)
    return qc_passed_id, case.id


def _case_with_competing_initials(db: SubmissionDb) -> tuple[str, str, int]:
    """A case with two competing initial submissions linked to it, basic QC still pending for
    both; returns (id of the submission that will pass basic QC, id of the other initial
    submission, case id).
    """
    qc_passed_id = _sid(SUBMITTER_A, "0000000a")
    duplicate_id = _sid(SUBMITTER_A, "0000000b")
    _add(db, qc_passed_id, SubmissionType.initial)
    _add(db, duplicate_id, SubmissionType.initial)
    case = db.assign_case(
        qc_passed_id, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.initial
    )
    db.assign_case(
        duplicate_id, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.initial
    )
    return qc_passed_id, duplicate_id, case.id


def test_modify_submission_type_cannot_smuggle_in_a_second_initial(db: SubmissionDb):
    """Changing submission_type to 'initial' on an already QC-passed submission must not create a
    second initial submission that passed basic QC for the case.
    """
    _case_with_qc_passed_initial(db)
    addition = _sid(SUBMITTER_A, "0000000b")
    _add(db, addition, SubmissionType.addition)
    db.assign_case(addition, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.addition)
    _record_basic_qc(db, addition, True)  # fine: addition is not an initial submission

    with pytest.raises(DuplicateInitialSubmissionError):
        db.modify_submission(addition, "submission_type", SubmissionType.initial.value)
    assert db.get_submission(addition).submission_type == SubmissionType.addition


def test_one_initial_violation_is_recognized_in_the_backend_error(db: SubmissionDb):
    """The partial unique index enforcing the rule must compile on every supported backend.

    Raw SQL bypasses the ORM path, so the backend's own error reaches the test untranslated.
    """
    qc_passed_id, duplicate_id, _case_id = _case_with_competing_initials(db)
    _record_basic_qc(db, qc_passed_id, True)

    with pytest.raises(IntegrityError) as excinfo:
        with db.transaction() as session:
            session.execute(
                sqlalchemy.text("UPDATE submissions SET basic_qc_passed = :passed WHERE id = :sid"),
                {"passed": True, "sid": duplicate_id},
            )
            session.commit()

    assert _rejected_by(excinfo.value) is _Index.ONE_INITIAL


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


def test_assign_case_with_psn_resolver_links_by_pseudonym(db: SubmissionDb, migrated_db_connection, test_author):
    """The pluggable resolver can locate a case by psn (the future PSN-based flow)."""
    psn_db = SubmissionDb(db_url=migrated_db_connection, author=test_author, case_resolver=PsnResolver())
    pre = db.create_case(submitter_id=SUBMITTER_A, local_case_id="caseX", psn="PSN-1")
    sid = _sid(SUBMITTER_A, "0000000a")
    _add(db, sid, SubmissionType.initial)

    linked = psn_db.assign_case(sid, psn="PSN-1", submission_type=SubmissionType.initial)

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


def test_relink_rejects_second_qc_passed_initial(db: SubmissionDb):
    """Relinking must not give a case a second initial submission that passed basic QC; an initial
    submission whose basic QC is still pending may still be relinked freely.
    """
    _qc_passed_id, case_a_id = _case_with_qc_passed_initial(db)
    sid_b = _sid(SUBMITTER_B, "0000000b")
    _add(db, sid_b, SubmissionType.initial)
    case_b = db.assign_case(
        sid_b, submitter_id=SUBMITTER_B, local_case_id="caseB", submission_type=SubmissionType.initial
    )
    _record_basic_qc(db, sid_b, True)

    with pytest.raises(DuplicateInitialSubmissionError):
        db.set_submission_case(sid_b, case_a_id)

    # an initial submission whose basic QC is still pending is not yet a duplicate, so it may be relinked
    sid_c = _sid(SUBMITTER_A, "0000000c")
    _add(db, sid_c, SubmissionType.initial)
    db.set_submission_case(sid_c, case_a_id)
    assert db.get_submission(sid_c).case_id == case_a_id
    assert db.get_initial_submission(case_b.id).id == sid_b


def test_unlink_leaves_the_case_behind(db: SubmissionDb):
    """Unlinking frees the submission, not the case: the case stays, possibly empty."""
    sid = _sid(SUBMITTER_A, "0000001a")
    _add(db, sid, SubmissionType.initial)
    case = db.assign_case(sid, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.initial)

    unlinked = db.clear_submission_case(sid)
    assert unlinked.case_id is None
    assert case.id is not None
    assert db.get_case(case.id) is not None
    assert db.list_submissions_for_case(case.id) == []


def test_unlink_frees_the_one_initial_slot(db: SubmissionDb):
    """The repair this exists for: the wrong initial submission took a case's one QC-passed slot.

    Relinking cannot fix that, since it can only move the wrong submission to some other case.
    Unlinking it lets the right one take the slot.
    """
    wrong, right = _sid(SUBMITTER_A, "0000002a"), _sid(SUBMITTER_A, "0000002b")
    _add(db, wrong, SubmissionType.initial)
    case = db.assign_case(
        wrong, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.initial
    )
    _record_basic_qc(db, wrong, True)

    _add(db, right, SubmissionType.initial)
    assert case.id is not None
    db.set_submission_case(right, case.id)
    with pytest.raises(DuplicateInitialSubmissionError):
        _record_basic_qc(db, right, True)

    db.clear_submission_case(wrong)
    _record_basic_qc(db, right, True)
    holder = db.get_initial_submission(case.id)
    assert holder is not None and holder.id == right


def test_unlink_is_a_no_op_on_an_unlinked_submission(db: SubmissionDb):
    sid = _sid(SUBMITTER_A, "0000003a")
    _add(db, sid, SubmissionType.initial)
    assert db.clear_submission_case(sid).case_id is None


def test_unlink_rejects_an_unknown_submission(db: SubmissionDb):
    with pytest.raises(SubmissionNotFoundError):
        db.clear_submission_case(_sid(SUBMITTER_A, "0000004a"))


def test_resolve_case_previews_without_linking(db: SubmissionDb):
    sid = _sid(SUBMITTER_A, "000000aa")
    _add(db, sid, SubmissionType.initial)

    resolved = db.resolve_case(
        sid, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.initial
    )
    assert resolved is None
    # a second resolve still finds nothing, so the first one created no case
    assert (
        db.resolve_case(sid, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.initial)
        is None
    )
    submission = db.get_submission(sid)
    assert submission is not None and submission.case_id is None

    case = db.assign_case(sid, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.initial)
    resolved = db.resolve_case(
        sid, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.initial
    )
    assert resolved is not None and resolved.id == case.id


def test_diff_reports_pending_case_link_and_commit_applies_it(db: SubmissionDb, metadata: GrzSubmissionMetadata):
    initial_metadata = _with_submission_type(metadata, "initial")
    submission_id = initial_metadata.submission_id
    db.add_submission(submission_id)

    changes = db.diff(submission_id, initial_metadata, submission_uploaded_date=None)
    link = changes.case_link
    assert link is not None and link.before is None and link.after is None  # would create a new case
    submission = db.get_submission(submission_id)
    assert submission is not None and submission.case_id is None  # diff performed no writes

    db.commit_changes(submission_id, changes)
    submission = db.get_submission(submission_id)
    assert submission is not None and submission.case_id is not None

    # link is now in sync, so a fresh diff carries no pending case link
    changes = db.diff(submission_id, initial_metadata, submission_uploaded_date=None)
    assert changes.case_link is None


#: A read of the submissions table. SQLAlchemy renders SQL over several lines, so the
#: statement is collapsed to one before matching.
_SUBMISSION_READ = re.compile(r"^SELECT\b.*\bFROM submissions\b")


@contextmanager
def _counting_reads(engine: sqlalchemy.Engine) -> Iterator[dict[str, int]]:
    """Count the transactions opened and the submission rows read on *engine* inside the block."""
    counts = {"transactions": 0, "submission_reads": 0}

    def on_begin(_connection) -> None:
        counts["transactions"] += 1

    def on_statement(_connection, _cursor, statement: str, _parameters, _context, _executemany) -> None:
        if _SUBMISSION_READ.match(" ".join(statement.split())):
            counts["submission_reads"] += 1

    sqlalchemy.event.listen(engine, "begin", on_begin)
    sqlalchemy.event.listen(engine, "before_cursor_execute", on_statement)
    try:
        yield counts
    finally:
        sqlalchemy.event.remove(engine, "begin", on_begin)
        sqlalchemy.event.remove(engine, "before_cursor_execute", on_statement)


def test_diff_reads_one_snapshot(db: SubmissionDb, metadata: GrzSubmissionMetadata):
    """A change set must answer for one submission at one moment: it opens one transaction and
    reads the submission row once, so a write committed mid-diff cannot leave fields, donors
    and the case link describing three different snapshots.

    Two reads of the submissions table rather than one: the row itself, plus the aggregate
    behind :func:`case_key_denotes_one_patient`, which decides whether this key may open a case
    at all. That one counts rows instead of fetching the submission again, so the snapshot the
    change set is built from is still a single read.
    """
    initial_metadata = _with_submission_type(metadata, "initial")
    submission_id = initial_metadata.submission_id
    db.add_submission(submission_id)

    with _counting_reads(db.engine) as counts:
        changes = db.diff(submission_id, initial_metadata, submission_uploaded_date=None)

    assert changes.case_link is not None, "case resolution must have run, or this pins nothing"
    assert changes.donors.added, "donor resolution must have run, or this pins nothing"
    assert counts == {"transactions": 1, "submission_reads": 2}


def test_diff_resolves_by_psn_when_given_one(
    db: SubmissionDb, migrated_db_connection, test_author, metadata: GrzSubmissionMetadata
):
    """The whole seam, end to end: a psn-keyed deployment links through diff and commit_changes.

    The resolver is the deployment's, set once on the SubmissionDb, so diff and commit_changes
    cannot be given different ones. Only the psn is per submission, nothing in the metadata
    carrying one. The submitter's local case ID names a different case here, which must not be
    the one chosen.
    """
    psn_db = SubmissionDb(db_url=migrated_db_connection, author=test_author, case_resolver=PsnResolver())
    initial_metadata = _with_submission_type(metadata, "initial")
    decoy = db.create_case(initial_metadata.submission.submitter_id, initial_metadata.submission.local_case_id)
    wanted = db.create_case(SUBMITTER_B, "by-psn", psn="PSN-0001")

    submission_id = initial_metadata.submission_id
    db.add_submission(submission_id)
    changes = psn_db.diff(submission_id, initial_metadata, submission_uploaded_date=None, psn="PSN-0001")
    assert changes.case_link is not None
    assert changes.case_link.after == wanted.id != decoy.id
    assert changes.case_link.psn == "PSN-0001"

    psn_db.commit_changes(submission_id, changes)
    submission = db.get_submission(submission_id)
    assert submission is not None and submission.case_id == wanted.id


def test_diff_records_the_psn_on_a_case_it_opens(db: SubmissionDb, metadata: GrzSubmissionMetadata):
    """A psn passed to diff reaches the case commit_changes creates, not just its resolution."""
    initial_metadata = _with_submission_type(metadata, "initial")
    submission_id = initial_metadata.submission_id
    db.add_submission(submission_id)

    changes = db.diff(submission_id, initial_metadata, submission_uploaded_date=None, psn="PSN-0002")
    # the default resolver reads the local case ID and ignores the psn, so this opens a case
    assert changes.case_link is not None and changes.case_link.after is None

    db.commit_changes(submission_id, changes)
    submission = db.get_submission(submission_id)
    assert submission is not None and submission.case_id is not None
    opened = db.get_case(submission.case_id)
    assert opened is not None and opened.psn == "PSN-0002"


def test_diff_reports_relink_to_existing_case(db: SubmissionDb, metadata: GrzSubmissionMetadata):
    initial_metadata = _with_submission_type(metadata, "initial")
    submission_id = initial_metadata.submission_id
    db.add_submission(submission_id)
    db.commit_changes(submission_id, db.diff(submission_id, initial_metadata, submission_uploaded_date=None))
    submission = db.get_submission(submission_id)
    assert submission is not None and submission.case_id is not None
    original_case_id = submission.case_id

    other = db.create_case(initial_metadata.submission.submitter_id, "some-other-case")
    db.set_submission_case(submission_id, other.id)

    link = db.diff(submission_id, initial_metadata, submission_uploaded_date=None).case_link
    assert link is not None and link.before == other.id and link.after == original_case_id


def test_populate_does_not_revert_a_relink_without_force(db: SubmissionDb, metadata: GrzSubmissionMetadata):
    """Undoing a deliberate relink to another case is destructive, so populate must refuse it without force."""
    initial_metadata = _with_submission_type(metadata, "initial")
    submission_id = initial_metadata.submission_id
    db.add_submission(submission_id)
    db.commit_changes(submission_id, db.diff(submission_id, initial_metadata, submission_uploaded_date=None))
    original_case_id = db.get_submission(submission_id).case_id

    other = db.create_case(initial_metadata.submission.submitter_id, "some-other-case")
    db.set_submission_case(submission_id, other.id)

    with pytest.raises(RuntimeError, match="force"):
        db.populate(submission_id, initial_metadata, submission_date=None)
    assert db.get_submission(submission_id).case_id == other.id, "the relink must survive"

    db.populate(submission_id, initial_metadata, submission_date=None, force=True)
    assert db.get_submission(submission_id).case_id == original_case_id


def test_diff_records_a_duplicate_initial_against_its_case(db: SubmissionDb, metadata: GrzSubmissionMetadata):
    """A duplicate initial submission is recorded in full, linked to the case it duplicates,
    rather than nothing being written at all.
    """
    initial_metadata = _with_submission_type(metadata, "initial")
    submission_id = initial_metadata.submission_id
    db.add_submission(submission_id)
    db.commit_changes(submission_id, db.diff(submission_id, initial_metadata, submission_uploaded_date=None))
    db.modify_submission(submission_id, "basic_qc_passed", True)
    case_id = db.get_submission(submission_id).case_id

    second = _sid(initial_metadata.submission.submitter_id, "000000fe")
    db.add_submission(second)
    # the rival carries the same tanG in this fixture, which is a separate uniqueness rule
    changes = db.diff(second, initial_metadata, submission_uploaded_date=None, ignore_fields={"tan_g"})
    db.commit_changes(second, changes)

    row = db.get_submission(second)
    assert row.case_id == case_id
    assert row.submission_type == SubmissionType.initial
    assert row.submission_size is not None and row.submission_metadata is not None
    assert db.get_donors(second)
    # the case still has exactly one QC-passed initial submission
    assert db.get_initial_submission(case_id).id == submission_id


def test_resolution_rejects_ambiguous_case_key(db: SubmissionDb):
    # duplicate (submitter, local case id) pairs are storable, but resolving through them must
    # refuse rather than silently pick one patient's case
    case_a = db.create_case(SUBMITTER_A, "caseX")
    _second_case_for_the_same_key(db, SUBMITTER_A, "caseX")
    case_b = next(case for case, _n in db.list_cases() if case.id != case_a.id)
    sid = _sid(SUBMITTER_A, "000000b0")
    _add(db, sid, SubmissionType.initial)

    # the error names the offending cases so the operator knows what to clean up
    id_lo, id_hi = sorted([case_a.id, case_b.id])
    expected = rf"(?s)\b{id_lo}\b.*\b{id_hi}\b"
    with pytest.raises(AmbiguousCaseError, match=expected):
        db.resolve_case(sid, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.initial)
    with pytest.raises(AmbiguousCaseError, match=expected):
        db.assign_case(sid, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.initial)


def test_diff_skips_case_link_when_case_id_ignored(db: SubmissionDb, metadata: GrzSubmissionMetadata):
    initial_metadata = _with_submission_type(metadata, "initial")
    submission_id = initial_metadata.submission_id
    db.add_submission(submission_id)

    changes = db.diff(submission_id, initial_metadata, submission_uploaded_date=None, ignore_fields={"case_id"})
    assert changes.case_link is None


def test_diff_links_case_while_pseudonym_column_stays_ignored(db: SubmissionDb, metadata: GrzSubmissionMetadata):
    """Ignoring the pseudonym skips the column write but not case resolution."""
    initial_metadata = _with_submission_type(metadata, "initial")
    submission_id = initial_metadata.submission_id
    db.add_submission(submission_id)

    changes = db.diff(submission_id, initial_metadata, submission_uploaded_date=None, ignore_fields={"pseudonym"})
    assert changes.case_link is not None
    assert all(d.key != "pseudonym" for d in changes.fields.pending)

    db.commit_changes(submission_id, changes)
    submission = db.get_submission(submission_id)
    assert submission is not None and submission.case_id is not None
    assert submission.pseudonym is None


def test_linked_failed_initial_resolves_without_error(db: SubmissionDb):
    """A duplicate initial submission already linked to the case keeps its link when resolved
    again: resolve_case is a no-op for a submission already linked to the target case, not a
    fresh attempt to link it.
    """
    qc_passed_id, duplicate_id, _case_id = _case_with_competing_initials(db)
    _record_basic_qc(db, qc_passed_id, True)

    case = db.resolve_case(
        duplicate_id, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.initial
    )
    assert case is not None and db.get_submission(duplicate_id).case_id == case.id


def test_case_tracking_rejects_test_submission(db: SubmissionDb):
    sid = _sid(SUBMITTER_A, "000000b1")
    _add(db, sid, SubmissionType.test)

    with pytest.raises(SubmissionTypeInvalidForCaseError):
        db.resolve_case(sid, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.test)
    with pytest.raises(SubmissionTypeInvalidForCaseError):
        db.assign_case(sid, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.test)
    assert db.list_cases() == []


def test_diff_skips_case_link_for_test_submission(db: SubmissionDb, metadata: GrzSubmissionMetadata):
    # the shipped example is a test submission; populate must not case-track it
    submission_id = metadata.submission_id
    db.add_submission(submission_id)

    changes = db.diff(submission_id, metadata, submission_uploaded_date=None)
    assert changes.case_link is None

    db.commit_changes(submission_id, changes)
    submission = db.get_submission(submission_id)
    assert submission is not None and submission.case_id is None
    assert db.list_cases() == []


@pytest.mark.parametrize("placeholder", ["", REDACTED_LOCAL_CASE_ID], ids=["empty", "sentinel"])
def test_diff_skips_case_link_for_a_redacted_local_case_id(
    db: SubmissionDb, metadata: GrzSubmissionMetadata, placeholder: str
):
    """A redaction placeholder identifies no patient, so it must never key a case."""
    initial_metadata = _with_submission_type(metadata, "initial")
    raw = json.loads(initial_metadata.model_dump_json(by_alias=True))
    raw["submission"]["localCaseId"] = placeholder
    redacted = GrzSubmissionMetadata.model_validate(raw)

    submission_id = initial_metadata.submission_id
    db.add_submission(submission_id)

    assert db.diff(submission_id, redacted, submission_uploaded_date=None).case_link is None
    assert db.list_cases() == []


def test_two_submissions_with_a_redacted_local_case_id_do_not_merge(db: SubmissionDb, metadata: GrzSubmissionMetadata):
    """Without the guard, every redacted submission of a submitter shares the key ("", submitter)
    and collapses onto one case, so the second patient's initial is rejected as a duplicate.
    """
    initial_metadata = _with_submission_type(metadata, "initial")
    raw = json.loads(initial_metadata.model_dump_json(by_alias=True))
    raw["submission"]["localCaseId"] = ""
    redacted = GrzSubmissionMetadata.model_validate(raw)

    for suffix in ("000000d1", "000000d2"):
        sid = _sid(initial_metadata.submission.submitter_id, suffix)
        db.add_submission(sid)
        changes = db.diff(sid, redacted, submission_uploaded_date=None, ignore_fields={"tan_g"})
        assert changes.case_link is None
        db.commit_changes(sid, changes)
        _record_basic_qc(db, sid, True)

    assert db.list_cases() == []


def test_resolver_ignores_an_empty_case_key(db: SubmissionDb):
    """``db case create`` accepts an empty local case ID, but resolution must never match it."""
    db.create_case(SUBMITTER_A, "")
    sid = _sid(SUBMITTER_A, "000000d3")
    _add(db, sid, SubmissionType.initial)

    assert (
        db.resolve_case(sid, submitter_id=SUBMITTER_A, local_case_id="", submission_type=SubmissionType.initial) is None
    )


def _reuse_key_across_patients(
    db: SubmissionDb, submitter_id: str, local_case_id: str, suffixes: tuple[str, str]
) -> None:
    """Leave two QC-passed initial submissions sharing one key, which is what marks it reused.

    The state the cases migration leaves behind for a submitter that reused one local case ID:
    both rows unlinked, and both having passed basic QC on their own merits, since nothing was
    grouping them at the time.
    """
    for suffix in suffixes:
        sid = _sid(submitter_id, suffix)
        _add(db, sid, SubmissionType.initial)
        db.modify_submission(sid, "submitter_id", submitter_id)
        db.modify_submission(sid, "pseudonym", local_case_id)
        _record_basic_qc(db, sid, True)


def test_resolver_refuses_a_key_reused_across_patients(db: SubmissionDb):
    """A key carrying two QC-passed initial submissions denotes two patients, so it names no case.

    The case exists here, so this pins a refusal rather than a mere absence of a match: without
    the guard the newcomer would be linked to whichever patient reached the key first.
    """
    db.create_case(SUBMITTER_A, "ready")
    _reuse_key_across_patients(db, SUBMITTER_A, "ready", ("000000e1", "000000e2"))

    newcomer = _sid(SUBMITTER_A, "000000e3")
    _add(db, newcomer, SubmissionType.initial)
    assert (
        db.resolve_case(
            newcomer, submitter_id=SUBMITTER_A, local_case_id="ready", submission_type=SubmissionType.initial
        )
        is None
    )


def test_a_retry_that_failed_basic_qc_does_not_make_a_key_untrusted(db: SubmissionDb):
    """One patient may have several initial submissions; only one of them passes basic QC.

    Counting every initial submission instead would read a submitter who fixed a failed upload as
    one who reused the key, and leave a perfectly good case unlinked.
    """
    case = db.create_case(SUBMITTER_A, "caseX")
    for suffix, qc_passed in (("000000f1", False), ("000000f2", True)):
        sid = _sid(SUBMITTER_A, suffix)
        _add(db, sid, SubmissionType.initial)
        db.modify_submission(sid, "submitter_id", SUBMITTER_A)
        db.modify_submission(sid, "pseudonym", "caseX")
        _record_basic_qc(db, sid, qc_passed)

    followup = _sid(SUBMITTER_A, "000000f3")
    _add(db, followup, SubmissionType.followup)
    resolved = db.resolve_case(
        followup, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.followup
    )
    assert resolved is not None and resolved.id == case.id


def test_diff_leaves_a_key_reused_across_patients_unlinked(db: SubmissionDb, metadata: GrzSubmissionMetadata):
    """Populate must neither link to a case on such a key nor create one from it.

    Reporting "no match" instead would have ``commit_changes`` mint a case from the very key that
    names no patient, which is the shape this guard exists to avoid.
    """
    initial_metadata = _with_submission_type(metadata, "initial")
    _reuse_key_across_patients(
        db,
        initial_metadata.submission.submitter_id,
        initial_metadata.submission.local_case_id,
        ("000000e4", "000000e5"),
    )

    submission_id = initial_metadata.submission_id
    db.add_submission(submission_id)
    changes = db.diff(submission_id, initial_metadata, submission_uploaded_date=None)
    assert changes.case_link is None and changes.case_link_error is None

    db.commit_changes(submission_id, changes)
    submission = db.get_submission(submission_id)
    assert submission is not None and submission.case_id is None
    assert db.list_cases() == []


def test_no_duplicate_initial_is_reported_for_a_key_reused_across_patients(db: SubmissionDb):
    """The rejection this guard prevents: one patient failing basic QC for an unrelated one's sake.

    ``assert_no_duplicate_initial`` resolves through the same resolver, so a key that names no
    case cannot answer the duplicate question either.
    """
    linked = _sid(SUBMITTER_A, "0000000a")
    _add(db, linked, SubmissionType.initial)
    case = db.assign_case(
        linked, submitter_id=SUBMITTER_A, local_case_id="ready", submission_type=SubmissionType.initial
    )
    for sid in (linked, _sid(SUBMITTER_A, "0000000b")):
        if sid != linked:
            _add(db, sid, SubmissionType.initial)
        db.modify_submission(sid, "submitter_id", SUBMITTER_A)
        db.modify_submission(sid, "pseudonym", "ready")
        _record_basic_qc(db, sid, True)

    # without the guard the newcomer would resolve to this case and be failed for its sake
    assert case.id is not None
    holder = db.get_initial_submission(case.id)
    assert holder is not None and holder.id == linked

    newcomer = _sid(SUBMITTER_A, "0000000c")
    _add(db, newcomer, SubmissionType.initial)
    db.assert_no_duplicate_initial(
        newcomer, submitter_id=SUBMITTER_A, local_case_id="ready", submission_type=SubmissionType.initial
    )


def test_diff_reports_an_ambiguous_case_key_instead_of_raising(db: SubmissionDb, metadata: GrzSubmissionMetadata):
    """Whether a case can be identified says nothing about the other diffs, so diff reports it.

    ``resolve_case`` still raises, since a caller that asked only about the case has no answer.
    """
    initial_metadata = _with_submission_type(metadata, "initial")
    submitter = initial_metadata.submission.submitter_id
    local_case_id = initial_metadata.submission.local_case_id
    db.create_case(submitter, local_case_id)
    _second_case_for_the_same_key(db, submitter, local_case_id)
    submission_id = initial_metadata.submission_id
    db.add_submission(submission_id)

    changes = db.diff(submission_id, initial_metadata, submission_uploaded_date=None)

    assert isinstance(changes.case_link_error, AmbiguousCaseError)
    assert changes.case_link is None
    assert changes.has_pending, "every other diff survives an unresolvable case"

    db.commit_changes(submission_id, changes)
    row = db.get_submission(submission_id)
    assert row.case_id is None
    assert row.submission_size is not None and row.submission_metadata is not None


def test_create_case_rejects_a_malformed_submitter_id(db: SubmissionDb):
    """A submitter ID no submission can carry would leave the case permanently unresolvable."""
    with pytest.raises(ValidationError):
        db.create_case("not-a-submitter-id", "caseX")
    assert db.list_cases() == []


def test_modify_case_rejects_a_malformed_submitter_id(db: SubmissionDb):
    case = db.create_case(SUBMITTER_A, "caseX")

    with pytest.raises(ValidationError):
        db.modify_case(case.id, "submitter_id", "nope")

    stored = db.get_case(case.id)
    assert stored.submitter_id == SUBMITTER_A


def test_relink_refuses_a_test_submission(db: SubmissionDb):
    """Relinking is the one door an operator drives by hand, so it needs the rule most."""
    initial = _sid(SUBMITTER_A, "0000000a")
    test_submission = _sid(SUBMITTER_A, "0000000e")
    _add(db, initial, SubmissionType.initial)
    _add(db, test_submission, SubmissionType.test)
    case = db.assign_case(
        initial, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.initial
    )

    with pytest.raises(SubmissionTypeInvalidForCaseError):
        db.set_submission_case(test_submission, case.id)

    assert db.get_submission(test_submission).case_id is None
    assert {s.id for s in db.list_submissions_for_case(case.id)} == {initial}


def test_relink_refuses_a_submission_whose_type_is_unknown(db: SubmissionDb):
    """An unpopulated row may yet turn out to be a test submission.

    Linking it would decide that before the metadata does, and populate resolves the link
    itself, so there is nothing to gain by allowing it.
    """
    initial = _sid(SUBMITTER_A, "0000000a")
    _add(db, initial, SubmissionType.initial)
    case = db.assign_case(
        initial, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.initial
    )
    unpopulated = _sid(SUBMITTER_A, "0000000e")
    db.add_submission(unpopulated)

    with pytest.raises(SubmissionTypeInvalidForCaseError):
        db.set_submission_case(unpopulated, case.id)

    assert db.get_submission(unpopulated).case_id is None


def test_concurrent_assign_case_settles_on_one_case(db: SubmissionDb, migrated_db_connection, test_author):
    """Two processes populating one patient's submissions must not mint two cases.

    Both look the case up before either commits, so both try to create one. The loser joins the
    winner's case instead of leaving the key ambiguous for good.
    """
    sids = [_sid(SUBMITTER_A, f"0000000{suffix}") for suffix in "ab"]
    for sid in sids:
        _add(db, sid, SubmissionType.initial)

    services = [SubmissionDb(db_url=migrated_db_connection, author=test_author, debug=False) for _ in sids]
    start = threading.Barrier(len(sids))
    outcomes: dict[str, object] = {}

    def link(service: SubmissionDb, sid: str) -> None:
        start.wait(timeout=10)
        try:
            outcomes[sid] = service.assign_case(
                sid, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.initial
            ).id
        except Exception as exc:
            outcomes[sid] = f"{type(exc).__name__}: {exc}"

    threads = [threading.Thread(target=link, args=(service, sid)) for service, sid in zip(services, sids, strict=True)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join(timeout=30)
    for service in services:
        service.engine.dispose()

    assert len(db.list_cases()) == 1, f"one patient, one case: {outcomes}"
    case_id = db.list_cases()[0][0].id
    assert set(outcomes.values()) == {case_id}, outcomes
    assert {s.id for s in db.list_submissions_for_case(case_id)} == set(sids)


def test_assign_case_joins_a_case_that_appeared_mid_flight(
    db: SubmissionDb, migrated_db_connection, test_author, monkeypatch: pytest.MonkeyPatch
):
    """``test_concurrent_assign_case_settles_on_one_case`` cannot guarantee the two actually
    interleave, so force it.
    """
    sid = _sid(SUBMITTER_A, "0000000a")
    _add(db, sid, SubmissionType.initial)
    other = SubmissionDb(db_url=migrated_db_connection, author=test_author, debug=False)
    _lose_the_race_to(monkeypatch, other.create_case)
    try:
        case = db.assign_case(
            sid, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.initial
        )
    finally:
        other.engine.dispose()

    assert [c.id for c, _n in db.list_cases()] == [case.id], "the loser must join, not create a second case"
    assert db.get_submission(sid).case_id == case.id


def test_modify_case_names_the_index_that_rejected_the_write(db: SubmissionDb):
    """The key being edited says what was attempted; the index says what happened."""
    kept = db.create_case(SUBMITTER_A, "caseX", psn="PSN-1")
    other = db.create_case(SUBMITTER_A, "caseY", psn="PSN-2")

    with pytest.raises(DuplicatePsnError):
        db.modify_case(other.id, "psn", "PSN-1")
    with pytest.raises(DuplicateCaseError):
        db.modify_case(other.id, "local_case_id", "caseX")

    assert db.get_case(other.id).psn == "PSN-2"
    assert db.get_case(other.id).local_case_id == "caseY"
    assert db.get_case(kept.id).local_case_id == "caseX"


def test_every_unique_index_is_classified(db: SubmissionDb):
    """The classifier matches names the migrations own, so drift would otherwise be silent.

    A renamed index simply stops being recognised, and a rejected write surfaces as a raw
    IntegrityError instead of the error saying what happened. A new unique index nobody
    classified fails here too, which is the moment to decide what it should say.
    """
    inspector = sqlalchemy.inspect(db.engine)
    in_the_schema = {
        index["name"]: ", ".join(f"{table}.{column}" for column in index["column_names"])
        for table in ("cases", "submissions")
        for index in inspector.get_indexes(table)
        if index.get("unique")
    }

    assert {index.index_name: index.sqlite_columns for index in _Index} == in_the_schema


def test_a_rejected_change_set_writes_nothing(db: SubmissionDb, metadata: GrzSubmissionMetadata):
    """A change set is one answer to one metadata file, so it is applied whole or not at all:
    the case link and every field write share one transaction, and a rejection undoes all of
    it.
    """
    initial_metadata = _with_submission_type(metadata, "initial")
    submitter = initial_metadata.submission.submitter_id
    local_case_id = initial_metadata.submission.local_case_id

    winner = _sid(submitter, "0000000a")
    _add(db, winner, SubmissionType.initial)
    db.assign_case(winner, submitter_id=submitter, local_case_id=local_case_id, submission_type=SubmissionType.initial)
    _record_basic_qc(db, winner, True)

    # a rival that already passed basic QC, but whose type is not recorded yet: it links, and
    # only writing submission_type brings it under the one-initial index
    rival = _sid(submitter, "0000000b")
    db.add_submission(rival)
    _record_basic_qc(db, rival, True)

    changes = db.diff(rival, initial_metadata, submission_uploaded_date=None, ignore_fields={"tan_g"})
    assert changes.case_link is not None and any(d.key == "submission_type" for d in changes.fields.pending)

    with pytest.raises(DuplicateInitialSubmissionError):
        db.commit_changes(rival, changes)

    row = db.get_submission(rival)
    assert row.case_id is None, "the case link must not survive a rejected change set"
    assert row.submission_type is None and row.submission_size is None
    assert db.get_donors(rival) == ()


def test_a_change_set_survives_a_case_that_appeared_mid_flight(
    db: SubmissionDb,
    migrated_db_connection,
    test_author,
    metadata: GrzSubmissionMetadata,
    monkeypatch: pytest.MonkeyPatch,
):
    """Joining a case that appeared mid-flight rolls the change set's transaction back.

    That is survivable only because the link is applied before anything else: anything written
    ahead of it would be discarded by that rollback and the change set would commit anyway,
    which is the partial write the ordering exists to prevent.
    """
    initial_metadata = _with_submission_type(metadata, "initial")
    submitter = initial_metadata.submission.submitter_id

    sid = _sid(submitter, "0000000a")
    db.add_submission(sid)
    changes = db.diff(sid, initial_metadata, submission_uploaded_date=None, ignore_fields={"tan_g"})
    assert changes.case_link is not None and changes.donors.added, "the rollback needs something to lose"

    other = SubmissionDb(db_url=migrated_db_connection, author=test_author, debug=False)
    _lose_the_race_to(monkeypatch, other.create_case)
    try:
        db.commit_changes(sid, changes)
    finally:
        other.engine.dispose()

    row = db.get_submission(sid)
    assert [c.id for c, _n in db.list_cases()] == [row.case_id], "the loser must join, not create a second case"
    assert row.submission_type == SubmissionType.initial, "the writes after the link must survive"
    assert db.get_donors(sid) != ()


def test_joining_a_full_case_names_the_index_that_rejected_the_link(
    db: SubmissionDb, migrated_db_connection, test_author, monkeypatch: pytest.MonkeyPatch
):
    """Losing the race to create a case does not mean the case will take the submission.

    The winner links its own QC-passed initial alongside the case it creates, so the case is
    already full when the loser joins it, and that rejection has a name of its own.
    """
    winner = _sid(SUBMITTER_A, "0000000a")
    _add(db, winner, SubmissionType.initial)
    _record_basic_qc(db, winner, True)

    loser = _sid(SUBMITTER_A, "0000000b")
    _add(db, loser, SubmissionType.initial)
    _record_basic_qc(db, loser, True)

    other = SubmissionDb(db_url=migrated_db_connection, author=test_author, debug=False)
    _lose_the_race_to(
        monkeypatch,
        lambda submitter_id, local_case_id: other.assign_case(
            winner,
            submitter_id=submitter_id,
            local_case_id=local_case_id,
            submission_type=SubmissionType.initial,
        ),
    )
    try:
        with pytest.raises(DuplicateInitialSubmissionError):
            db.assign_case(
                loser, submitter_id=SUBMITTER_A, local_case_id="caseX", submission_type=SubmissionType.initial
            )
    finally:
        other.engine.dispose()

    assert db.get_submission(loser).case_id is None, "a rejected link must not be stored"


def test_assign_case_names_the_psn_that_rejected_the_case_it_opened(db: SubmissionDb):
    """Opening a case carries the pseudonym, so the psn index can reject the write that opens it.

    Only the insert can break that index; the link that follows writes ``submissions.case_id``
    and nothing a case is keyed on.
    """
    db.create_case(SUBMITTER_B, "caseY", psn="RKI-PSN-1")

    sid = _sid(SUBMITTER_A, "0000000a")
    _add(db, sid, SubmissionType.initial)

    with pytest.raises(DuplicatePsnError):
        db.assign_case(
            sid,
            submitter_id=SUBMITTER_A,
            local_case_id="caseX",
            psn="RKI-PSN-1",
            submission_type=SubmissionType.initial,
        )

    assert [c.local_case_id for c, _n in db.list_cases()] == ["caseY"], "a rejected case must not be stored"
    assert db.get_submission(sid).case_id is None, "and must not leave the submission linked"


def test_assert_no_duplicate_initial_never_keys_on_a_redaction_placeholder(db: SubmissionDb, metadata):
    """A redaction placeholder keys no case, so it cannot collapse two submissions onto one.

    The placeholder-keyed case is written directly, which no application path does, so a resolver
    that stopped rejecting placeholders would match it here and report a duplicate.
    """
    initial_metadata = _with_submission_type(metadata, "initial")
    submitter = initial_metadata.submission.submitter_id

    stray = db.create_case(submitter_id=submitter, local_case_id=REDACTED_LOCAL_CASE_ID)
    holder = _sid(submitter, "0000ee01")
    _add(db, holder, SubmissionType.initial)
    db.set_submission_case(holder, stray.id)
    _record_basic_qc(db, holder, True)

    rival = _sid(submitter, "0000ee02")
    _add(db, rival, SubmissionType.initial)

    db.assert_no_duplicate_initial(
        rival,
        submitter_id=submitter,
        local_case_id=REDACTED_LOCAL_CASE_ID,
        submission_type=SubmissionType.initial,
    )


def test_assert_no_duplicate_initial_consults_the_configured_resolver(
    db: SubmissionDb, migrated_db_connection, test_author, metadata
):
    """The db's own resolver answers here, on identifiers the default resolver cannot key.

    ``local_case_id`` is absent, so the default resolver has nothing to match; only a psn-keyed
    one can see the case. Screening on the default key before consulting the configured resolver
    would report this rival as clear.
    """
    psn_db = SubmissionDb(db_url=migrated_db_connection, author=test_author, case_resolver=PsnResolver())
    initial_metadata = _with_submission_type(metadata, "initial")
    submitter = initial_metadata.submission.submitter_id

    case = db.create_case(submitter_id=submitter, local_case_id="caseZ", psn="RKI-PSN-9")
    holder = _sid(submitter, "0000ff01")
    _add(db, holder, SubmissionType.initial)
    db.set_submission_case(holder, case.id)
    _record_basic_qc(db, holder, True)

    rival = _sid(submitter, "0000ff02")
    _add(db, rival, SubmissionType.initial)

    with pytest.raises(DuplicateInitialSubmissionError) as excinfo:
        psn_db.assert_no_duplicate_initial(rival, psn="RKI-PSN-9", submission_type=SubmissionType.initial)
    assert holder in str(excinfo.value)


def test_assert_no_duplicate_initial_names_the_submission_holding_the_slot(db: SubmissionDb, metadata):
    """The rival is told which submission already holds the case's one QC-passed initial slot.

    The holder is not a duplicate of itself, which is what lets a submission be re-validated.
    """
    initial_metadata = _with_submission_type(metadata, "initial")
    submitter = initial_metadata.submission.submitter_id
    local_case_id = initial_metadata.submission.local_case_id

    holder = _sid(submitter, "0000dd01")
    _add(db, holder, SubmissionType.initial)
    db.assign_case(holder, submitter_id=submitter, local_case_id=local_case_id, submission_type=SubmissionType.initial)
    _record_basic_qc(db, holder, True)

    db.assert_no_duplicate_initial(
        holder, submitter_id=submitter, local_case_id=local_case_id, submission_type=SubmissionType.initial
    )

    rival = _sid(submitter, "0000dd02")
    _add(db, rival, SubmissionType.initial)
    with pytest.raises(DuplicateInitialSubmissionError) as excinfo:
        db.assert_no_duplicate_initial(
            rival, submitter_id=submitter, local_case_id=local_case_id, submission_type=SubmissionType.initial
        )
    assert holder in str(excinfo.value)


def test_assert_no_duplicate_initial_does_not_swallow_a_resolution_failure(db: SubmissionDb, metadata):
    """Answering "not a duplicate" when the question could not be answered is how one gets through.

    Nothing downstream treats an unresolvable case as an error, so a swallow here would be
    the last chance to notice.
    """
    initial_metadata = _with_submission_type(metadata, "initial")
    submitter = initial_metadata.submission.submitter_id
    local_case_id = initial_metadata.submission.local_case_id
    db.create_case(submitter, local_case_id)
    _second_case_for_the_same_key(db, submitter, local_case_id)
    submission_id = initial_metadata.submission_id
    _add(db, submission_id, SubmissionType.initial)

    with pytest.raises(AmbiguousCaseError):
        db.assert_no_duplicate_initial(
            submission_id,
            submitter_id=submitter,
            local_case_id=local_case_id,
            submission_type=SubmissionType.initial,
        )
