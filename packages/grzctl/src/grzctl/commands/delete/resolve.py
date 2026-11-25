"""Find the submissions a BfArM deletion request names, and what deleting them involves."""

from dataclasses import dataclass
from enum import StrEnum

import rich.table
from grz_db.models.submission import Submission, SubmissionDb

from ...models.deletion_request import DeletionRequest


class DeletionOutcome(StrEnum):
    """What ``grzctl delete run`` does for a tanG of a request."""

    DELETE = "will be deleted"
    RESUME = "deletion interrupted, will be finished"
    ALREADY_DELETED = "already deleted"
    NOT_PRESENT = "not present"


@dataclass(frozen=True)
class ResolvedTanG:
    """A tanG of a request, the submission that has it, and what happens to that submission.

    :param tan_g: The tanG as the request names it.
    :param submission: The submission with that tanG, ``None`` if the database has none.
    :param outcome: What ``grzctl delete run`` does for it.
    """

    tan_g: str
    submission: Submission | None
    outcome: DeletionOutcome


@dataclass(frozen=True)
class ResolvedDeletionRequest:
    """A deletion request with the submission of each of its tanGs.

    :param tan_gs: Each tanG of the request, in the request's order.
    :param case_warnings: One line for each submission whose case has submissions the request
        does not name.
    """

    tan_gs: list[ResolvedTanG]
    case_warnings: list[str]

    @property
    def to_delete(self) -> list[Submission]:
        """The submissions whose deletion ``grzctl delete run`` starts or finishes."""
        return [
            resolved.submission
            for resolved in self.tan_gs
            if resolved.submission is not None and resolved.outcome in (DeletionOutcome.DELETE, DeletionOutcome.RESUME)
        ]


class DeletionRequestRefusedError(Exception):
    """A deletion request that must not be carried out at all.

    :param problems: Why, one line each.
    """

    def __init__(self, problems: list[str]):
        super().__init__("\n".join(problems))
        self.problems = problems


_SUBMISSION_ID_LENGTH = len("000000000_2025-01-01_00000000")


def new_table(*columns: str) -> rich.table.Table:
    """A table with *columns* that never cuts or wraps a submission ID, which is copied from it."""
    table = rich.table.Table(box=None)
    for column in columns:
        if column == "Submission":
            table.add_column(column, no_wrap=True, min_width=_SUBMISSION_ID_LENGTH)
        else:
            table.add_column(column)
    return table


def short_tan_g(tan_g: str) -> str:
    """A tanG cut to its first and last four characters, to show it without writing it out."""
    return f"{tan_g[:4]}…{tan_g[-4:]}"


def _outcome(submission: Submission | None) -> DeletionOutcome:
    if submission is None:
        return DeletionOutcome.NOT_PRESENT
    if submission.is_deleted:
        return DeletionOutcome.ALREADY_DELETED
    if submission.deletion_started:
        return DeletionOutcome.RESUME
    return DeletionOutcome.DELETE


def resolve_deletion_request(request: DeletionRequest, db: SubmissionDb, grz_id: str) -> ResolvedDeletionRequest:
    """Find the submission of each tanG of a deletion request and check that the request may be carried out.

    :param request: The checked deletion request.
    :param db: The submission database.
    :param grz_id: ID of this GRZ.
    :returns: The submission and outcome of each tanG, and the case warnings.
    :raises DeletionRequestRefusedError: if the request is addressed to another GRZ, or a submission
        it names belongs to another LE or another GRZ.
    """
    problems: list[str] = []
    if request.data_node_id != grz_id:
        problems.append(f"The request is addressed to {request.data_node_id}, but this is {grz_id}.")

    resolved = []
    for tan_g in request.tan_g:
        submission = db.get_submission_by_tan_g(tan_g)
        resolved.append(ResolvedTanG(tan_g=tan_g, submission=submission, outcome=_outcome(submission)))
        if submission is None:
            continue
        if submission.submitter_id != request.submitter_id:
            problems.append(
                f"tanG {short_tan_g(tan_g)}: submission {submission.id} was submitted by "
                f"{submission.submitter_id}, but the request is from {request.submitter_id}."
            )
        if submission.data_node_id is not None and submission.data_node_id != grz_id:
            problems.append(
                f"tanG {short_tan_g(tan_g)}: submission {submission.id} names {submission.data_node_id} as its GRZ, "
                f"but this is {grz_id}."
            )
    if problems:
        raise DeletionRequestRefusedError(problems)

    named = {r.submission.id for r in resolved if r.submission is not None}
    case_warnings = []
    for r in resolved:
        if r.submission is None or r.outcome == DeletionOutcome.ALREADY_DELETED or r.submission.case_id is None:
            continue
        unnamed = [s for s in db.list_submissions_for_case(r.submission.case_id) if s.id not in named]
        if unnamed:
            others = ", ".join(f"{s.id} ({s.submission_type})" for s in unnamed)
            case_warnings.append(
                f"{r.submission.id}: case {r.submission.case_id} also has {others}, which the request does not name."
            )
    return ResolvedDeletionRequest(tan_gs=resolved, case_warnings=case_warnings)
