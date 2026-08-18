from collections.abc import Sequence


class GrzDbError(Exception):
    """Base class for all grz_db errors."""


class SubmissionError(GrzDbError):
    """Base class for errors related to submissions."""


class DatabaseError(GrzDbError):
    """Base class for database-related errors."""


class SubmissionNotFoundError(SubmissionError):
    """Exception for when a submission is not found in the database."""

    def __init__(self, submission_id: str):
        super().__init__(f"Submission not found for ID {submission_id}")


class DuplicateSubmissionError(SubmissionError):
    """Exception for when a submission ID already exists in the database."""

    def __init__(self, submission_id: str):
        super().__init__(f"Duplicate submission ID {submission_id}")


class DuplicateTanGError(SubmissionError):
    """Exception for when a tanG is already in use."""

    def __init__(self):
        super().__init__("Duplicate tanG")


class SubmissionDateIsNoneError(SubmissionError):
    """Exception for when a submission date is None."""

    def __init__(self):
        super().__init__("Submission date is None")


class SubmissionTypeIsNoneError(SubmissionError):
    """Exception for when a submission type is None."""

    def __init__(self):
        super().__init__("Submission type is None")


class SubmissionBasicQCNotPassedError(SubmissionError):
    """Exception for when a submission has not passed basic QC."""

    def __init__(self, submission_id: str):
        super().__init__(f"Submission with ID {submission_id} has not passed basic QC")


class DatabaseConfigurationError(DatabaseError):
    """Exception for database configuration issues."""


class CaseError(GrzDbError):
    """Base class for errors related to cases."""


class CaseNotFoundError(CaseError):
    """Exception for when a case is not found in the database."""

    def __init__(self, case_id: int):
        super().__init__(f"Case not found for ID {case_id}")


class DuplicatePsnError(CaseError):
    """Exception for when a case with the given RKI pseudonym already exists."""

    def __init__(self, psn: str):
        super().__init__(f"A case with pseudonym '{psn}' already exists")


class DuplicateCaseError(CaseError):
    """Exception for when a case with the given submitter and local case ID already exists."""

    def __init__(self, submitter_id: str | None, local_case_id: str | None):
        super().__init__(f"A case for submitter '{submitter_id}' and local case '{local_case_id}' already exists")


class CaseHasLinkedSubmissionsError(CaseError):
    """Exception for when a case cannot be deleted because submissions are still linked to it."""

    def __init__(self, case_id: int, count: int):
        super().__init__(f"Case {case_id} still has {count} linked submission(s) and cannot be deleted")


class SubmissionTypeInvalidForCaseError(SubmissionError):
    """Exception for when a submission's type is incompatible with its case state."""


class DuplicateInitialSubmissionError(SubmissionTypeInvalidForCaseError):
    """Exception for when a case would end up with a second initial submission that passed basic QC."""

    def __init__(self, case_id: int | None, winning_submission_id: str | None = None):
        qc_passed = f" '{winning_submission_id}'" if winning_submission_id else ""
        super().__init__(
            f"Case {case_id} already has an initial submission that passed basic QC{qc_passed}; "
            "at most one initial submission per case may pass basic QC."
        )


class AmbiguousCaseError(SubmissionError):
    """Exception for when a case resolution key matches more than one case."""

    def __init__(self, submitter_id: str | None, local_case_id: str | None, case_ids: Sequence[int]):
        ids = ", ".join(str(case_id) for case_id in sorted(case_ids))
        super().__init__(
            f"Cases [{ids}] share (submitter '{submitter_id}', local case '{local_case_id}'); cannot resolve "
            "automatically. Relink the cases' submissions and delete the duplicate cases to disambiguate."
        )
