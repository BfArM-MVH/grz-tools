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


class CaseHasLinkedSubmissionsError(CaseError):
    """Exception for when a case cannot be deleted because submissions are still linked to it."""

    def __init__(self, case_id: int, count: int):
        super().__init__(f"Case {case_id} still has {count} linked submission(s) and cannot be deleted")


class SubmissionTypeInvalidForCaseError(SubmissionError):
    """Exception for when a submission's type is incompatible with its case state."""

    def __init__(self, message: str):
        super().__init__(message)
