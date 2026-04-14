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
