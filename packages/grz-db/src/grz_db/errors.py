class SubmissionNotFoundError(ValueError):
    """Exception for when a submission is not found in the database."""

    def __init__(self, submission_id: str):
        super().__init__(f"Submission not found for ID {submission_id}")


class DuplicateSubmissionError(ValueError):
    """Exception for when a submission ID already exists in the database."""

    def __init__(self, submission_id: str):
        super().__init__(f"Duplicate submission ID {submission_id}")


class DuplicateTanGError(ValueError):
    """Exception for when a tanG is already in use."""

    def __init__(self):
        super().__init__("Duplicate tanG")


class DatabaseConfigurationError(Exception):
    """Exception for database configuration issues."""

    pass
