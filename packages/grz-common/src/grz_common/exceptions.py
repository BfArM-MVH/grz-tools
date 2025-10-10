class VersionFileError(Exception):
    """Base exception for version file related errors."""


class VersionFileNotFoundError(VersionFileError):
    """Raised when the version file is missing from S3."""


class VersionFileAccessError(VersionFileError):
    """Raised when the version file cannot be accessed due to permissions or network issues."""


class VersionFileValidationError(VersionFileError):
    """Raised when the version file content is invalid or cannot be parsed."""

