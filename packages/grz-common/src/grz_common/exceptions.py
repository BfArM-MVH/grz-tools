class VersionFileError(Exception):
    """Base exception for version file related errors."""


class VersionFileNotFoundError(VersionFileError):
    """Raised when the version file is missing from S3."""


class VersionFileAccessError(VersionFileError):
    """Raised when the version file cannot be accessed due to permissions or network issues."""


class VersionFileValidationError(VersionFileError):
    """Raised when the version file content is invalid or cannot be parsed."""


class DecryptionError(Exception):
    """Raised when decryption of a submission file fails."""


class EncryptionError(Exception):
    """Raised when encryption of a submission file fails."""


class NetworkError(Exception):
    """Raised when a network-related operation fails."""


class UploadError(Exception):
    """Raised when uploading a submission file fails."""


class IncompleteSubmissionError(Exception):
    """Raised when a submission is missing required files or metadata."""
