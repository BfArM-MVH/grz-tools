class VersionFileError(Exception):
    """Base exception for version file related errors."""


class VersionFileNotFoundError(VersionFileError):
    """Raised when the version file is missing from S3."""


class VersionFileAccessError(VersionFileError):
    """Raised when the version file cannot be accessed due to permissions or network issues."""


class VersionFileValidationError(VersionFileError):
    """Raised when the version file content is invalid or cannot be parsed."""


class GrzError(Exception):
    """Base of every failure that grz-tools accounts for, as opposed to a bug.

    grzctl records a failure reason for each subclass, so a new subclass needs one.
    """


class SubmissionRejectedError(GrzError):
    """The submission itself is at fault, so the submitter has to send a corrected one."""


class MissingSubmissionFileError(SubmissionRejectedError):
    """The inbox lacks the metadata or a file that the metadata lists."""


class SubmissionValidationError(SubmissionRejectedError):
    """The metadata or the content of a file breaks the specification."""


class DecryptionError(SubmissionRejectedError):
    """A submission file cannot be decrypted."""


class DuplicateUploadError(SubmissionRejectedError):
    """The bucket already holds a submission with this ID, so it was uploaded before."""


class IncompleteSubmissionError(GrzError):
    """A step ran before an earlier step had passed for every file of the submission."""


class SubmissionCleanedError(GrzError):
    """``grzctl clean`` has started on the submission, so the inbox no longer holds it."""


class ConfigurationError(GrzError):
    """The setup is wrong or incomplete, such as a missing key or credentials that a service rejects."""


class TransferError(GrzError):
    """Moving data to or from S3 or BfArM failed."""


class DownloadError(TransferError):
    """Reading from S3 failed."""


class MissingObjectError(DownloadError):
    """S3 holds no object under the key.

    Only the caller knows whose object it is, so it turns a missing submission file into a
    :class:`MissingSubmissionFileError`.
    """


class UploadError(TransferError):
    """Writing to S3 failed."""


class NetworkError(TransferError):
    """A request to S3 or BfArM did not get through, such as for a lost connection, a timeout, or a server error."""


class EncryptionError(GrzError):
    """Encrypting a submission file failed."""


class DetailedQCError(GrzError):
    """The detailed QC workflow failed."""


class PruefberichtGenerationError(GrzError):
    """The Prüfbericht cannot be generated, for example because the database lacks a field that it needs."""


class PruefberichtRejectedError(GrzError):
    """BfArM rejected the Prüfbericht."""
