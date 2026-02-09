"""S3 key construction utilities for submissions."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from grz_pydantic_models.submission.metadata.v1 import File as SubmissionFileMetadata


class S3KeyBuilder:
    """Utility for constructing S3 keys for submission files."""

    @staticmethod
    def get_encrypted_path(file_path: str | Path) -> str:
        """
        Get the encrypted file path by adding .c4gh extension.

        :param file_path: Original file path
        :returns: Path with .c4gh extension
        """
        path_obj = Path(file_path) if isinstance(file_path, str) else file_path
        return str(path_obj) + ".c4gh"

    @staticmethod
    def source_key(submission_id: str, file_metadata: SubmissionFileMetadata) -> str:
        """
        Construct the source S3 key (inbox) for a file.

        :param submission_id: Submission ID
        :param file_metadata: File metadata
        :returns: S3 key in format: {submission_id}/files/{encrypted_path}
        """
        encrypted_path = S3KeyBuilder.get_encrypted_path(file_metadata.file_path)
        return f"{submission_id}/files/{encrypted_path}"

    @staticmethod
    def target_key(submission_id: str, file_metadata: SubmissionFileMetadata) -> str:
        """
        Construct the target S3 key (archive) for a file.

        :param submission_id: Submission ID
        :param file_metadata: File metadata
        :returns: S3 key in format: {submission_id}/files/{encrypted_path}
        """
        # Target uses the same structure as source
        return S3KeyBuilder.source_key(submission_id, file_metadata)

    @staticmethod
    def metadata_key(submission_id: str) -> str:
        """
        Construct the metadata S3 key.

        :param submission_id: Submission ID
        :returns: S3 key in format: {submission_id}/metadata/metadata.json
        """
        return f"{submission_id}/metadata/metadata.json"

    @staticmethod
    def log_key(submission_id: str, log_filename: str) -> str:
        """
        Construct a log file S3 key.

        :param submission_id: Submission ID
        :param log_filename: Log file name
        :returns: S3 key in format: {submission_id}/logs/{log_filename}
        """
        return f"{submission_id}/logs/{log_filename}"
