"""
This module contains the type definitions for the progress logging states.
"""

from typing import TypedDict


class State(TypedDict, total=False):
    """
    Base state class for progress logging.
    """

    errors: list[str]


class ValidationState(State):
    """
    State for validation progress.
    """

    validation_passed: bool


class UploadState(State):
    """
    State for upload progress.
    """

    upload_successful: bool


class MultipartUploadState(State):
    """
    State for a potentially resumable multipart upload.
    Stores the parameters of the upload session and the progress.
    """

    upload_id: str
    chunk_size: int
    total_parts: int
    completed_parts: dict[int, str]  # Maps part number to ETag


class EncryptionState(State):
    """
    State for encryption progress.
    """

    encryption_successful: bool


class DecryptionState(State):
    """
    State for decryption progress.
    """

    decryption_successful: bool


class DownloadState(State):
    """
    State for download progress.
    """

    download_successful: bool
