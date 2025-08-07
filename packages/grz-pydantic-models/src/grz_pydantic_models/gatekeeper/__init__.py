"""Pydantic models for the grz-gatekeeper API."""

from pydantic import Field

from ..common import StrictBaseModel
from ..submission.metadata import GrzSubmissionMetadata


class EncryptedFileInfo(StrictBaseModel):
    """Size information for a single encrypted file."""

    file_path: str = Field(..., description="The relative path of the file as in the metadata.")
    encrypted_size_in_bytes: int = Field(..., description="The size of the encrypted .c4gh file in bytes.")


class GatekeeperInitiateUploadRequest(StrictBaseModel):
    """The request model for the gatekeeper's /upload endpoint."""

    metadata: GrzSubmissionMetadata
    encrypted_files: list[EncryptedFileInfo]


class PresignedUrlPart(StrictBaseModel):
    """A pre-signed URL for a part of a file in a multipart upload."""

    part_number: int = Field(..., description="The part number for the multipart upload.")
    presigned_url: str = Field(..., description="The pre-signed URL for uploading this part.")
    size: int = Field(..., description="The size of this part in bytes.")
    offset: int = Field(..., description="The byte offset from the beginning of the file for this part.")


class FileUploadInfo(StrictBaseModel):
    """All information required to upload a single file."""

    file_path: str = Field(..., description="The relative path of the file as in the metadata.")
    encrypted_object_key: str = Field(..., description="The full S3 object key for the encrypted file.")
    upload_id: str = Field(..., description="The ID for the multipart upload on S3.")
    parts: list[PresignedUrlPart] = Field(..., description="A list of pre-signed URLs for each part of the file.")


class InitiateUploadResponse(StrictBaseModel):
    """The response from the gatekeeper after initiating an upload."""

    submission_id: str = Field(..., description="The unique identifier for this submission.")
    upload_info: list[FileUploadInfo] = Field(..., description="Upload information for all files in the submission.")


class CompletedPart(StrictBaseModel):
    """Information about a successfully uploaded part."""

    part_number: int
    etag: str = Field(..., description="The ETag returned by S3 for the uploaded part.")


class CompletedFileUpload(StrictBaseModel):
    """Information about a fully uploaded file to be sent for completion."""

    file_path: str = Field(..., description="The relative path of the file, for logging and verification.")
    upload_id: str = Field(..., description="The ID for the multipart upload on S3.")
    parts: list[CompletedPart]


class CompleteUploadRequest(StrictBaseModel):
    """The request to the gatekeeper to finalize a submission upload."""

    files: list[CompletedFileUpload]
