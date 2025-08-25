import logging
import math
from datetime import UTC, datetime, timedelta

from fastapi import HTTPException, status
from grz_common.transfer import S3Client
from grz_db.models.submission import SubmissionDb, SubmissionStateEnum
from grz_pydantic_models.gatekeeper import FileUploadInfo, PresignedUrlPart
from grz_pydantic_models.submission.metadata.v1 import File, GrzSubmissionMetadata, SubmissionType
from sqlmodel import Session, select

from grz_gatekeeper.session_store import UploadSession

log = logging.getLogger(__name__)


async def gatekeep(metadata: GrzSubmissionMetadata, db: SubmissionDb):
    submission_id = metadata.submission_id
    if submission := db.get_submission(submission_id):
        latest_state = submission.get_latest_state()
        if not latest_state:
            return
        if latest_state.state not in [SubmissionStateEnum.UPLOADING, None]:
            raise HTTPException(
                status_code=status.HTTP_409_CONFLICT,
                detail=f"Submission with ID '{submission_id}' is already in state: {latest_state.state}.",
            )
    elif metadata.submission.submission_type not in [SubmissionType.initial, SubmissionType.test]:
        raise HTTPException(
            status_code=status.HTTP_409_CONFLICT,
            detail=f"Submission with ID '{submission_id}' is not an initial or test submission "
            f"but found no records of matching submissions in database.",
        )


def abort_s3_upload_session(
    s3_client: S3Client,
    bucket_name: str,
    upload_info_list: list[FileUploadInfo] | list[UploadSession],
):
    log.warning("Aborting stale or inconsistent upload session.")
    for file_info in upload_info_list:
        try:
            log.info(f"Aborting upload for key: {file_info.encrypted_object_key} (UploadId: {file_info.upload_id})")
            s3_client.abort_multipart_upload(
                Bucket=bucket_name,
                Key=file_info.encrypted_object_key,
                UploadId=file_info.upload_id,
            )
        except s3_client.exceptions.NoSuchUpload:
            log.warning(f"Upload {file_info.upload_id} for {file_info.encrypted_object_key} not found.")
        except Exception:
            log.error(f"Failed to abort multipart upload {file_info.upload_id}", exc_info=True)


def regenerate_presigned_urls(
    s3_client: S3Client, bucket_name: str, file_info: FileUploadInfo | UploadSession, valid_for: int = 60 * 60 * 2
) -> FileUploadInfo:
    log.info(f"Regenerating URLs for {file_info.encrypted_object_key} (UploadId: {file_info.upload_id})")
    new_parts = []
    parts = map(lambda p: PresignedUrlPart(**p), file_info.parts)
    for part in parts:
        new_url = s3_client.generate_presigned_url(
            "upload_part",
            Params={
                "Bucket": bucket_name,
                "Key": file_info.encrypted_object_key,
                "UploadId": file_info.upload_id,
                "PartNumber": part.part_number,
            },
            ExpiresIn=valid_for,
        )
        new_parts.append(
            PresignedUrlPart(part_number=part.part_number, presigned_url=new_url, size=part.size, offset=part.offset)
        )
    file_info.parts = new_parts
    file_info.last_modified_at = datetime.now(UTC)
    return file_info


def generate_presigned_urls_for_file(  # noqa: PLR0913
    s3_client: S3Client, bucket_name: str, submission_id: str, part_size: int, file_metadata: File, encrypted_size: int
) -> FileUploadInfo:
    encrypted_file_path = f"{file_metadata.file_path}.c4gh"
    object_key = f"{submission_id}/files/{encrypted_file_path}"
    log.info(f"Generating multipart upload for object key: {object_key}")

    try:
        multipart_upload = s3_client.create_multipart_upload(Bucket=bucket_name, Key=object_key)
        upload_id = multipart_upload["UploadId"]
        num_parts = math.ceil(encrypted_size / part_size) if encrypted_size > 0 else 0
        urls = []
        for i in range(num_parts):
            part_number = i + 1
            offset = i * part_size
            current_part_size = min(part_size, encrypted_size - offset)
            url = s3_client.generate_presigned_url(
                "upload_part",
                Params={"Bucket": bucket_name, "Key": object_key, "UploadId": upload_id, "PartNumber": part_number},
                ExpiresIn=60 * 60 * 2,
            )
            urls.append(
                PresignedUrlPart(part_number=part_number, presigned_url=url, size=current_part_size, offset=offset)
            )
        return FileUploadInfo(
            file_path=file_metadata.file_path, encrypted_object_key=object_key, upload_id=upload_id, parts=urls
        )
    except Exception as e:
        log.error(f"Failed to create multipart upload for {object_key}: {e}", exc_info=True)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail="Could not generate S3 pre-signed URLs."
        ) from e


def session_cleanup(session: Session, s3_client: S3Client, bucket_name: str, older_than_days: int = 10):
    """
    Find and remove stale upload sessions from the database and aborts them.
    """
    stale_threshold = datetime.now(UTC) - timedelta(days=older_than_days)
    log.info(f"Searching for sessions not modified since {stale_threshold.isoformat()}")

    statement = select(UploadSession).where(UploadSession.last_modified_at < stale_threshold)
    stale_sessions = session.exec(statement).all()

    if not stale_sessions:
        log.info("No stale sessions found.")
        return

    log.warning(f"Found {len(stale_sessions)} stale session(s) to clean up.")
    abort_s3_upload_session(s3_client, bucket_name, list(stale_sessions))

    for stale_session in stale_sessions:
        session.delete(stale_session)

    session.commit()
    log.info("Finished cleaning stale sessions.")
