import logging
from typing import Annotated

from fastapi import Depends, FastAPI, HTTPException, status
from fastapi.security import OAuth2PasswordRequestForm
from grz_common.models.s3 import S3Options
from grz_common.transfer import S3Client, init_s3_client
from grz_db.models.submission import SubmissionDb, SubmissionStateEnum
from grz_pydantic_models.gatekeeper import (
    CompleteUploadRequest,
    FileUploadInfo,
    GatekeeperInitiateUploadRequest,
    InitiateUploadResponse,
    PresignedUrlPart,
)
from sqlmodel import Session, select

from .auth import (
    AuthUser,
    Token,
    authenticate_user,
    create_access_token,
    get_current_active_user,
)
from .config import GatekeeperConfig, get_gatekeeper_config
from .dependencies import get_gatekeeper_db_session, get_s3_options, get_submission_db
from .rate_limiter import RateLimiter
from .services import (
    abort_s3_upload_session,
    gatekeep,
    generate_presigned_urls_for_file,
    regenerate_presigned_urls,
)
from .session_store import UploadSession

log = logging.getLogger(__name__)

app = FastAPI(
    title="GRZ Gatekeeper API",
    description="API for validating submission metadata and managing multipart uploads.",
    version="0.1.0",
)


@app.post("/v1/token", response_model=Token, dependencies=[Depends(RateLimiter(times=10, seconds=60))])
async def login_for_access_token(
    form_data: Annotated[OAuth2PasswordRequestForm, Depends()],
    gatekeeper_config: GatekeeperConfig = Depends(get_gatekeeper_config),
) -> Token:
    user = authenticate_user(gatekeeper_config, form_data.username, form_data.password)
    if not user:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Incorrect username or password",
            headers={"WWW-Authenticate": "Bearer"},
        )
    access_token = create_access_token(data={"sub": user.username}, config=gatekeeper_config)
    return Token(access_token=access_token, token_type="bearer")  # noqa: S106


@app.post(
    "/v1/submissions/upload",
    response_model=InitiateUploadResponse,
    status_code=status.HTTP_201_CREATED,
    dependencies=[Depends(RateLimiter(times=10, seconds=60))],
)
async def initiate_upload(  # noqa: C901, PLR0912, PLR0913, PLR0915
    request: GatekeeperInitiateUploadRequest,
    current_user: Annotated[AuthUser, Depends(get_current_active_user)],
    submission_db: SubmissionDb = Depends(get_submission_db),
    s3_options: S3Options = Depends(get_s3_options),
    gatekeeper_config: GatekeeperConfig = Depends(get_gatekeeper_config),
    session: Session = Depends(get_gatekeeper_db_session),
):
    log.info(
        f"Authenticated user '{current_user.username}' initiated upload for submission '{request.metadata.submission_id}'"
    )
    metadata = request.metadata
    # run basic checks
    await gatekeep(metadata, submission_db)

    # if checks pass, generate pre-signed URLs for all files
    s3_client: S3Client = init_s3_client(s3_options)
    submission_id = metadata.submission_id
    bucket_name = s3_options.bucket
    part_size = s3_options.multipart_chunksize

    statement = select(UploadSession).where(UploadSession.submission_id == submission_id)
    old_upload_sessions = session.exec(statement).all()

    if old_upload_sessions:
        log.info(f"Found existing upload session for submission_id: {submission_id}. Checking for consistency.")
        old_files = {f.file_path: f for f in old_upload_sessions}
        new_files_req = {f.file_path: f for f in request.encrypted_files}
        is_consistent = old_files.keys() == new_files_req.keys()
        if is_consistent:
            for file_path, upload_info in old_files.items():
                parts = map(lambda p: PresignedUrlPart(**p), upload_info.parts)
                old_size = sum(p.size for p in parts)
                if old_size != new_files_req[file_path].encrypted_size_in_bytes:
                    is_consistent = False
                    break
        if is_consistent:
            log.info(f"Session for {submission_id} is consistent. Regenerating URLs and resuming.")
            refreshed_info = [regenerate_presigned_urls(s3_client, bucket_name, f) for f in old_upload_sessions]
            return InitiateUploadResponse(submission_id=submission_id, upload_info=refreshed_info)
        else:
            log.warning(f"Session for {submission_id} is inconsistent. Aborting old session and starting new.")
            abort_s3_upload_session(s3_client, bucket_name, list(old_upload_sessions))
            for old_session in old_upload_sessions:
                session.delete(old_session)

    log.info(f"Creating new upload session for submission_id: {submission_id}")
    path_to_metadata = {
        file.file_path: file
        for donor in metadata.donors
        for lab_datum in donor.lab_data
        if lab_datum.sequence_data
        for file in lab_datum.sequence_data.files
    }
    path_to_encrypted_size = {info.file_path: info.encrypted_size_in_bytes for info in request.encrypted_files}

    new_upload_info_list: list[FileUploadInfo] = []
    for file_path, encrypted_size in path_to_encrypted_size.items():
        if file_path in path_to_metadata:
            file_meta = path_to_metadata[file_path]
            info = generate_presigned_urls_for_file(
                s3_client, bucket_name, submission_id, part_size, file_meta, encrypted_size
            )
            new_upload_info_list.append(info)
        else:
            log.warning(
                f"File path '{file_path}' from encrypted_files list not found in submission metadata. Skipping."
            )

    # session cleanup and creation
    count_statement = select(UploadSession)
    num_sessions = len(session.exec(count_statement).all())
    max_sessions = gatekeeper_config.max_active_sessions
    if num_sessions >= max_sessions:
        limit = max((num_sessions - max_sessions) + 1, 1)
        exceeding_sessions = session.exec(select(UploadSession).order_by(UploadSession.created_at).limit(limit)).all()
        if exceeding_sessions:
            log.warning(f"Max sessions ({max_sessions}) reached. Deleting {limit} oldest sessions.")
            for exceeding_session in exceeding_sessions:
                session.delete(exceeding_session)

    for file_info in new_upload_info_list:
        db_session_object = UploadSession(**file_info.model_dump(), submission_id=submission_id)
        session.add(db_session_object)

    submission = submission_db.get_submission(submission_id)
    if not submission:
        _ = submission_db.add_submission(submission_id)
    submission_db.update_submission_state(submission_id, SubmissionStateEnum.UPLOADING)

    session.commit()
    return InitiateUploadResponse(submission_id=submission_id, upload_info=new_upload_info_list)


@app.post(
    "/v1/submissions/{submission_id}/complete",
    status_code=status.HTTP_200_OK,
    dependencies=[Depends(RateLimiter(times=10, seconds=60))],
)
async def complete_upload(  # noqa: PLR0913
    submission_id: str,
    completed_request: CompleteUploadRequest,
    current_user: Annotated[AuthUser, Depends(get_current_active_user)],
    db: SubmissionDb = Depends(get_submission_db),
    s3_options: S3Options = Depends(get_s3_options),
    session: Session = Depends(get_gatekeeper_db_session),
):
    log.info(f"Authenticated user '{current_user.username}' is completing upload for submission '{submission_id}'")
    s3_client: S3Client = init_s3_client(s3_options)
    bucket_name = s3_options.bucket

    submission = db.get_submission(submission_id)
    if not submission:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Submission not found.")

    latest_state = submission.get_latest_state()
    if not latest_state or latest_state.state != SubmissionStateEnum.UPLOADING:
        raise HTTPException(
            status_code=status.HTTP_409_CONFLICT, detail=f"Submission is not in state {SubmissionStateEnum.UPLOADING}."
        )

    for file_info in completed_request.files:
        object_key = f"{submission_id}/files/{file_info.file_path}.c4gh"

        parts = {"Parts": [{"ETag": part.etag, "PartNumber": part.part_number} for part in file_info.parts]}

        log.info(f"Completing multipart upload for object: {object_key}")
        try:
            s3_client.complete_multipart_upload(
                Bucket=bucket_name,
                Key=object_key,
                UploadId=file_info.upload_id,
                MultipartUpload=parts,
            )
        except Exception as e:
            db.update_submission_state(
                submission_id, SubmissionStateEnum.ERROR, data={"reason": "Failed to complete S3 upload"}
            )
            log.error(f"Failed to complete S3 multipart upload for {object_key}: {e}", exc_info=True)
            raise HTTPException(
                status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                detail=f"Failed to finalize file: {file_info.file_path}",
            ) from e

    db.update_submission_state(submission_id, SubmissionStateEnum.UPLOADED)

    statement = select(UploadSession).where(UploadSession.submission_id == submission_id)
    results = session.exec(statement).all()
    for db_info in results:
        session.delete(db_info)
    session.commit()
    log.info(f"Removed completed upload session from DB for submission_id: {submission_id}")

    return {"status": "ok", "message": "Submission successfully uploaded and finalized."}
