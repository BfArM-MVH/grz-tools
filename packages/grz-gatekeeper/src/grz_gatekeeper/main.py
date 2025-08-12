import logging
import math
import threading
from datetime import UTC, datetime, timedelta
from functools import lru_cache
from typing import TYPE_CHECKING, Annotated

import click
import jwt
import uvicorn
from fastapi import Depends, FastAPI, HTTPException, status
from fastapi.security import OAuth2PasswordBearer, OAuth2PasswordRequestForm
from grz_common.models.base import IgnoringBaseSettings
from grz_common.models.s3 import S3ConfigModel, S3Options
from grz_common.transfer import init_s3_client
from grz_db.models.author import Author
from grz_db.models.submission import SubmissionDb, SubmissionStateEnum
from grz_pydantic_models.gatekeeper import (
    CompleteUploadRequest,
    FileUploadInfo,
    GatekeeperInitiateUploadRequest,
    InitiateUploadResponse,
    PresignedUrlPart,
)
from grz_pydantic_models.submission.metadata.v1 import File, GrzSubmissionMetadata
from grzctl.models.config import DbConfig
from jwt import PyJWTError
from passlib.context import CryptContext
from pydantic import BaseModel, Field

if TYPE_CHECKING:
    from types_boto3_s3 import S3Client
else:
    S3Client = object

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

CONFIG_FILE_PATH: str | None = None

SubmissionId = str
# In-memory store for active upload sessions, protected by a lock
UPLOAD_SESSIONS: dict[SubmissionId, list[FileUploadInfo]] = {}
SESSIONS_LOCK = threading.Lock()


app = FastAPI(
    title="GRZ Gatekeeper API",
    description="API for validating submission metadata and managing multipart uploads.",
    version="0.1.0",
)


class AuthUserConfig(IgnoringBaseSettings):
    hashed_password: str
    disabled: bool = False


class AuthConfig(IgnoringBaseSettings):
    secret_key: str
    algorithm: str
    access_token_expire_minutes: int
    users: dict[str, AuthUserConfig] = Field(default_factory=dict)


class GatekeeperConfig(DbConfig, S3ConfigModel):
    auth: AuthConfig


@lru_cache
def get_gatekeeper_config() -> GatekeeperConfig:
    if CONFIG_FILE_PATH is None:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Server configuration path not set.",
        )
    try:
        return GatekeeperConfig.from_path(CONFIG_FILE_PATH)
    except FileNotFoundError as e:
        log.error(f"Configuration file not found at: {CONFIG_FILE_PATH}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail="Server configuration is missing."
        ) from e


pwd_context = CryptContext(schemes=["bcrypt"], deprecated="auto")
oauth2_scheme = OAuth2PasswordBearer(tokenUrl="/v1/token")


class TokenData(BaseModel):
    username: str | None = None


class AuthUser(BaseModel):
    username: str
    disabled: bool | None = None


class UserInDB(AuthUser):
    hashed_password: str


def get_user_from_config(config: GatekeeperConfig, username: str) -> UserInDB | None:
    """Fetches user details from the loaded configuration."""
    if username in config.auth.users:
        user_config = config.auth.users[username]
        return UserInDB(username=username, hashed_password=user_config.hashed_password, disabled=user_config.disabled)
    return None


def verify_password(plain_password: str, hashed_password: str) -> bool:
    return pwd_context.verify(plain_password, hashed_password)


def authenticate_user(config: GatekeeperConfig, username: str, password: str) -> AuthUser | None:
    user = get_user_from_config(config, username)
    if not user:
        return None
    if not verify_password(password, user.hashed_password):
        return None
    return AuthUser(username=user.username, disabled=user.disabled)


def create_access_token(data: dict, config: GatekeeperConfig) -> str:
    to_encode = data.copy()
    expire = datetime.now(UTC) + timedelta(minutes=config.auth.access_token_expire_minutes)
    to_encode.update({"exp": expire})
    encoded_jwt = jwt.encode(to_encode, config.auth.secret_key, algorithm=config.auth.algorithm)
    return encoded_jwt


async def get_current_user(
    token: Annotated[str, Depends(oauth2_scheme)], config: GatekeeperConfig = Depends(get_gatekeeper_config)
) -> AuthUser:
    credentials_exception = HTTPException(
        status_code=status.HTTP_401_UNAUTHORIZED,
        detail="Could not validate credentials",
        headers={"WWW-Authenticate": "Bearer"},
    )
    try:
        payload = jwt.decode(token, config.auth.secret_key, algorithms=[config.auth.algorithm])
        username: str | None = payload.get("sub")
        if username is None:
            raise credentials_exception
        token_data = TokenData(username=username)
    except PyJWTError as e:
        raise credentials_exception from e

    user_config = get_user_from_config(config, token_data.username)
    if user_config is None:
        raise credentials_exception
    return AuthUser(username=user_config.username, disabled=user_config.disabled)


async def get_current_active_user(current_user: Annotated[AuthUser, Depends(get_current_user)]) -> AuthUser:
    if current_user.disabled:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Inactive user")
    return current_user


class Token(BaseModel):
    access_token: str
    token_type: str


@app.post("/v1/token", response_model=Token)
async def login_for_access_token(
    form_data: Annotated[OAuth2PasswordRequestForm, Depends()],
    config: GatekeeperConfig = Depends(get_gatekeeper_config),
) -> Token:
    user = authenticate_user(config, form_data.username, form_data.password)
    if not user:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Incorrect username or password",
            headers={"WWW-Authenticate": "Bearer"},
        )
    access_token = create_access_token(data={"sub": user.username}, config=config)
    return Token(access_token=access_token, token_type="bearer")  # noqa: S106


def get_submission_db(config: GatekeeperConfig = Depends(get_gatekeeper_config)) -> SubmissionDb:
    if path := config.db.author.private_key_path:
        with open(path, "rb") as f:
            private_key_bytes = f.read()
    elif key := config.db.author.private_key:
        private_key_bytes = key.encode("utf-8")
    else:
        raise ValueError("Either private_key or private_key_path must be provided.")

    author = Author(
        name=config.db.author.name,
        private_key_bytes=private_key_bytes,
        private_key_passphrase=config.db.author.private_key_passphrase,
    )
    return SubmissionDb(db_url=config.db.database_url, author=author)


def get_s3_options(config: GatekeeperConfig = Depends(get_gatekeeper_config)) -> S3Options:
    return config.s3


async def gatekeep(metadata: GrzSubmissionMetadata, db: SubmissionDb):
    submission_id = metadata.submission_id
    if submission := db.get_submission(submission_id):
        latest_state = submission.get_latest_state()
        if not latest_state:
            # Submission exists but has no state yet. This is fine.
            return
        match latest_state.state:
            case SubmissionStateEnum.UPLOADING:
                # This is fine, we can either continue via `UPLOAD_SESSIONS` or abort the existing upload session.
                return
            case _:
                raise HTTPException(
                    status_code=status.HTTP_409_CONFLICT,
                    detail=f"Submission with ID '{submission_id}' is already in state: {latest_state.state}.",
                )


def abort_s3_upload_session(
    s3_client: S3Client,
    bucket_name: str,
    upload_info_list: list[FileUploadInfo],
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
            log.warning(
                f"Upload {file_info.upload_id} for {file_info.encrypted_object_key} not found. It may have been completed or aborted already."
            )
        except Exception:
            log.error(
                f"Failed to abort multipart upload {file_info.upload_id} for {file_info.encrypted_object_key}",
                exc_info=True,
            )


def regenerate_presigned_urls(
    s3_client: S3Client, bucket_name: str, file_info: FileUploadInfo, valid_for: int = 60 * 60 * 2
) -> FileUploadInfo:
    log.info(f"Regenerating URLs for {file_info.encrypted_object_key} (UploadId: {file_info.upload_id})")
    new_parts = []
    for part in file_info.parts:
        new_url = s3_client.generate_presigned_url(
            ClientMethod="upload_part",
            Params={
                "Bucket": bucket_name,
                "Key": file_info.encrypted_object_key,
                "UploadId": file_info.upload_id,
                "PartNumber": part.part_number,
            },
            ExpiresIn=valid_for,
        )
        new_parts.append(
            PresignedUrlPart(
                part_number=part.part_number,
                presigned_url=new_url,
                size=part.size,
                offset=part.offset,
            )
        )
    file_info.parts = new_parts
    return file_info


def generate_presigned_urls_for_file(  # noqa: PLR0913
    s3_client: S3Client,
    bucket_name: str,
    submission_id: str,
    part_size: int,
    file_metadata: File,
    encrypted_size: int,
) -> FileUploadInfo:
    encrypted_file_path = f"{file_metadata.file_path}.c4gh"
    object_key = f"{submission_id}/files/{encrypted_file_path}"

    log.info(f"Generating multipart upload for object key: {object_key}")

    try:
        multipart_upload = s3_client.create_multipart_upload(Bucket=bucket_name, Key=object_key)
        upload_id = multipart_upload["UploadId"]

        file_size = encrypted_size
        num_parts = math.ceil(file_size / part_size) if file_size > 0 else 0

        urls = []
        for i in range(num_parts):
            part_number = i + 1
            offset = i * part_size
            current_part_size = min(part_size, file_size - (i * part_size))

            url = s3_client.generate_presigned_url(
                ClientMethod="upload_part",
                Params={
                    "Bucket": bucket_name,
                    "Key": object_key,
                    "UploadId": upload_id,
                    "PartNumber": part_number,
                },
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


@app.post("/v1/submissions/upload", response_model=InitiateUploadResponse, status_code=status.HTTP_201_CREATED)
async def initiate_upload(
    request: GatekeeperInitiateUploadRequest,
    current_user: Annotated[AuthUser, Depends(get_current_active_user)],
    db: SubmissionDb = Depends(get_submission_db),
    s3_options: S3Options = Depends(get_s3_options),
):
    log.info(
        f"Authenticated user '{current_user.username}' initiated upload for submission '{request.metadata.submission_id}'"
    )
    metadata = request.metadata
    # run basic checks
    await gatekeep(metadata, db)

    # if checks pass, generate pre-signed URLs for all files
    s3_client: S3Client = init_s3_client(s3_options)
    submission_id = metadata.submission_id
    bucket_name = s3_options.bucket
    part_size = s3_options.multipart_chunksize

    with SESSIONS_LOCK:
        if submission_id in UPLOAD_SESSIONS:
            log.info(f"Found existing upload session for submission_id: {submission_id}. Checking for consistency.")
            old_upload_info_list = UPLOAD_SESSIONS[submission_id]

            old_files = {f.file_path: f for f in old_upload_info_list}
            new_files_req = {f.file_path: f for f in request.encrypted_files}

            is_consistent = old_files.keys() == new_files_req.keys()
            if is_consistent:
                for file_path, upload_info in old_files.items():
                    old_size = sum(p.size for p in upload_info.parts)
                    if old_size != new_files_req[file_path].encrypted_size_in_bytes:
                        is_consistent = False
                        break

            if is_consistent:
                log.info(f"Session for {submission_id} is consistent. Regenerating URLs and resuming.")
                refreshed_upload_info = [
                    regenerate_presigned_urls(s3_client, bucket_name, file_info) for file_info in old_upload_info_list
                ]
                return InitiateUploadResponse(submission_id=submission_id, upload_info=refreshed_upload_info)
            else:
                log.warning(f"Session for {submission_id} is inconsistent. Aborting old session and starting new.")
                abort_s3_upload_session(s3_client, bucket_name, old_upload_info_list)
                # Fall through to create a new session

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

        UPLOAD_SESSIONS[submission_id] = new_upload_info_list

        submission = db.get_submission(submission_id)
        if not submission:
            _ = db.add_submission(submission_id)
        db.update_submission_state(submission_id, SubmissionStateEnum.UPLOADING)

        return InitiateUploadResponse(submission_id=submission_id, upload_info=new_upload_info_list)


@app.post("/v1/submissions/{submission_id}/complete", status_code=status.HTTP_200_OK)
async def complete_upload(
    submission_id: str,
    completed_request: CompleteUploadRequest,
    current_user: Annotated[AuthUser, Depends(get_current_active_user)],
    db: SubmissionDb = Depends(get_submission_db),
    s3_options: S3Options = Depends(get_s3_options),
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

    with SESSIONS_LOCK:
        if submission_id in UPLOAD_SESSIONS:
            del UPLOAD_SESSIONS[submission_id]
            log.info(f"Removed completed upload session for submission_id: {submission_id}")

    return {"status": "ok", "message": "Submission successfully uploaded and finalized."}


@click.command()
@click.option("--config-file", help="Path to the configuration file.", required=True)
@click.option("--host", default="127.0.0.1")
@click.option("--port", default=54321)
def run(config_file, host, port):
    global CONFIG_FILE_PATH
    CONFIG_FILE_PATH = config_file
    # validate config file on startup
    _ = get_gatekeeper_config()
    log.info(f"Starting server with config from: {CONFIG_FILE_PATH}")
    uvicorn.run(app, host=host, port=port)


if __name__ == "__main__":
    run()
