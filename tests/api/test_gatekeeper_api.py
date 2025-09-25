import importlib.resources
from unittest.mock import MagicMock, call, patch

from grz_db.models.submission import SubmissionStateEnum
from grz_pydantic_models.gatekeeper import (
    CompletedFileUpload,
    CompletedPart,
    CompleteUploadRequest,
    EncryptedFileInfo,
    GatekeeperInitiateUploadRequest,
    InitiateUploadResponse,
)
from grz_pydantic_models.submission.metadata import GrzSubmissionMetadata

from .. import mock_files


def test_get_token_success(test_app_client):
    """Test successful authentication and token retrieval."""
    response = test_app_client.post(
        "/v1/token",
        data={"username": "le-123456789", "password": "a_very_secret_password"},
    )
    assert response.status_code == 200
    token_data = response.json()
    assert "access_token" in token_data
    assert token_data["token_type"] == "bearer"


def test_get_token_failure(test_app_client):
    """Test authentication failure with incorrect credentials."""
    response = test_app_client.post("/v1/token", data={"username": "le-123456789", "password": "wrong_password"})
    assert response.status_code == 401
    assert response.json()["detail"] == "Incorrect username or password"


@patch("grz_gatekeeper.main.init_s3_client")
def test_upload_success(mock_init_s3, test_app_client, mock_submission_db):
    """Test full, successful upload."""
    # mock client
    mock_s3_client = mock_init_s3.return_value
    mock_s3_client.create_multipart_upload.return_value = {"UploadId": "mock-upload-id"}
    mock_s3_client.generate_presigned_url.return_value = "https://s3.mock/presigned-url"

    # the submission is not registered in the submission DB _yet_
    mock_submission_db.get_submission.return_value = None

    # token authentication
    token_res = test_app_client.post(
        "/v1/token",
        data={"username": "le-123456789", "password": "a_very_secret_password"},
    )
    token = token_res.json()["access_token"]
    headers = {"Authorization": f"Bearer {token}"}

    # prepare and initiate upload
    metadata_str = (
        importlib.resources.files(mock_files)
        .joinpath("submissions", "valid_submission", "metadata", "metadata.json")
        .read_text()
    )
    metadata = GrzSubmissionMetadata.model_validate_json(metadata_str)
    submission_id = metadata.submission_id

    encrypted_files = [
        EncryptedFileInfo(
            file_path="aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_normal.read1.fastq.gz",
            encrypted_size_in_bytes=165912,
        )
    ]
    initiate_payload = GatekeeperInitiateUploadRequest(metadata=metadata, encrypted_files=encrypted_files)

    initiate_response = test_app_client.post(
        "/v1/submissions/upload",
        headers=headers,
        content=initiate_payload.model_dump_json(by_alias=True),
    )
    assert initiate_response.status_code == 201, initiate_response.text

    initiate_data = InitiateUploadResponse.model_validate(initiate_response.json())

    assert initiate_data.submission_id == submission_id
    assert len(initiate_data.upload_info) == 1

    upload_info = initiate_data.upload_info[0]
    assert upload_info.upload_id == "mock-upload-id"
    assert upload_info.parts[0].presigned_url == "https://s3.mock/presigned-url"

    # the submission should now exist in the DB and be in the "Uploading" state
    mock_submission_state = MagicMock()
    mock_submission_state.state = SubmissionStateEnum.UPLOADING
    mock_submission = MagicMock()
    mock_submission.get_latest_state.return_value = mock_submission_state
    mock_submission_db.get_submission.return_value = mock_submission

    # complete upload
    completion_payload = CompleteUploadRequest(
        files=[
            CompletedFileUpload(
                file_path=upload_info.file_path,
                upload_id=upload_info.upload_id,
                parts=[CompletedPart(part_number=1, etag="mock-etag-for-part-1")],
            )
        ]
    )

    complete_response = test_app_client.post(
        f"/v1/submissions/{submission_id}/complete",
        headers=headers,
        content=completion_payload.model_dump_json(by_alias=True),
    )

    # check for successful status
    assert complete_response.status_code == 200
    assert complete_response.json() == {"status": "ok", "message": "Submission successfully uploaded and finalized."}

    # verify that the submission db was updated correctly
    mock_submission_db.add_submission.assert_called_once_with(submission_id)

    # check that the state was updated to 'Uploading' and then to 'Uploaded'
    state_update_calls = [
        call(submission_id, "Uploading"),
        call(submission_id, "Uploaded"),
    ]
    mock_submission_db.update_submission_state.assert_has_calls(state_update_calls)
