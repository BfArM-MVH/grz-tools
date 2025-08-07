import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from os import PathLike
from os.path import getsize
from pathlib import Path
from typing import override

import requests
from grz_pydantic_models.gatekeeper import (
    CompletedFileUpload,
    CompletedPart,
    CompleteUploadRequest,
    EncryptedFileInfo,
    FileUploadInfo,
    GatekeeperInitiateUploadRequest,
    InitiateUploadResponse,
)
from grz_pydantic_models.submission.metadata import File as SubmissionFileMetadata
from grz_pydantic_models.submission.metadata import GrzSubmissionMetadata
from tqdm.auto import tqdm

from ...constants import TQDM_DEFAULTS
from ...progress import FileProgressLogger, MultipartUploadState
from ...workers.upload import UploadError, UploadWorker, log
from ..submission import EncryptedSubmission


class GrzGatekeeperUploadWorker(UploadWorker):
    """
    Implementation of an upload worker that uses the grz-gatekeeper API to get
    pre-signed URLs and manages a parallel multipart upload.
    """

    __log = log.getChild("GrzGatekeeperUploadWorker")

    def __init__(
        self,
        api_base_url: str,
        status_file_path: str | PathLike,
        threads: int = 1,
    ):
        super().__init__()
        self._api_base_url = api_base_url
        self._threads = threads
        self._progress_logger = FileProgressLogger[MultipartUploadState](status_file_path)

    @override
    def upload_file(self, local_file_path: str | PathLike, s3_object_id: str):
        raise NotImplementedError("GrzGatekeeperUploadWorker does not support direct single-file uploads.")

    @override
    def archive(self, encrypted_submission: EncryptedSubmission):
        raise NotImplementedError("archive is not yet implemented.")

    @classmethod
    def _upload_chunk_worker(cls, chunk_data: bytes, part_info: dict) -> dict:
        try:
            response = requests.put(part_info["presigned_url"], data=chunk_data, timeout=300)
            response.raise_for_status()
            etag = response.headers.get("ETag")
            if not etag:
                raise ValueError("ETag not found in S3 response header.")

            return {"part_number": part_info["part_number"], "etag": etag.strip('"'), "size": len(chunk_data)}
        except Exception as e:
            log.error(f"Upload failed for part {part_info.get('part_number', 'N/A')}: {e}")
            raise

    def _initiate_upload(self, encrypted_submission: EncryptedSubmission) -> InitiateUploadResponse:
        self.__log.info("Initiating upload via grz-gatekeeper API…")
        metadata: GrzSubmissionMetadata = encrypted_submission.metadata.content

        encrypted_files_info = [
            EncryptedFileInfo(file_path=original_meta.file_path, encrypted_size_in_bytes=getsize(encrypted_path))
            for encrypted_path, original_meta in encrypted_submission.encrypted_files.items()
        ]

        initiate_request = GatekeeperInitiateUploadRequest(metadata=metadata, encrypted_files=encrypted_files_info)
        initiate_url = f"{self._api_base_url}/v1/submissions/upload"

        response = requests.post(
            initiate_url,
            data=initiate_request.model_dump_json(by_alias=True),
            headers={"Content-Type": "application/json"},
            timeout=60,
        )

        if not response.ok:
            self.__log.error(f"API request to {initiate_url} failed with status {response.status_code}.")
            try:
                error_detail = response.json()
                self.__log.error(f"API Error Detail: {error_detail.get('detail', response.text)}")
            except json.JSONDecodeError:
                self.__log.error(f"API Response: {response.text}")
            response.raise_for_status()

        upload_data = InitiateUploadResponse.model_validate(response.json())
        self.__log.info(f"Upload initiated for Submission ID: {upload_data.submission_id}")
        return upload_data

    def _get_or_create_session_state(
        self, local_file_path: Path, original_file_metadata: SubmissionFileMetadata, file_to_upload: FileUploadInfo
    ) -> MultipartUploadState:
        current_total_parts = len(file_to_upload.parts)
        current_chunk_size = file_to_upload.parts[0].size if current_total_parts > 0 else -1
        current_upload_id = file_to_upload.upload_id

        saved_state = self._progress_logger.get_state(local_file_path, original_file_metadata)
        if saved_state:
            # Check if upload parameters have changed, invalidating the saved state
            params_match = (
                saved_state.get("upload_id") == current_upload_id
                and saved_state["chunk_size"] == current_chunk_size
                and saved_state["total_parts"] == current_total_parts
            )
            if not params_match:
                self.__log.warning(
                    f"Upload parameters for {local_file_path.name} have changed. Restarting upload for this file."
                )
                saved_state = None

        if saved_state:
            saved_state["completed_parts"] = {int(k): v for k, v in saved_state["completed_parts"].items()}
            return saved_state
        else:
            return MultipartUploadState(
                upload_id=current_upload_id,
                chunk_size=current_chunk_size,
                total_parts=current_total_parts,
                completed_parts={},
            )

    def _upload_file_parts(  # noqa: PLR0913
        self,
        local_file_path: Path,
        parts_to_upload: list[dict],
        session_state: MultipartUploadState,
        original_file_metadata: SubmissionFileMetadata,
        pbar_total: tqdm,
        pbar_file: tqdm,
    ):
        parts_to_upload.sort(key=lambda p: p["offset"])

        with ThreadPoolExecutor(max_workers=self._threads) as executor, open(local_file_path, "rb") as f:
            futures = {
                executor.submit(GrzGatekeeperUploadWorker._upload_chunk_worker, f.read(part["size"]), part): part
                for part in parts_to_upload
                if f.seek(part["offset"]) == part["offset"]
            }

            try:
                for future in as_completed(futures):
                    result = future.result()
                    part_num = result["part_number"]
                    session_state["completed_parts"][part_num] = result["etag"]
                    self._progress_logger.set_state(local_file_path, original_file_metadata, session_state)
                    pbar_file.update(result["size"])
                    pbar_total.update(result["size"])
            except Exception as e:
                self.__log.error(f"An upload worker failed for {local_file_path.name}. Aborting submission.")
                raise UploadError(f"Upload failed for file {local_file_path.name}") from e

    def _complete_upload(self, submission_id: str, completed_files: list[CompletedFileUpload]) -> None:
        self.__log.info("All parts uploaded. Finalizing submission with the API…")
        complete_payload = CompleteUploadRequest(files=completed_files)
        complete_url = f"{self._api_base_url}/v1/submissions/{submission_id}/complete"
        response = requests.post(
            complete_url,
            data=complete_payload.model_dump_json(by_alias=True),
            timeout=60,
            headers={"Content-Type": "application/json"},
        )
        response.raise_for_status()

    @override
    def upload(self, encrypted_submission: EncryptedSubmission):
        upload_data = self._initiate_upload(encrypted_submission)
        submission_id = upload_data.submission_id

        relative_path_to_metadata: dict[str, SubmissionFileMetadata] = {
            meta.file_path: meta for meta in encrypted_submission.encrypted_files.values()
        }
        all_completed_files: list[CompletedFileUpload] = []
        total_upload_size = sum(
            getsize(encrypted_submission.encrypted_files_dir / f"{f.file_path}.c4gh") for f in upload_data.upload_info
        )

        pbar_total = None
        pbar_file = None
        try:
            pbar_total = tqdm(total=total_upload_size, desc="OVERALL ", position=0, **TQDM_DEFAULTS)  # type: ignore[call-overload]
            pbar_file = tqdm(total=0, desc="UPLOAD  ", position=1, **TQDM_DEFAULTS)  # type: ignore[call-overload]

            for file_to_upload in upload_data.upload_info:
                local_file_path = encrypted_submission.encrypted_files_dir / f"{file_to_upload.file_path}.c4gh"
                original_meta = relative_path_to_metadata[file_to_upload.file_path]
                session_state = self._get_or_create_session_state(local_file_path, original_meta, file_to_upload)

                all_parts_info = [p.model_dump() for p in file_to_upload.parts]
                completed_part_numbers = set(session_state["completed_parts"].keys())

                already_uploaded_bytes = sum(
                    p["size"] for p in all_parts_info if p["part_number"] in completed_part_numbers
                )
                pbar_total.update(already_uploaded_bytes)

                parts_to_upload = [p for p in all_parts_info if p["part_number"] not in completed_part_numbers]

                if parts_to_upload:
                    pbar_file.reset(total=getsize(local_file_path))
                    pbar_file.set_postfix_str(f"{local_file_path.name}", refresh=True)
                    pbar_file.update(already_uploaded_bytes)
                    try:
                        self._upload_file_parts(
                            local_file_path, parts_to_upload, session_state, original_meta, pbar_total, pbar_file
                        )
                        pbar_file.set_postfix_str(f"✓ OK    {local_file_path.name}", refresh=True)
                    except UploadError as e:
                        pbar_file.set_postfix_str(f"✗ ERROR {local_file_path.name}", refresh=True)
                        raise e

                final_parts = [
                    CompletedPart(part_number=pn, etag=etag) for pn, etag in session_state["completed_parts"].items()
                ]
                all_completed_files.append(
                    CompletedFileUpload(
                        file_path=file_to_upload.file_path,
                        upload_id=file_to_upload.upload_id,
                        parts=sorted(final_parts, key=lambda p: p.part_number),
                    )
                )

            pbar_total.set_postfix_str("~ Finalizing upload…", refresh=True)
            try:
                self._complete_upload(submission_id, all_completed_files)
                pbar_total.set_postfix_str("✓ Upload complete!", refresh=True)
            except Exception as e:
                pbar_total.set_postfix_str("✗ Upload failed!", refresh=True)
                raise e

        finally:
            if pbar_file:
                pbar_file.leave = False
                pbar_file.close()
            if pbar_total:
                pbar_total.close()

        self.__log.info("Upload finished successfully for submission %s!", submission_id)
