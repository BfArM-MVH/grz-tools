"""Module for uploading encrypted submissions to a remote storage"""

from __future__ import annotations

import abc
import json
import logging
import math
import re
import shutil
from os import PathLike
from os.path import getsize
from pathlib import Path
from typing import TYPE_CHECKING, override

import botocore.handlers
from boto3.s3.transfer import S3Transfer, TransferConfig  # type: ignore[import-untyped]
from tqdm.auto import tqdm

from ..constants import TQDM_SMOOTHING
from ..models.s3 import S3Options
from ..progress import FileProgressLogger, UploadState
from ..transfer import init_s3_client

MULTIPART_THRESHOLD = 8 * 1024**2  # 8MiB, boto3 default, largely irrelevant
MULTIPART_MAX_CHUNKS = 1000  # CEPH S3 limit, AWS limit is 10000

if TYPE_CHECKING:
    from .submission import EncryptedSubmission

log = logging.getLogger(__name__)

# see discussion: https://github.com/boto/boto3/discussions/4251 for acception bucketnames with : in the name
botocore.handlers.VALID_BUCKET = re.compile(r"^[:a-zA-Z0-9.\-_]{1,255}$")  # type: ignore[import-untyped]


class UploadError(Exception):
    """Exception raised when an upload fails"""

    pass


class UploadWorker(metaclass=abc.ABCMeta):
    """Worker baseclass for uploading encrypted submissions"""

    @abc.abstractmethod
    def upload(self, encrypted_submission: EncryptedSubmission):
        """
        Upload an encrypted submission to a GRZ inbox

        :param encrypted_submission: The encrypted submission to upload
        :raises UploadError: when the upload failed
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def upload_file(self, local_file_path: str | PathLike, s3_object_id: str):
        """
        Upload a single file to the specified object ID
        :param local_file_path: Path to the file to upload
        :param s3_object_id: Remote S3 object ID under which the file should be stored
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def archive(self, encrypted_submission: EncryptedSubmission):
        """
        Archive an encrypted submission within a GRZ

        :param encrypted_submission: The encrypted submission to archive
        :raises UploadError: when archival failed
        """
        raise NotImplementedError()


class S3BotoUploadWorker(UploadWorker):
    """Implementation of an upload worker using boto3 for S3"""

    __log = log.getChild("S3BotoUploadWorker")

    def __init__(
        self,
        s3_options: S3Options,
        status_file_path: str | PathLike,
        threads: int = 1,
    ):
        """
        An upload manager for S3 storage

        :param config: instance of `ConfigModel`
        :param status_file_path: file for storing upload state. Can be used for e.g. resumable uploads.
        """
        super().__init__()

        self._status_file_path = Path(status_file_path)
        self._s3_options = s3_options
        self._threads = threads

        self._s3_client = init_s3_client(s3_options)

    @override
    def upload_file(self, local_file_path: str | PathLike, s3_object_id: str):
        """
        Upload a single file to the specified object ID
        :param local_file_path: Path to the file to upload
        :param s3_object_id: Remote S3 object ID under which the file should be stored
        """
        self.__log.info(f"Uploading {local_file_path} to {s3_object_id}...")

        filesize = getsize(local_file_path)
        multipart_chunksize = self._s3_options.multipart_chunksize

        chunksize = (
            math.ceil(filesize / MULTIPART_MAX_CHUNKS)
            if filesize / multipart_chunksize > MULTIPART_MAX_CHUNKS
            else multipart_chunksize
        )
        self.__log.debug(
            f"Using a chunksize of: {chunksize / 1024**2}MiB, results in {math.ceil(filesize / chunksize)} chunk(s)"
        )

        config = TransferConfig(
            multipart_threshold=MULTIPART_THRESHOLD,
            multipart_chunksize=chunksize,
            max_concurrency=self._threads,
            use_threads=self._threads > 1,
        )

        transfer = S3Transfer(self._s3_client, config)  # type: ignore[arg-type]
        progress_bar = tqdm(total=filesize, unit="B", unit_scale=True, unit_divisor=1024, smoothing=TQDM_SMOOTHING)
        transfer.upload_file(
            str(local_file_path),
            self._s3_options.bucket,
            s3_object_id,
            callback=lambda bytes_transferred: progress_bar.update(bytes_transferred),
        )

    def _remote_id_exists(self, s3_object_id: str) -> bool:
        """
        Determine if a remote ID already exists
        :param s3_object_id: Remote S3 object ID under which the file should be stored
        """
        exists = True
        try:
            self._s3_client.head_object(Bucket=self._s3_options.bucket, Key=s3_object_id)
        except self._s3_client.exceptions.NoSuchKey:
            exists = False
        except botocore.exceptions.ClientError as error:
            if error.response["Error"]["Code"] in {"403", "404"}:
                # backend can return forbidden instead if user has no ListBucket permission
                exists = False
            else:
                raise error

        return exists

    def _upload_logged_files(self, encrypted_submission, progress_logger, files_to_upload):
        for file_path in files_to_upload:
            if not Path(file_path).exists():
                raise UploadError(f"File {file_path} does not exist")

        for file_path, file_metadata in encrypted_submission.encrypted_files.items():
            logged_state = progress_logger.get_state(file_path, file_metadata)
            self.__log.debug("state for %s: %s", file_path, logged_state)

            s3_object_id = files_to_upload[file_path]

            if (logged_state is None) or not logged_state.get("upload_successful", False):
                self.__log.info(
                    "Uploading file: '%s' -> '%s'",
                    str(file_path),
                    str(s3_object_id),
                )

                try:
                    self.upload_file(file_path, s3_object_id)

                    self.__log.info(f"Upload complete for {str(file_path)}. ")
                    progress_logger.set_state(
                        file_path,
                        file_metadata,
                        state=UploadState(upload_successful=True),
                    )
                except Exception as e:
                    self.__log.error("Upload failed for '%s'", str(file_path))

                    progress_logger.set_state(
                        file_path,
                        file_metadata,
                        state=UploadState(upload_successful=False, errors=[str(e)]),
                    )

                    raise e
            else:
                self.__log.info(
                    "File '%s' already uploaded (at '%s')",
                    str(file_path),
                    str(s3_object_id),
                )

    def _upload_metadata(self, metadata_file_path, metadata_s3_object_id):
        # upload metadata unconditionally
        try:
            self.upload_file(metadata_file_path, metadata_s3_object_id)
            self.__log.info(f"Upload complete for {str(metadata_file_path)}. ")
        except Exception as e:
            self.__log.error("Upload failed for '%s'", str(metadata_file_path))
            raise e

    @override
    def upload(self, encrypted_submission: EncryptedSubmission):
        """
        Upload an encrypted submission
        :param encrypted_submission: The encrypted submission to upload
        """
        progress_logger = FileProgressLogger[UploadState](self._status_file_path)
        metadata_file_path, metadata_s3_object_id = encrypted_submission.get_metadata_file_path_and_object_id()

        if self._remote_id_exists(metadata_s3_object_id):
            raise UploadError("Submission already uploaded. Corrections, additions, and followups require a new tanG.")

        files_to_upload = encrypted_submission.get_encrypted_files_and_object_id()
        files_to_upload[metadata_file_path] = metadata_s3_object_id

        self._upload_logged_files(encrypted_submission, progress_logger, files_to_upload)

        self._upload_metadata(metadata_file_path, metadata_s3_object_id)

    @override
    def archive(self, encrypted_submission: EncryptedSubmission):
        """
        Archive an encrypted submission
        :param encrypted_submission: The encrypted submission to upload
        """
        progress_logger = FileProgressLogger[UploadState](self._status_file_path)
        metadata_file_path, metadata_s3_object_id = encrypted_submission.get_metadata_file_path_and_object_id()

        if self._remote_id_exists(metadata_s3_object_id):
            raise UploadError("Submission already archived.")

        files_to_upload = encrypted_submission.get_encrypted_files_and_object_id()
        files_to_upload[metadata_file_path] = metadata_s3_object_id

        log_file_to_object_id = encrypted_submission.get_log_files_and_object_id()
        overlapping_files = files_to_upload.keys() & log_file_to_object_id.keys()
        if overlapping_files:
            raise ValueError(
                f"Conflict in files specified for archive: {overlapping_files}. This is a bug. Please report this to the developers."
            )
        files_to_upload.update(log_file_to_object_id)

        self._upload_logged_files(encrypted_submission, progress_logger, files_to_upload)

        # do not track upload state for logs, instead just reupload in case of a failure
        for file_path, s3_object_id in encrypted_submission.get_log_files_and_object_id().items():
            try:
                self.upload_file(file_path, s3_object_id)
                self.__log.info(f"Upload complete for {str(file_path)}.")
            except Exception as e:
                self.__log.error("Upload failed for '%s'", str(file_path))
                raise e

        # make a back up copy of metadata before editing
        shutil.copy(metadata_file_path, metadata_file_path.with_suffix(".orig.json"))

        with open(metadata_file_path, mode="r+") as metadata_file:
            metadata = json.load(metadata_file)
            # redact tanG as all zeros
            metadata["submission"]["tanG"] = "".join(["0"] * 64)

            # redact local case ID
            metadata["submission"]["localCaseId"] = ""

            for donor in metadata["donors"]:
                if donor["relation"] == "index":
                    # redact index donorPseudonym (which can be the tanG)
                    donor["donorPseudonym"] = "index"

            metadata_file.seek(0)
            json.dump(metadata, metadata_file, indent=2)
            metadata_file.truncate()

        self._upload_metadata(metadata_file_path, metadata_s3_object_id)

    def _check_for_completed_submission(self, s3_object_id: str) -> bool:
        try:
            return s3_object_id in self._list_keys(self._s3_options.bucket, prefix=str(Path(s3_object_id).parent))
        except Exception as e:
            self.__log.warning(
                "Exception occured during check for completed submission; assuming submission is incomplete.",
                exc_info=e,
            )
            return False

    # https://stackoverflow.com/a/54014862
    def _list_keys(self, bucket_name, prefix="/", delimiter="/", start_after=""):
        s3_paginator = self._s3_client.get_paginator("list_objects_v2")
        prefix = prefix.lstrip(delimiter)
        start_after = (start_after or prefix) if prefix.endswith(delimiter) else start_after
        for page in s3_paginator.paginate(Bucket=bucket_name, Prefix=prefix, StartAfter=start_after):
            for content in page.get("Contents", ()):
                yield content["Key"]
