"""Module for uploading encrypted submissions to a remote storage"""

from __future__ import annotations

import abc
import logging
import math
import re
from os import PathLike
from os.path import getsize
from pathlib import Path
from typing import TYPE_CHECKING, override

import botocore.handlers  # type: ignore[import-untyped]
from boto3.s3.transfer import S3Transfer, TransferConfig  # type: ignore[import-untyped]
from tqdm.auto import tqdm

from .constants import TQDM_SMOOTHING
from .models.config import ConfigModel
from .progress_logging import FileProgressLogger
from .states import UploadState
from .transfer import init_s3_client

MULTIPART_THRESHOLD = 8 * 1024**2  # 8MiB, boto3 default, largely irrelevant
MULTIPART_MAX_CHUNKS = 1000  # CEPH S3 limit, AWS limit is 10000

if TYPE_CHECKING:
    from .parser import EncryptedSubmission

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
        Upload an encrypted submission

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


class S3BotoUploadWorker(UploadWorker):
    """Implementation of an upload worker using boto3 for S3"""

    __log = log.getChild("S3BotoUploadWorker")

    def __init__(
        self,
        config: ConfigModel,
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
        self._config = config
        self._threads = threads

        self._s3_client = init_s3_client(config)

    @override
    def upload_file(self, local_file_path, s3_object_id):
        """
        Upload a single file to the specified object ID
        :param local_file_path: Path to the file to upload
        :param s3_object_id: Remote S3 object ID under which the file should be stored
        """
        self.__log.info(f"Uploading {local_file_path} to {s3_object_id}...")

        filesize = getsize(local_file_path)
        multipart_chunksize = self._config.s3_options.multipart_chunksize

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

        transfer = S3Transfer(self._s3_client, config)
        progress_bar = tqdm(total=filesize, unit="B", unit_scale=True, unit_divisor=1024, smoothing=TQDM_SMOOTHING)
        transfer.upload_file(
            local_file_path,
            self._config.s3_options.bucket,
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
            self._s3_client.head_object(Bucket=self._config.s3_options.bucket, Key=s3_object_id)
        except self._s3_client.exceptions.NoSuchKey:
            exists = False
        except botocore.exceptions.ClientError as error:
            if error.response["Error"]["Code"] in {"403", "404"}:
                # backend can return forbidden instead if user has no ListBucket permission
                exists = False
            else:
                raise error

        return exists

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

        # finally upload the metadata.json file, unconditionally:
        try:
            self.upload_file(metadata_file_path, metadata_s3_object_id)
            self.__log.info(f"Upload complete for {str(metadata_file_path)}. ")
        except Exception as e:
            self.__log.error("Upload failed for '%s'", str(metadata_file_path))
            raise e

    def _check_for_completed_submission(self, s3_object_id: str) -> bool:
        try:
            return s3_object_id in self._list_keys(
                self._config.s3_options.bucket, prefix=str(Path(s3_object_id).parent)
            )
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
