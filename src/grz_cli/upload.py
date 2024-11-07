"""Module for uploading encrypted submissions to a remote storage"""

from __future__ import annotations

import abc
import logging
import multiprocessing
from os import PathLike
from os.path import getsize
from pathlib import Path
from typing import TYPE_CHECKING, override

import boto3  # type: ignore[import-untyped]
from boto3 import client as boto3_client  # type: ignore[import-untyped]
from boto3.s3.transfer import S3Transfer, TransferConfig  # type: ignore[import-untyped]
from botocore.config import Config as Boto3Config  # type: ignore[import-untyped]
from tqdm.auto import tqdm

from .models.config import ConfigModel

if TYPE_CHECKING:
    from .submission import EncryptedSubmission

log = logging.getLogger(__name__)


class UploadError(Exception):
    """Exception raised when an upload fails"""

    pass


class UploadWorker(metaclass=abc.ABCMeta):
    """Worker baseclass for uploading encrypted submissions"""

    def upload(self, encrypted_submission: EncryptedSubmission):
        """
        Upload an encrypted submission

        :param encrypted_submission: The encrypted submission to upload
        :raises UploadError: when the upload failed
        """
        submission_id = encrypted_submission.metadata.index_case_id
        for (
            local_file_path,
            file_metadata,
        ) in encrypted_submission.encrypted_files.items():
            relative_file_path = file_metadata.file_path
            s3_object_id = Path(submission_id) / "files" / relative_file_path

            try:
                self.upload_file(local_file_path, str(s3_object_id))
            except Exception as e:
                raise UploadError(
                    f"Failed to upload {local_file_path} (object id: {s3_object_id})"
                ) from e

        # upload the metadata file
        s3_object_id = Path(submission_id) / "metadata" / "metadata.json"
        try:
            self.upload_file(encrypted_submission.metadata.file_path, str(s3_object_id))
        except Exception as e:
            raise UploadError(
                f"Failed to upload metadata: {encrypted_submission.metadata.file_path} (object id: {s3_object_id})"
            ) from e

    @abc.abstractmethod
    def upload_file(self, local_file_path: str | PathLike, s3_object_id: str):
        """
        Upload a single file to the specified object ID
        :param local_file_path: Path to the file to upload
        :param s3_object_id: Remote S3 object ID under which the file should be stored
        """
        raise NotImplementedError()


class S3BotoUploadWorker(UploadWorker):
    """Implementation of an upload operations using boto3 for S3"""

    __log = log.getChild("S3BotoUploadWorker")

    MULTIPART_CHUNK_SIZE = 64 * 1024 * 1024  # 64 MB
    MAX_SINGLEPART_UPLOAD_SIZE = 5 * 1024 * 1024  # 5 MB

    def __init__(
        self,
        config: ConfigModel,
        status_file_path: str | PathLike,
        threads: int | None = None,
    ):
        """
        An upload manager for S3 storage

        :param config: instance of `ConfigModel`
        :param status_file_path: file for storing upload state. Can be used for e.g. resumable uploads.
        """
        super().__init__()

        self._status_file_path = Path(status_file_path)
        self._config = config
        self._threads = threads or multiprocessing.cpu_count()

        self._init_s3_client()

    def _init_s3_client(self):
        # if user specifies empty strings, this might be an issue
        def empty_str_to_none(string: str | None) -> str | None:
            if string == "" or string is None:
                return None
            else:
                return string

        # configure proxies if proxy_url is defined
        proxy_url = empty_str_to_none(self._config.s3_options.proxy_url)
        if proxy_url is not None:
            config = Boto3Config(proxies={"http": proxy_url, "https": proxy_url})
        else:
            config = None

        # Initialize S3 client for uploading
        self._s3_client: boto3.session.Session.client = boto3_client(
            service_name="s3",
            region_name=empty_str_to_none(self._config.s3_options.region_name),
            api_version=empty_str_to_none(self._config.s3_options.api_version),
            use_ssl=empty_str_to_none(str(self._config.s3_options.use_ssl).lower()),
            endpoint_url=empty_str_to_none(str(self._config.s3_options.endpoint_url)),
            aws_access_key_id=empty_str_to_none(self._config.s3_options.access_key),
            aws_secret_access_key=empty_str_to_none(self._config.s3_options.secret),
            aws_session_token=empty_str_to_none(self._config.s3_options.session_token),
            config=config,
        )

    # def show_information(self):
    #     self.__log.info(f"total files in metafile: {self.__file_total}")
    #     self.__log.info(f"uploaded files: {self.__file_done}")
    #     self.__log.info(f"failed files: {self.__file_failed}")
    #     self.__log.info(
    #         f"already finished files before current upload: {self.__file_prefinished}"
    #     )

    def _multipart_upload(self, local_file, s3_object_id):
        """
        Upload the file in chunks to S3.

        :param local_file: pathlib.Path()
        :param s3_object_id: string
        :return: sha256 value for uploaded file
        """
        config = TransferConfig(
            multipart_threshold=5 * 1024 * 1024,
            max_concurrency=self._threads,
        )
        transfer = S3Transfer(self._s3_client, config)
        progress_bar = tqdm(
            total=getsize(local_file), unit="B", unit_scale=True, unit_divisor=1024
        )
        transfer.upload_file(
            local_file,
            self._config.s3_options.bucket,
            s3_object_id,
            callback=lambda bytes_transferred: progress_bar.update(bytes_transferred),
        )

    def _upload(self, local_file, s3_object_id):
        """
        Upload the file to S3.

        :param local_file: pathlib.Path()
        :param s3_object_id: string
        :return: sha256 values for original file
        """
        try:
            with open(local_file, "rb") as fd:
                # calculate sha256sum
                # original_sha256 = sha256(data)

                # Upload data
                self._s3_client.put_object(
                    Bucket=self._config.s3_options.bucket, Key=s3_object_id, Body=fd
                )
            # return original_sha256.hexdigest()
        except Exception as e:
            raise e

    @override
    def upload_file(self, local_file_path, s3_object_id):
        """
        Upload a single file to the specified object ID
        :param local_file_path: Path to the file to upload
        :param s3_object_id: Remote S3 object ID under which the file should be stored
        """
        self.__log.info(f"Uploading {local_file_path} to {s3_object_id}...")

        # Get the file size to decide whether to use multipart upload
        file_size = getsize(local_file_path)
        if file_size > S3BotoUploadWorker.MAX_SINGLEPART_UPLOAD_SIZE:
            # do multipart upload
            _sha256sums = self._multipart_upload(local_file_path, s3_object_id)
        else:
            _sha256sums = self._upload(local_file_path, s3_object_id)
        # TODO: check sha256sums
