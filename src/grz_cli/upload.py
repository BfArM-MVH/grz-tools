"""Module for uploading encrypted submissions to a remote storage"""

from __future__ import annotations

import abc
import logging
import multiprocessing
from hashlib import sha256
from multiprocessing import Pool, Queue
from os import PathLike
from os.path import getsize
from pathlib import Path
from traceback import format_exc
from typing import TYPE_CHECKING, override

import boto3  # type: ignore[import-untyped]
from boto3 import client as boto3_client  # type: ignore[import-untyped]
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

    def __init__(self, config: ConfigModel, status_file_path: str | PathLike):
        """
        An upload manager for S3 storage

        :param config: instance of `ConfigModel`
        :param status_file_path: file for storing upload state. Can be used for e.g. resumable uploads.
        """
        super().__init__()

        self._status_file_path = Path(status_file_path)
        self._config = config

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

    def _multipart_upload(self, local_file, s3_object_id):  # noqa: C901
        """
        Upload the file in chunks to S3.

        :param local_file: pathlib.Path()
        :param s3_object_id: string
        :return: sha256 value for uploaded file
        """
        multipart_upload = self._s3_client.create_multipart_upload(
            Bucket=self._config.s3_options.bucket, Key=s3_object_id
        )
        upload_id = multipart_upload["UploadId"]

        # Get the file size for progress bar
        file_size = getsize(local_file)

        # initialize progress bar
        progress_bar = tqdm(
            total=file_size, unit="B", unit_scale=True, unit_divisor=1024
        )
        self._threads = 1

        def _chunk_file(file_path: Path, chunk_size: int, queue: Queue):
            with open(file_path, "rb") as infile:
                part_number = 1
                while True:
                    chunk_data = infile.read(chunk_size)
                    if not chunk_data:
                        break
                    queue.put((part_number, chunk_data))
                    part_number += 1
            # Signal to workers that there's no more data
            for _ in range(self._threads):
                queue.put((None, None))

        def _upload_part(
            bucket_name: str,
            key: str,
            part_number: int,
            chunk_data: bytes,
            upload_id: str,
        ) -> dict | None:
            try:
                response = self._s3_client.upload_part(
                    Bucket=bucket_name,
                    Key=key,
                    PartNumber=part_number,
                    UploadId=upload_id,
                    Body=chunk_data,
                )
                return {"PartNumber": part_number, "ETag": response["ETag"]}
            except Exception as e:
                print(f"Error uploading part {part_number}: {e}")
                return None

        def _upload_part_worker(  # noqa: PLR0913
            queue: Queue,
            progress_bar,
            bucket_name: str,
            key: str,
            upload_id: str,
            parts: list,
        ):
            while True:
                part_number, chunk_data = queue.get()
                if part_number is None:  # Exit signal
                    break

                result = _upload_part(
                    bucket_name, key, part_number, chunk_data, upload_id
                )
                if result:
                    parts.append(result)
                    progress_bar.update(len(chunk_data))

        # Initialize sha256 calculations
        original_sha256 = sha256()

        chunk_size = S3BotoUploadWorker.MULTIPART_CHUNK_SIZE

        # TODO upper limit for queue size / number of parts in queue to avoid memory issues
        queue = Queue(maxsize=self._threads * 4)

        manager = multiprocessing.Manager()
        parts = manager.list()

        worker_config = (
            queue,
            progress_bar,
            self._config.s3_options.bucket,
            s3_object_id,
            upload_id,
            parts,
        )
        with Pool(
            self._threads,
            _upload_part_worker,
            worker_config,
        ) as pool:
            # Read and enqueue chunks
            _chunk_file(local_file, chunk_size, queue)

            pool.close()
            pool.join()

        # Complete the multipart upload
        try:
            log.info("Completing multipart upload...")
            self._s3_client.complete_multipart_upload(
                Bucket=self._config.s3_options.bucket,
                Key=s3_object_id,
                UploadId=upload_id,
                MultipartUpload={"Parts": list(parts)},
            )
        except Exception as e:
            # TODO check whether this is the correct way to handle this
            for i in format_exc().split("\n"):
                log.error(i)
            self._s3_client.abort_multipart_upload(
                Bucket=self._config.s3_options.bucket,
                Key=s3_object_id,
                UploadId=upload_id,
            )
            raise e

        progress_bar.close()  # close progress bar
        return original_sha256.hexdigest()

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
