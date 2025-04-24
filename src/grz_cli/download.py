"""Module for uploading encrypted submissions to a remote storage"""

from __future__ import annotations

import itertools
import logging
import math
import re
from operator import itemgetter
from os import PathLike
from pathlib import Path
from typing import TYPE_CHECKING

import boto3  # type: ignore[import-untyped]
import botocore.handlers  # type: ignore[import-untyped]
from boto3 import client as boto3_client  # type: ignore[import-untyped]
from boto3.s3.transfer import S3Transfer, TransferConfig  # type: ignore[import-untyped]
from botocore.config import Config as Boto3Config  # type: ignore[import-untyped]
from tqdm.auto import tqdm

from .constants import TQDM_SMOOTHING
from .models.config import ConfigModel
from .progress_logging import FileProgressLogger
from .states import DownloadState

MULTIPART_THRESHOLD = 8 * 1024 * 1024  # 8MiB, boto3 default
MULTIPART_CHUNKSIZE = 8 * 1024 * 1024  # 8MiB, boto3 default
MULTIPART_MAX_CHUNKS = 1000  # CEPH S3 limit, AWS limit is 10000

if TYPE_CHECKING:
    from .parser import EncryptedSubmission

log = logging.getLogger(__name__)

# see discussion: https://github.com/boto/boto3/discussions/4251 to accept bucket names with ":" in the name
botocore.handlers.VALID_BUCKET = re.compile(r"^[:a-zA-Z0-9.\-_]{1,255}$")  # type: ignore[import-untyped]


# if user specifies empty strings, this might be an issue
def _empty_str_to_none(string: str | None) -> str | None:
    if string == "" or string is None:
        return None
    else:
        return string


def _init_s3_client(config: ConfigModel) -> boto3.session.Session.client:
    # configure proxies if proxy_url is defined
    proxy_url = config.s3_options.proxy_url
    boto3_config = Boto3Config(
        proxies={"http": str(proxy_url), "https": str(proxy_url)} if proxy_url is not None else None,
        request_checksum_calculation=config.s3_options.request_checksum_calculation,
    )

    # Initialize S3 client for downloading
    s3_client: boto3.session.Session.client = boto3_client(
        service_name="s3",
        region_name=_empty_str_to_none(config.s3_options.region_name),
        api_version=_empty_str_to_none(config.s3_options.api_version),
        use_ssl=config.s3_options.use_ssl,
        endpoint_url=_empty_str_to_none(str(config.s3_options.endpoint_url)),
        aws_access_key_id=_empty_str_to_none(config.s3_options.access_key),
        aws_secret_access_key=_empty_str_to_none(config.s3_options.secret),
        aws_session_token=_empty_str_to_none(config.s3_options.session_token),
        config=boto3_config,
    )

    return s3_client


def get_next_submission_id(config: ConfigModel) -> str | None:
    """Finds submission ID for oldest complete submission across all inboxes"""
    logger = log.getChild("SubmissionIDScanner")
    boto3_client = _init_s3_client(config)
    paginator = boto3_client.get_paginator("list_objects_v2")
    metadata_objects_all_inboxes = []

    if (config.internal is None) or (config.internal.inboxes is None):
        logger.warning("internal.inboxes not specified in config, falling back to bucket defined in s3_options")
        inboxes = [config.s3_options.bucket]
    else:
        inboxes = config.internal.inboxes

    for inbox in inboxes:
        logger.info(f"Scanning inbox {inbox}")
        objects = itertools.chain.from_iterable(page["Contents"] for page in paginator.paginate(Bucket=inbox))
        metadata_objects = list(filter(lambda obj: obj["Key"].endswith("metadata/metadata.json"), objects))
        if metadata_objects:
            oldest_metadata_object = sorted(metadata_objects, key=itemgetter("LastModified"))[0]
            metadata_objects_all_inboxes.append(oldest_metadata_object)

    if not metadata_objects_all_inboxes:
        return None

    oldest_metadata_object_all_inboxes = sorted(metadata_objects_all_inboxes, key=itemgetter("LastModified"))[0]
    oldest_submission_id_all_inboxes = oldest_metadata_object_all_inboxes["Key"].split("/")[0]

    return oldest_submission_id_all_inboxes


class DownloadError(Exception):
    """Exception raised when an upload fails"""

    pass


class S3BotoDownloadWorker:
    """Implementation of a download worker using boto3 for S3"""

    __log = log.getChild("S3BotoDownloadWorker")

    def __init__(
        self,
        config: ConfigModel,
        status_file_path: str | PathLike,
        threads: int = 1,
    ):
        """
        An download manager for S3 storage

        :param config: The configuration model
        :param status_file_path: The path to the status file
        """
        super().__init__()

        self._status_file_path = Path(status_file_path)
        self._config = config
        self._threads = threads

        self._init_s3_client()

    def _init_s3_client(self):
        self._s3_client = _init_s3_client(self._config)

    def prepare_download(
        self,
        metadata_dir: Path,
        encrypted_files_dir: Path,
        log_dir: Path,
    ):
        """
        Prepare the download of an encrypted submission

        :param metadata_dir: Path to the metadata directory
        :param encrypted_files_dir: Path to the encrypted_files directory
        :param log_dir: Path to the logs directory
        """
        for dir_path in [metadata_dir, encrypted_files_dir, log_dir]:
            if not dir_path.exists():
                self.__log.debug("Creating directory: %s", dir_path)
                dir_path.mkdir(parents=False, exist_ok=False)  # Create the directories
            else:
                self.__log.debug("Directory exists: %s", dir_path)

    def download_metadata(
        self,
        submission_id: str,
        metadata_dir: Path,
        metadata_file_name: str = "metadata.json",
    ):
        """
        Download the metadata.json

        :param submission_id: submission folder on S3 structure
        :param metadata_dir: Path of the metadir folder
        :param metadata_file_name: name of the metadata.json
        """
        metadata_key = str(Path(submission_id) / metadata_dir.name / metadata_file_name)
        metadata_file_path = metadata_dir / metadata_file_name

        try:
            self.download_file(str(metadata_file_path), metadata_key, 10000)
        except Exception as e:
            self.__log.error("Download failed for '%s'", str(metadata_key))

            raise e

    def download_file(self, local_file_path, s3_object_id, filesize):
        """
        Upload a single file to the specified object ID
        :param local_file_path: Path to the file to upload
        :param s3_object_id: Remote S3 object ID under which the file should be stored
        :param filesize: size of the file
        """
        # self.__log.info(f"Download {s3_object_id} to {local_file_path}...")

        chunksize = (
            math.ceil(filesize / MULTIPART_MAX_CHUNKS)
            if filesize / MULTIPART_CHUNKSIZE > MULTIPART_MAX_CHUNKS
            else MULTIPART_CHUNKSIZE
        )
        self.__log.debug(f"Using a chunksize of: {chunksize / 1024**2}MiB, results in {filesize / chunksize} chunks")

        config = TransferConfig(
            multipart_threshold=MULTIPART_THRESHOLD,
            multipart_chunksize=chunksize,
            max_concurrency=self._threads,
        )

        transfer = S3Transfer(self._s3_client, config)
        progress_bar = tqdm(total=filesize, unit="B", unit_scale=True, unit_divisor=1024, smoothing=TQDM_SMOOTHING)
        transfer.download_file(
            self._config.s3_options.bucket,
            s3_object_id,
            local_file_path,
            callback=lambda bytes_transferred: progress_bar.update(bytes_transferred),
        )

    def download(self, submission_id: str, encrypted_submission: EncryptedSubmission):
        """
        Upload an encrypted submission
        :param encrypted_submission: The encrypted submission to upload
        """
        progress_logger = FileProgressLogger[DownloadState](self._status_file_path)

        encrypted_files_key_prefix = f"{submission_id}/files/"

        response = self._s3_client.list_objects_v2(
            Bucket=self._config.s3_options.bucket, Prefix=encrypted_files_key_prefix
        )

        if "Contents" in response:
            for file in response["Contents"]:
                if file["Size"] == 0:
                    continue
                else:
                    file_key = file["Key"]
                    file_path = Path(file_key).relative_to(encrypted_files_key_prefix)
                    full_path = encrypted_submission.encrypted_files_dir / file_path
                    if full_path not in encrypted_submission.encrypted_files:
                        raise DownloadError(f"File {file_path} not listed in metadata.json")
                    file_metadata = encrypted_submission.encrypted_files[full_path]
                    logged_state = progress_logger.get_state(full_path, file_metadata)

                    if (logged_state is None) or not logged_state.get("download_successful", False):
                        self.__log.info(
                            "Download file: '%s' -> '%s'",
                            file_key,
                            str(full_path),
                        )

                        try:
                            full_path.parent.mkdir(mode=0o770, parents=True, exist_ok=True)
                            self.download_file(str(full_path), file_key, file["Size"])

                            self.__log.info(f"Download complete for {str(full_path)}. ")
                            progress_logger.set_state(
                                full_path,
                                file_metadata,
                                state=DownloadState(download_successful=True),
                            )
                        except Exception as e:
                            self.__log.error("Download failed for '%s'", str(full_path))

                            progress_logger.set_state(
                                full_path,
                                file_metadata,
                                state=DownloadState(download_successful=False, errors=[str(e)]),
                            )

                            raise e
                    else:
                        self.__log.info(
                            "File '%s' already downloaded (at '%s')",
                            file_key,
                            str(full_path),
                        )
