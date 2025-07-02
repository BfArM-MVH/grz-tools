"""Module for downloading encrypted submissions to local storage"""

from __future__ import annotations

import datetime
import enum
import itertools
import logging
import math
import re
from collections import OrderedDict
from operator import attrgetter, itemgetter
from os import PathLike
from pathlib import Path
from typing import TYPE_CHECKING

import botocore.handlers
from boto3.s3.transfer import S3Transfer, TransferConfig  # type: ignore[import-untyped]
from grz_pydantic_models.submission.metadata.v1 import File as SubmissionFileMetadata
from pydantic import BaseModel
from tqdm.auto import tqdm

from ..constants import TQDM_SMOOTHING
from ..models.s3 import S3Options
from ..progress import DownloadState, FileProgressLogger
from ..transfer import init_s3_client

MULTIPART_THRESHOLD = 8 * 1024 * 1024  # 8MiB, boto3 default
MULTIPART_CHUNKSIZE = 8 * 1024 * 1024  # 8MiB, boto3 default
MULTIPART_MAX_CHUNKS = 1000  # CEPH S3 limit, AWS limit is 10000

if TYPE_CHECKING:
    from .submission import EncryptedSubmission

log = logging.getLogger(__name__)

# see discussion: https://github.com/boto/boto3/discussions/4251 to accept bucket names with ":" in the name
botocore.handlers.VALID_BUCKET = re.compile(r"^[:a-zA-Z0-9.\-_]{1,255}$")


class DownloadError(Exception):
    """Exception raised when an upload fails"""

    pass


class S3BotoDownloadWorker:
    """Implementation of a download worker using boto3 for S3"""

    __log = log.getChild("S3BotoDownloadWorker")

    def __init__(
        self,
        s3_options: S3Options,
        status_file_path: str | PathLike,
        threads: int = 1,
    ):
        """
        A download manager for S3 storage

        :param s3_options: The S3 configuration options
        :param status_file_path: The path to the status file
        :param threads: The number of concurrent download threads
        """
        super().__init__()

        self._status_file_path = Path(status_file_path)
        self._s3_options = s3_options
        self._threads = threads

        self._s3_client = init_s3_client(s3_options)

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
                dir_path.mkdir(parents=False, exist_ok=False)
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

        self.__log.info("Downloading metadata file: '%s'", metadata_key)
        try:
            # Ensure the local target directory exists
            metadata_file_path.parent.mkdir(mode=0o770, parents=True, exist_ok=True)

            self._s3_client.download_file(self._s3_options.bucket, metadata_key, str(metadata_file_path))
            self.__log.info("Metadata download complete.")
        except botocore.exceptions.ClientError as e:
            if e.response.get("Error", {}).get("Code") == "404":
                error_msg = f"Metadata file '{metadata_key}' not found in S3 bucket '{self._s3_options.bucket}'."
                self.__log.error(error_msg)
                raise DownloadError(error_msg) from e
            raise e
        except Exception as e:
            self.__log.error("Download failed for metadata '%s'", metadata_key)
            raise e

    def _download_with_progress(self, local_file_path: str, s3_object_id: str):
        """
        Download a single file from S3 to local storage.

        :param local_file_path: Path to the local target file.
        :param s3_object_id: The S3 object key to download.
        """
        s3_object_meta = self._s3_client.head_object(Bucket=self._s3_options.bucket, Key=s3_object_id)
        filesize = s3_object_meta["ContentLength"]

        chunksize = (
            math.ceil(filesize / MULTIPART_MAX_CHUNKS)
            if filesize / MULTIPART_CHUNKSIZE > MULTIPART_MAX_CHUNKS
            else MULTIPART_CHUNKSIZE
        )
        self.__log.debug(
            f"Using a chunksize of: {chunksize / 1024**2}MiB, results in {math.ceil(filesize / chunksize)} chunks"
        )

        config = TransferConfig(
            multipart_threshold=MULTIPART_THRESHOLD,
            multipart_chunksize=chunksize,
            max_concurrency=self._threads,
        )

        transfer = S3Transfer(self._s3_client, config)  # type: ignore[arg-type]
        with tqdm(
            total=filesize, unit="B", unit_scale=True, unit_divisor=1024, smoothing=TQDM_SMOOTHING
        ) as progress_bar:
            transfer.download_file(
                self._s3_options.bucket,
                s3_object_id,
                local_file_path,
                callback=lambda bytes_transferred: progress_bar.update(bytes_transferred),
            )

    def download_file(
        self,
        local_file_path: Path,
        s3_object_id: str,
        progress_logger: FileProgressLogger[DownloadState],
        file_metadata: SubmissionFileMetadata,
    ):
        """
        Download a single file from S3 to the specified local_file_path.

        :param local_file_path: Path to the local file.
        :param s3_object_id: S3 key of the file to download.
        :param progress_logger: The progress logger instance.
        :param file_metadata: The metadata for the file.
        """
        try:
            local_file_path.parent.mkdir(mode=0o770, parents=True, exist_ok=True)

            self._download_with_progress(str(local_file_path), s3_object_id)

            self.__log.info(f"Download complete for {str(local_file_path)}.")
            progress_logger.set_state(local_file_path, file_metadata, state=DownloadState(download_successful=True))

        except botocore.exceptions.ClientError as e:
            if e.response.get("Error", {}).get("Code") == "404":
                error_msg = f"File '{s3_object_id}' not found in S3 bucket '{self._s3_options.bucket}'."
                exc = DownloadError(error_msg)
            else:
                error_msg = f"S3 client error for '{s3_object_id}': {e}"
                exc = e  # type: ignore[assignment]
            self.__log.error(error_msg)
            progress_logger.set_state(
                local_file_path, file_metadata, state=DownloadState(download_successful=False, errors=[str(exc)])
            )
            raise exc from e
        except Exception as e:
            self.__log.error("Download failed for '%s': %s", str(local_file_path), e)
            progress_logger.set_state(
                local_file_path, file_metadata, state=DownloadState(download_successful=False, errors=[str(e)])
            )
            raise e

    def download(self, submission_id: str, encrypted_submission: EncryptedSubmission):
        """
        Download an encrypted submission.

        This method iterates through the files listed in the submission's metadata,
        constructs their S3 object keys, and downloads them.

        :param submission_id: The ID of the submission, used as a prefix in S3.
        :param encrypted_submission: The encrypted submission to download.
        """
        progress_logger = FileProgressLogger[DownloadState](self._status_file_path)

        for local_file_path, file_metadata in encrypted_submission.encrypted_files.items():
            relative_encrypted_path = file_metadata.encrypted_file_path()
            file_key = f"{submission_id}/files/{relative_encrypted_path}"

            logged_state = progress_logger.get_state(local_file_path, file_metadata)
            if logged_state and logged_state.get("download_successful"):
                self.__log.info(
                    "File '%s' already downloaded (at '%s'), skipping.",
                    file_key,
                    str(local_file_path),
                )
                continue

            self.__log.info("Downloading file: '%s' -> '%s'", file_key, str(local_file_path))
            self.download_file(local_file_path, file_key, progress_logger, file_metadata)


class InboxSubmissionState(enum.StrEnum):
    INCOMPLETE = "incomplete"
    COMPLETE = "complete"
    CLEANING = "cleaning"
    CLEANED = "cleaned"
    ERROR = "error"


class InboxSubmissionSummary(BaseModel):
    """A summary of the state of a submission in an inbox"""

    submission_id: str
    state: InboxSubmissionState
    oldest_upload: datetime.datetime
    newest_upload: datetime.datetime


def query_submissions(s3_options: S3Options, show_cleaned: bool) -> list[InboxSubmissionSummary]:
    """Queries the state of all submissions in the configured bucket."""
    s3_client = init_s3_client(s3_options)
    paginator = s3_client.get_paginator("list_objects_v2")

    objects = itertools.chain.from_iterable(
        page["Contents"] for page in paginator.paginate(Bucket=s3_options.bucket) if "Contents" in page
    )
    objects_sorted = sorted(objects, key=itemgetter("Key"))
    submission2objects = {
        key: tuple(group) for key, group in itertools.groupby(objects_sorted, key=lambda o: o["Key"].split("/")[0])
    }

    submissions = []
    for submission_id, submission_objects in submission2objects.items():
        submission_objects_sorted = OrderedDict(
            (o["Key"], o) for o in sorted(submission_objects, key=itemgetter("LastModified"))
        )
        oldest_object = submission_objects_sorted[next(iter(submission_objects_sorted))]
        newest_object = submission_objects_sorted[next(reversed(submission_objects_sorted))]

        cleaning_key = f"{submission_id}/cleaning"
        cleaned_key = f"{submission_id}/cleaned"
        if (cleaning_key in submission_objects_sorted) and (cleaned_key in submission_objects_sorted):
            log.warning("Submission '{submission_id}' is in an incomplete cleaned state!")
            state = InboxSubmissionState.ERROR
        elif cleaning_key in submission_objects_sorted:
            state = InboxSubmissionState.CLEANING
        elif cleaned_key in submission_objects_sorted:
            state = InboxSubmissionState.CLEANED
        else:
            state = (
                InboxSubmissionState.COMPLETE
                if f"{submission_id}/metadata/metadata.json" in submission_objects_sorted
                else InboxSubmissionState.INCOMPLETE
            )

        submission = InboxSubmissionSummary(
            submission_id=submission_id,
            state=state,
            oldest_upload=oldest_object["LastModified"],
            newest_upload=newest_object["LastModified"],
        )

        if state in {InboxSubmissionState.CLEANING, InboxSubmissionState.CLEANED} and (not show_cleaned):
            # skip listing cleaning/cleaned submissions unless show_cleaned is true
            continue

        submissions.append(submission)

    return sorted(submissions, key=attrgetter("oldest_upload"))
