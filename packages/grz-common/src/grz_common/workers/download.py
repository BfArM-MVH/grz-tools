"""Module for downloading encrypted submissions to local storage"""

from __future__ import annotations

import datetime
import enum
import itertools
import logging
import math
import re
import sys
import threading
from collections import OrderedDict
from collections.abc import Iterable
from operator import attrgetter, itemgetter
from os import PathLike
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

import botocore.handlers
from boto3.s3.transfer import S3Transfer, TransferConfig  # type: ignore[import-untyped]
from grz_pydantic_models.submission.metadata.v1 import File as SubmissionFileMetadata
from pydantic import BaseModel
from tqdm.auto import tqdm

from ..constants import TQDM_DEFAULTS
from ..models.s3 import S3Options
from ..progress import DownloadState, FileProgressLogger
from ..transfer import init_s3_client
from ..utils.concurrency import _run_parallel_with_progress

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

        pool_size = (self._threads * 10) + 1
        self._s3_client = init_s3_client(s3_options, max_pool_connections=pool_size)

    def prepare_download(self, metadata_dir: Path, encrypted_files_dir: Path, log_dir: Path):
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

    def download_metadata(self, submission_id: str, metadata_dir: Path, metadata_file_name: str = "metadata.json"):
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

    def download_file(  # noqa: C901, PLR0913
        self,
        local_file_path: Path,
        s3_object_id: str,
        progress_logger: FileProgressLogger[DownloadState],
        file_metadata: SubmissionFileMetadata,
        submission_id: str,
        filesize: int | None = None,
        pbar_local: Any | None = None,
        pbar_global: Any | None = None,
        global_lock: threading.Lock | None = None,
    ):
        """
        Download a single file from S3 to local storage.

        :param local_file_path: Path to the local target file.
        :param s3_object_id: The S3 object key to download.
        """
        try:
            local_file_path.parent.mkdir(mode=0o770, parents=True, exist_ok=True)

            if filesize is None:
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
                use_threads=self._threads > 1,
            )

            transfer = S3Transfer(self._s3_client, config)  # type: ignore[arg-type]

            if pbar_local:
                pbar_local.reset(total=filesize)
                pbar_local.set_description("DOWNLOAD")
                pbar_local.set_postfix({"file": local_file_path.name})
                pbar_local.refresh()

            def _progress_callback(bytes_transferred):
                if global_lock:
                    with global_lock:
                        if pbar_local:
                            pbar_local.update(bytes_transferred)
                        if pbar_global:
                            pbar_global.update(bytes_transferred)
                else:
                    if pbar_local:
                        pbar_local.update(bytes_transferred)
                    if pbar_global:
                        pbar_global.update(bytes_transferred)

            transfer.download_file(
                self._s3_options.bucket,
                s3_object_id,
                str(local_file_path),
                callback=_progress_callback,
            )

            self.__log.info(f"Download complete for {str(local_file_path)}.")
            progress_logger.set_state(
                local_file_path,
                file_metadata,
                state=DownloadState(download_successful=True, submission_id=submission_id),
            )

        except botocore.exceptions.ClientError as e:
            if e.response.get("Error", {}).get("Code") == "404":
                error_msg = f"File '{s3_object_id}' not found in S3 bucket '{self._s3_options.bucket}'."
                exc = DownloadError(error_msg)
            else:
                error_msg = f"S3 client error for '{s3_object_id}': {e}"
                exc = e  # type: ignore[assignment]
            self.__log.error(error_msg)
            progress_logger.set_state(
                local_file_path,
                file_metadata,
                state=DownloadState(download_successful=False, errors=[str(exc)], submission_id=submission_id),
            )
            raise exc from e
        except Exception as e:
            self.__log.error("Download failed for '%s': %s", str(local_file_path), e)
            progress_logger.set_state(
                local_file_path,
                file_metadata,
                state=DownloadState(download_successful=False, errors=[str(e)], submission_id=submission_id),
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

        self.__log.info("Fetching remote file sizes...")
        remote_sizes = {}

        for _local_file_path, file_metadata in encrypted_submission.encrypted_files.items():
            key = f"{submission_id}/files/{file_metadata.encrypted_file_path()}"
            meta = self._s3_client.head_object(Bucket=self._s3_options.bucket, Key=key)
            remote_sizes[key] = meta["ContentLength"]

        def _get_size(item: tuple[Path, SubmissionFileMetadata]) -> int:
            file_meta = item[1]
            key = f"{submission_id}/files/{file_meta.encrypted_file_path()}"
            return remote_sizes.get(key, 0)

        def _single_download_task(
            item: tuple[Path, SubmissionFileMetadata],
            ui_pos: int,
            lock: threading.Lock,
            pbar_global: tqdm,
        ):
            local_file_path, file_metadata = item
            relative_encrypted_path = file_metadata.encrypted_file_path()
            file_key = f"{submission_id}/files/{relative_encrypted_path}"

            logged_state = progress_logger.get_state(local_file_path, file_metadata)
            filesize = _get_size(item)

            if (
                logged_state
                and logged_state.get("download_successful")
                and logged_state.get("submission_id") == submission_id
            ):
                self.__log.info("File '%s' already downloaded (at '%s'), skipping.", file_key, str(local_file_path))

                with (
                    tqdm(  # type: ignore[arg-type]
                        total=filesize,
                        desc="SKIPPED ",
                        position=ui_pos,
                        file=sys.stderr,
                        postfix={"file": local_file_path.name},
                        **TQDM_DEFAULTS,
                    ) as pbar_local,
                    lock,
                ):
                    pbar_local.update(filesize)
                    pbar_global.update(filesize)
                return

            self.__log.info("Downloading file: '%s' -> '%s'", file_key, str(local_file_path))

            with tqdm(  # type: ignore[arg-type]
                total=filesize,
                desc="DOWNLOAD",
                position=ui_pos,
                file=sys.stderr,
                postfix={"file": local_file_path.name},
                **TQDM_DEFAULTS,
            ) as pbar_local:
                try:
                    self.download_file(
                        local_file_path=local_file_path,
                        s3_object_id=file_key,
                        progress_logger=progress_logger,
                        file_metadata=file_metadata,
                        submission_id=submission_id,
                        filesize=filesize,
                        pbar_local=pbar_local,
                        pbar_global=pbar_global,
                        global_lock=lock,
                    )
                except Exception as e:
                    raise e

        _run_parallel_with_progress(
            items=encrypted_submission.encrypted_files.items(),
            get_size_fn=_get_size,
            worker_fn=_single_download_task,
            threads=self._threads,
        )


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
    total_size_bytes: int


def query_submissions(s3_options: S3Options, show_cleaned: bool) -> list[InboxSubmissionSummary]:
    """Queries the state of all submissions in the configured bucket."""
    s3_client = init_s3_client(s3_options)
    paginator = s3_client.get_paginator("list_objects_v2")

    # cast the pages to a more specific type
    pages = paginator.paginate(Bucket=s3_options.bucket)

    # Use Dict[str, Any] for flexibility
    objects: Iterable[dict[str, Any]] = itertools.chain.from_iterable(
        cast(dict[str, Any], page).get("Contents", []) for page in pages
    )

    # filter safely
    objects = filter(lambda obj: "/" in obj.get("Key", ""), objects)
    objects_sorted = sorted(objects, key=lambda obj: obj.get("Key", ""))

    submission2objects = {
        key: tuple(group)
        for key, group in itertools.groupby(objects_sorted, key=lambda obj: obj.get("Key", "").split("/")[0])
    }

    submissions = []
    for submission_id, submission_objects in submission2objects.items():
        submission_objects_sorted = OrderedDict(
            (o["Key"], o) for o in sorted(submission_objects, key=itemgetter("LastModified"))
        )
        oldest_object = submission_objects_sorted[next(iter(submission_objects_sorted))]
        newest_object = submission_objects_sorted[next(reversed(submission_objects_sorted))]

        total_size_bytes = sum(map(itemgetter("Size"), submission_objects))

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
            total_size_bytes=total_size_bytes,
        )

        if state in {InboxSubmissionState.CLEANING, InboxSubmissionState.CLEANED} and (not show_cleaned):
            # skip listing cleaning/cleaned submissions unless show_cleaned is true
            continue

        submissions.append(submission)

    return sorted(submissions, key=attrgetter("oldest_upload"))
