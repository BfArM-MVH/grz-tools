"""Module for uploading encrypted submissions to a remote storage"""

from __future__ import annotations

import abc
import logging
import re
from os import PathLike
from typing import TYPE_CHECKING

import botocore.handlers

MULTIPART_THRESHOLD = 8 * 1024**2  # 8MiB, boto3 default, largely irrelevant
MULTIPART_MAX_CHUNKS = 1000  # CEPH S3 limit, AWS limit is 10000

if TYPE_CHECKING:
    from ..submission import EncryptedSubmission

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
