"""
Common methods for transferring data to and from GRZ buckets.
"""

import datetime
import logging
from collections.abc import Iterator
from contextlib import contextmanager
from http import HTTPStatus
from typing import TYPE_CHECKING, Any

import boto3  # type: ignore[import-untyped]
import grz_common.exceptions as grzexc
from boto3 import client as boto3_client  # type: ignore[import-untyped]
from boto3.exceptions import Boto3Error  # type: ignore[import-untyped]
from botocore.config import Config as Boto3Config
from botocore.exceptions import (
    BotoCoreError,
    ClientError,
    NoCredentialsError,
    ParamValidationError,
    PartialCredentialsError,
)

logger = logging.getLogger(__name__)


if TYPE_CHECKING:
    from types_boto3_s3 import S3Client
    from types_boto3_s3.service_resource import S3ServiceResource
else:
    # avoid undefined objects when not type checking
    S3Client = object
    S3ServiceResource = object

from .models.base import get_secret_value
from .models.s3 import S3Options


def _empty_str_to_none(string: str | None) -> str | None:
    # if user specifies empty strings, this might be an issue
    if string == "" or string is None:
        return None
    else:
        return string


def init_s3_client(s3_options: S3Options) -> S3Client:
    """Create a boto3 Client from a grz-cli configuration."""
    # configure proxies if proxy_url is defined
    proxy_url = s3_options.proxy_url
    proxies_config = s3_options.proxy_config.model_dump(exclude_none=True) if s3_options.proxy_config else None
    s3_config = Boto3Config(
        proxies={"http": str(proxy_url), "https": str(proxy_url)} if proxy_url is not None else None,
        proxies_config=proxies_config,  # type: ignore
        request_checksum_calculation=s3_options.request_checksum_calculation,
    )

    # Initialize S3 client for uploading
    s3_client: S3Client = boto3_client(
        service_name="s3",
        region_name=_empty_str_to_none(s3_options.region_name),
        api_version=_empty_str_to_none(s3_options.api_version),
        use_ssl=s3_options.use_ssl,
        endpoint_url=_empty_str_to_none(str(s3_options.endpoint_url)) if s3_options.endpoint_url else None,
        aws_access_key_id=_empty_str_to_none(s3_options.access_key),
        aws_secret_access_key=_empty_str_to_none(get_secret_value(s3_options.secret)),
        aws_session_token=_empty_str_to_none(get_secret_value(s3_options.session_token)),
        config=s3_config,
    )

    return s3_client


def init_s3_resource(s3_options: S3Options) -> S3ServiceResource:
    """Create a boto3 Resource from a grz-cli configuration."""
    proxy_url = s3_options.proxy_url
    proxies_config = s3_options.proxy_config.model_dump(exclude_none=True) if s3_options.proxy_config else None
    s3_config = Boto3Config(
        proxies={"http": str(proxy_url), "https": str(proxy_url)} if proxy_url is not None else None,
        proxies_config=proxies_config,  # type: ignore
        request_checksum_calculation=s3_options.request_checksum_calculation,
    )
    s3_resource = boto3.resource(
        service_name="s3",
        region_name=_empty_str_to_none(s3_options.region_name),
        api_version=_empty_str_to_none(s3_options.api_version),
        use_ssl=s3_options.use_ssl,
        endpoint_url=_empty_str_to_none(str(s3_options.endpoint_url)) if s3_options.endpoint_url else None,
        aws_access_key_id=_empty_str_to_none(s3_options.access_key),
        aws_secret_access_key=_empty_str_to_none(get_secret_value(s3_options.secret)),
        aws_session_token=_empty_str_to_none(get_secret_value(s3_options.session_token)),
        config=s3_config,
    )

    return s3_resource


_SETUP_ERROR_CODES = frozenset({"InvalidAccessKeyId", "SignatureDoesNotMatch", "NoSuchBucket"})
"""S3 error codes that only a faulty setup causes.

``AccessDenied`` is not among them: S3 also answers it for a missing object if the credentials
may not list the bucket.
"""


@contextmanager
def s3_errors(action: str, transfer_error: type[grzexc.TransferError] = grzexc.TransferError) -> Iterator[None]:
    """Raise an error of the S3 client in the body as the failure it stands for.

    :param action: What the body does, such as ``"Upload to s3://bucket/key"``.
    :param transfer_error: The class for a failed transfer.
    :raises ConfigurationError: If only a faulty setup causes the error, such as rejected credentials.
    :raises TransferError: As ``transfer_error``, for any other error of the S3 client.
    """
    try:
        yield
    except ClientError as e:
        error_class = grzexc.ConfigurationError if e.response["Error"]["Code"] in _SETUP_ERROR_CODES else transfer_error
        raise error_class(f"{action} failed: {e}") from e
    except (NoCredentialsError, PartialCredentialsError, ParamValidationError) as e:
        # botocore raises ParamValidationError for a bucket name from the config before it sends a request
        raise grzexc.ConfigurationError(f"{action} failed: {e}") from e
    except (BotoCoreError, Boto3Error) as e:
        raise transfer_error(f"{action} failed: {e}") from e


def head_object(
    s3_client: Any,
    bucket: str,
    key: str,
    transfer_error: type[grzexc.TransferError] = grzexc.DownloadError,
    missing_error: type[grzexc.GrzError] = grzexc.MissingObjectError,
) -> dict[str, Any]:
    """Return the ``head_object`` response of an S3 object.

    :param s3_client: boto3 S3 client.
    :param bucket: Name of the bucket.
    :param key: Key of the object.
    :param transfer_error: The class for a failed transfer.
    :param missing_error: The class for a missing object.
    :returns: The ``head_object`` response.
    :raises MissingObjectError: As ``missing_error``, if the object does not exist.
    :raises ConfigurationError: If only a faulty setup causes the error, see :func:`s3_errors`.
    :raises TransferError: As ``transfer_error``, for any other error of the S3 client.
    """
    with s3_errors(f"Reading s3://{bucket}/{key}", transfer_error):
        try:
            return s3_client.head_object(Bucket=bucket, Key=key)
        except ClientError as e:
            if e.response["ResponseMetadata"]["HTTPStatusCode"] not in {HTTPStatus.FORBIDDEN, HTTPStatus.NOT_FOUND}:
                raise
            head_error = e
        # a HEAD answer has no body, so its error code is only the HTTP status. A GET answers with the real code.
        try:
            s3_client.get_object(Bucket=bucket, Key=key, Range="bytes=0-0")["Body"].close()
        except ClientError as e:
            if e.response["Error"]["Code"] == "NoSuchKey":
                raise missing_error(f"s3://{bucket}/{key} does not exist") from e
            raise
        raise head_error


def raise_if_cleaned(s3_client: Any, bucket: str, submission_id: str) -> None:
    """Raise if ``grzctl clean`` has started on the submission in the inbox.

    Cleaning puts a ``<submission_id>/cleaning`` object before it deletes the files and empties
    the metadata.json, and replaces it with ``<submission_id>/cleaned`` at the end. So an emptied
    metadata.json always sits next to one of the two markers. An empty metadata.json without a
    marker was uploaded like that, so it is no concern of this check.

    :param s3_client: boto3 S3 client pointed at the inbox bucket.
    :param bucket: Name of the inbox bucket.
    :param submission_id: Submission identifier (the top-level S3 prefix).
    :raises SubmissionCleanedError: If a marker of ``grzctl clean`` exists.
    :raises ConfigurationError: If only a faulty setup causes the error of the S3 client.
    :raises DownloadError: For any other error of the S3 client.
    """
    for key in (f"{submission_id}/cleaning", f"{submission_id}/cleaned"):
        with s3_errors(f"Reading s3://{bucket}/{key}", grzexc.DownloadError):
            try:
                s3_client.head_object(Bucket=bucket, Key=key)
            except ClientError as e:
                # the caller has read the metadata.json, so the bucket exists and a 404 means a missing marker
                if e.response["ResponseMetadata"]["HTTPStatusCode"] != HTTPStatus.NOT_FOUND:
                    raise
            else:
                raise grzexc.SubmissionCleanedError(
                    f"grzctl clean has started on {submission_id}: s3://{bucket} holds {key}"
                )


def get_metadata_upload_timestamp(s3_client: S3Client, bucket: str, submission_id: str) -> datetime.datetime:
    """Return the S3 last-modified timestamp of a submission's ``metadata/metadata.json`` object.

    This is the authoritative "received at the inbox" timestamp: it cannot be forged by
    the submitter (unlike ``submission.submissionDate`` inside the JSON itself) and is
    only meaningful while the object still lives in the inbox bucket. Do **not** call
    this against the archive bucket: the archive's ``LastModified`` reflects the time
    of archival, not the time of submission.

    :param s3_client: boto3 S3 client pointed at the inbox bucket.
    :param bucket: Name of the inbox bucket.
    :param submission_id: Submission identifier (the top-level S3 prefix).
    :returns: ``LastModified`` (timezone-aware ``datetime``) for
        ``<submission_id>/metadata/metadata.json``. Callers that only need the date
        portion should call ``.date()`` themselves.
    :raises MissingSubmissionFileError: If the inbox lacks the metadata.
    :raises SubmissionCleanedError: If ``grzctl clean`` has started on the submission, see :func:`raise_if_cleaned`.
    :raises ConfigurationError: If only a faulty setup causes the error of the S3 client.
    :raises DownloadError: For any other error of the S3 client.
    """
    key = f"{submission_id}/metadata/metadata.json"
    response = head_object(s3_client, bucket, key, missing_error=grzexc.MissingSubmissionFileError)
    # Check if the submission is (being) cleaned from the inbox. If yes, the metadata.json's timestamp is invalid.
    raise_if_cleaned(s3_client, bucket, submission_id)
    return response["LastModified"]
