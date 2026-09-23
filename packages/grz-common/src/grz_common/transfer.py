"""
Common methods for transferring data to and from GRZ buckets.
"""

import datetime
import logging
from collections.abc import Iterator
from contextlib import contextmanager
from typing import TYPE_CHECKING, Any

import boto3
import grz_common.exceptions as grzexc
from boto3 import client as boto3_client  # type: ignore[import-untyped]
from boto3.exceptions import Boto3Error
from botocore.config import Config as Boto3Config
from botocore.exceptions import BotoCoreError, ClientError, NoCredentialsError, PartialCredentialsError

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


S3_CLIENT_ERRORS = (ClientError, BotoCoreError, Boto3Error)
"""What boto3 and botocore raise for a failed request."""

_SETUP_ERROR_CODES = frozenset({"InvalidAccessKeyId", "SignatureDoesNotMatch", "NoSuchBucket"})
"""S3 error codes that only a faulty setup causes.

``AccessDenied`` is not among them: S3 also answers it for a missing object if the credentials
may not list the bucket.
"""


def _chain(error: BaseException) -> Iterator[BaseException]:
    """Yield ``error``, then the exception that ``error`` was raised from or while handling, and so on.

    The chain is the one a traceback shows. It matters for ``S3Transfer.upload_file``, which raises
    ``S3UploadFailedError`` while it handles the ``ClientError`` that carries the error code.
    """
    seen: set[int] = set()
    current: BaseException | None = error
    # a cause can form a cycle, as in ``raise e from e``
    while current is not None and id(current) not in seen:
        seen.add(id(current))
        yield current
        current = current.__cause__ or (None if current.__suppress_context__ else current.__context__)


def _is_faulty_setup(error: BaseException) -> bool:
    """Whether only a faulty setup causes ``error``."""
    if isinstance(error, NoCredentialsError | PartialCredentialsError):
        return True
    return isinstance(error, ClientError) and error.response.get("Error", {}).get("Code") in _SETUP_ERROR_CODES


def s3_error(
    error: Exception, action: str, transfer_error: type[grzexc.TransferError] = grzexc.TransferError
) -> grzexc.GrzError:
    """Classify an error of the S3 client as the failure it stands for.

    The exceptions that ``error`` was raised from or while handling count as well.

    :param error: What the S3 client raised.
    :param action: What failed, such as ``"Upload to s3://bucket/key"``.
    :param transfer_error: The class for a failed transfer.
    :returns: A :class:`ConfigurationError` if only a faulty setup causes ``error``, otherwise a ``transfer_error``.
    """
    faulty_setup = any(_is_faulty_setup(e) for e in _chain(error))
    error_class = grzexc.ConfigurationError if faulty_setup else transfer_error
    return error_class(f"{action} failed: {error}")


@contextmanager
def s3_errors(action: str, transfer_error: type[grzexc.TransferError] = grzexc.TransferError) -> Iterator[None]:
    """Raise an error of the S3 client in the body as the failure it stands for.

    :param action: What the body does, such as ``"Upload to s3://bucket/key"``.
    :param transfer_error: The class for a failed transfer.
    :raises ConfigurationError: If only a faulty setup causes the error, see :func:`s3_error`.
    :raises TransferError: As ``transfer_error``, for any other error of the S3 client.
    """
    try:
        yield
    except S3_CLIENT_ERRORS as e:
        raise s3_error(e, action, transfer_error) from e


def _is_missing_object(error: Exception) -> bool:
    """Whether ``error`` is S3's response for an object that does not exist."""
    if not isinstance(error, ClientError):
        return False
    return error.response.get("Error", {}).get("Code") in {"404", "NoSuchKey", "NotFound"}


def _explain_head_error(s3_client: Any, bucket: str, key: str, error: ClientError) -> ClientError:
    """Return the error that tells why a HEAD request for an S3 object failed.

    S3 answers a HEAD request without a body, so botocore sets the HTTP status as the error code.
    A ``403`` can then mean rejected credentials, and a ``404`` a missing bucket. For these two
    codes, a GET request for the first byte of the object learns the error code.

    :param error: The error of the HEAD request.
    :returns: The error of the GET request, or ``error`` if the GET request succeeds or is not needed.
    """
    if error.response.get("Error", {}).get("Code") not in {"403", "404"}:
        return error
    try:
        response = s3_client.get_object(Bucket=bucket, Key=key, Range="bytes=0-0")
    except ClientError as e:
        return e
    response["Body"].close()
    return error


def head_object(
    s3_client: Any, bucket: str, key: str, transfer_error: type[grzexc.TransferError] = grzexc.DownloadError
) -> dict[str, Any]:
    """Return the ``head_object`` response of an S3 object.

    If S3 answers the HEAD request with ``403`` or ``404`` alone, the error of a GET request
    decides, see :func:`_explain_head_error`.

    :param s3_client: boto3 S3 client.
    :param bucket: Name of the bucket.
    :param key: Key of the object.
    :param transfer_error: The class for a failed transfer.
    :returns: The ``head_object`` response.
    :raises MissingObjectError: If the object does not exist.
    :raises ConfigurationError: If only a faulty setup causes the error, see :func:`s3_error`.
    :raises TransferError: As ``transfer_error``, for any other error of the S3 client.
    """
    with s3_errors(f"Reading s3://{bucket}/{key}", transfer_error):
        try:
            return s3_client.head_object(Bucket=bucket, Key=key)
        except ClientError as e:
            error = _explain_head_error(s3_client, bucket, key, e)
            if _is_missing_object(error):
                raise grzexc.MissingObjectError(f"s3://{bucket}/{key} does not exist") from error
            if error is e:
                raise
            raise error from e


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
    :raises ConfigurationError: If only a faulty setup causes the error of the S3 client.
    :raises DownloadError: For any other error of the S3 client.
    """
    key = f"{submission_id}/metadata/metadata.json"
    try:
        response = head_object(s3_client, bucket, key)
    except grzexc.MissingObjectError as e:
        raise grzexc.MissingSubmissionFileError(f"s3://{bucket}/{key} does not exist") from e
    return response["LastModified"]
