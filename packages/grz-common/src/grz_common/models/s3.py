from typing import Literal

from pydantic import AnyHttpUrl, AnyUrl

from .base import IgnoringBaseModel, IgnoringBaseSettings


class S3Options(IgnoringBaseModel):
    endpoint_url: AnyHttpUrl
    """
    The URL for the S3 service.
    """

    bucket: str
    """
    The name of the S3 bucket.
    """

    access_key: str | None = None
    """
    The access key for the S3 bucket.
    If undefined, it is read from the AWS_ACCESS_KEY_ID environment variable.
    """

    secret: str | None = None
    """
    The secret key for the S3 bucket.
    If undefined, it is read from the AWS_SECRET_ACCESS_KEY environment variable.
    """

    session_token: str | None = None
    """
    The session token for temporary credentials (optional).
    """

    region_name: str | None = None
    """
    The region name for the S3 bucket.
    """

    api_version: str | None = None
    """
    The S3 API version.
    """

    use_ssl: bool = True
    """
    Whether to use SSL for S3 operations.
    """

    proxy_url: AnyUrl | None = None
    """
    The proxy URL for S3 operations (optional).
    """

    request_checksum_calculation: Literal["when_supported", "when_required"] | None = None
    """
    Whether to calculate checksums for S3 request payloads (optional).
    Valid values are:

    * ``when_supported`` -- When set, a checksum will be calculated for
      all request payloads of operations modeled with the ``httpChecksum``
      trait where ``requestChecksumRequired`` is ``true`` or a
      ``requestAlgorithmMember`` is modeled.

    * ``when_required`` -- When set, a checksum will only be calculated
      for request payloads of operations modeled with the ``httpChecksum``
      trait where ``requestChecksumRequired`` is ``true`` or where a
      ``requestAlgorithmMember`` is modeled and supplied.

    Defaults to None.
    """

    multipart_chunksize: int = 256 * 1024**2
    """
    The size of the chunks to use for multipart uploads in bytes.
    """


class S3ConfigModel(IgnoringBaseSettings):
    s3: S3Options
