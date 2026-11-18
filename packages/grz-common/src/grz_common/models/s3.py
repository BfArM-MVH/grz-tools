from typing import Literal

from pydantic import AnyHttpUrl, AnyUrl, SecretStr

from .base import IgnoringBaseModel, IgnoringBaseSettings


class ProxyOptions(IgnoringBaseModel):
    proxy_ca_bundle: str | None = None
    """
    The path to a custom certificate bundle to use when establishing SSL/TLS connections with proxy (optional).
    """

    proxy_client_cert: str | None = None
    """
    The path to a certificate for proxy TLS client authentication  (optional).
    """

    proxy_use_forwarding_for_https: bool = False
    """
    For HTTPS proxies, forward your requests to HTTPS destinations with an absolute URI. We strongly recommend you only use this option with trusted or corporate proxies.
    """


class S3ConnectionBase(IgnoringBaseModel):
    """
    Base configuration for S3 connection details.
    """

    endpoint_url: AnyHttpUrl | None = None
    """
    The URL for the S3 service.
    If undefined, the default AWS endpoint is used.
    """

    access_key: str | None = None
    """
    The access key for the S3 bucket.
    If neither this field nor ``secret`` is set, boto3 looks for credentials itself,
    for example in ``AWS_ACCESS_KEY_ID`` and ``AWS_SECRET_ACCESS_KEY``.
    Setting only one of the two fails.
    """

    secret: SecretStr | None = None
    """
    The secret key for the S3 bucket.
    If neither this field nor ``access_key`` is set, boto3 looks for credentials itself,
    for example in ``AWS_ACCESS_KEY_ID`` and ``AWS_SECRET_ACCESS_KEY``.
    Setting only one of the two fails.
    """

    session_token: SecretStr | None = None
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

    proxy_config: ProxyOptions | None = None
    """Proxy configuration for S3 operations."""

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


class S3Options(S3ConnectionBase):
    """
    S3ConnectionBase with a mandatory bucket name.
    """

    bucket: str
    """
    The name of the S3 bucket.
    """


class S3ConfigModel(IgnoringBaseSettings):
    s3: S3Options
