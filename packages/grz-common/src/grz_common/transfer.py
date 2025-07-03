"""
Common methods for transferring data to and from GRZ buckets.
"""

from typing import TYPE_CHECKING

import boto3
from boto3 import client as boto3_client  # type: ignore[import-untyped]
from botocore.config import Config as Boto3Config  # type: ignore[import-untyped]

if TYPE_CHECKING:
    from types_boto3_s3 import S3Client
    from types_boto3_s3.service_resource import S3ServiceResource
else:
    # avoid undefined objects when not type checking
    S3Client = object
    S3ServiceResource = object

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
    s3_config = Boto3Config(
        proxies={"http": str(proxy_url), "https": str(proxy_url)} if proxy_url is not None else None,
        request_checksum_calculation=s3_options.request_checksum_calculation,
    )

    # Initialize S3 client for uploading
    s3_client: S3Client = boto3_client(
        service_name="s3",
        region_name=_empty_str_to_none(s3_options.region_name),
        api_version=_empty_str_to_none(s3_options.api_version),
        use_ssl=s3_options.use_ssl,
        endpoint_url=_empty_str_to_none(str(s3_options.endpoint_url)),
        aws_access_key_id=_empty_str_to_none(s3_options.access_key),
        aws_secret_access_key=_empty_str_to_none(s3_options.secret),
        aws_session_token=_empty_str_to_none(s3_options.session_token),
        config=s3_config,
    )

    return s3_client


def init_s3_resource(s3_options: S3Options) -> S3ServiceResource:
    """Create a boto3 Resource from a grz-cli configuration."""
    proxy_url = s3_options.proxy_url
    s3_config = Boto3Config(
        proxies={"http": str(proxy_url), "https": str(proxy_url)} if proxy_url is not None else None,
        request_checksum_calculation=s3_options.request_checksum_calculation,
    )
    s3_resource = boto3.resource(
        service_name="s3",
        region_name=_empty_str_to_none(s3_options.region_name),
        api_version=_empty_str_to_none(s3_options.api_version),
        use_ssl=s3_options.use_ssl,
        endpoint_url=_empty_str_to_none(str(s3_options.endpoint_url)),
        aws_access_key_id=_empty_str_to_none(s3_options.access_key),
        aws_secret_access_key=_empty_str_to_none(s3_options.secret),
        aws_session_token=_empty_str_to_none(s3_options.session_token),
        config=s3_config,
    )

    return s3_resource
