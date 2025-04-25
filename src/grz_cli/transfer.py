"""
Common methods for transfering data to and from GRZ buckets.
"""

import boto3  # type: ignore[import-untyped]
from boto3 import client as boto3_client  # type: ignore[import-untyped]
from botocore.config import Config as Boto3Config  # type: ignore[import-untyped]

from .models.config import ConfigModel


def _empty_str_to_none(string: str | None) -> str | None:
    # if user specifies empty strings, this might be an issue
    if string == "" or string is None:
        return None
    else:
        return string


def init_s3_client(config: ConfigModel) -> boto3.session.Session.client:
    """Create a boto3 Client from a grz-cli configuration."""
    # configure proxies if proxy_url is defined
    proxy_url = config.s3_options.proxy_url
    s3_config = Boto3Config(
        proxies={"http": str(proxy_url), "https": str(proxy_url)} if proxy_url is not None else None,
        request_checksum_calculation=config.s3_options.request_checksum_calculation,
    )

    # Initialize S3 client for uploading
    s3_client: boto3.session.Session.client = boto3_client(
        service_name="s3",
        region_name=_empty_str_to_none(config.s3_options.region_name),
        api_version=_empty_str_to_none(config.s3_options.api_version),
        use_ssl=config.s3_options.use_ssl,
        endpoint_url=_empty_str_to_none(str(config.s3_options.endpoint_url)),
        aws_access_key_id=_empty_str_to_none(config.s3_options.access_key),
        aws_secret_access_key=_empty_str_to_none(config.s3_options.secret),
        aws_session_token=_empty_str_to_none(config.s3_options.session_token),
        config=s3_config,
    )

    return s3_client
