"""
Common methods for transferring data to and from GRZ buckets.
"""

import json
import logging
from typing import TYPE_CHECKING

import boto3
import botocore
from boto3 import client as boto3_client  # type: ignore[import-untyped]
from botocore.config import Config as Boto3Config
from packaging import version
from pydantic import BaseModel, Field, ValidationError, ConfigDict

from .exceptions import (
    VersionFileAccessError,
    VersionFileNotFoundError,
    VersionFileValidationError,
)

logger = logging.getLogger(__name__)


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
        endpoint_url=_empty_str_to_none(str(s3_options.endpoint_url)),
        aws_access_key_id=_empty_str_to_none(s3_options.access_key),
        aws_secret_access_key=_empty_str_to_none(s3_options.secret),
        aws_session_token=_empty_str_to_none(s3_options.session_token),
        config=s3_config,
    )

    return s3_resource


class VersionFile(BaseModel):
    schema_version: int = Field(1, description="Version of this schema")
    minimal_version: version.Version = Field(..., description="Minimum supported version")
    recommended_version: version.Version = Field(..., description="Recommended version")

    # allow arbitrary types like packaging.version.Version
    model_config = ConfigDict(arbitrary_types_allowed=True)


def get_version_info(s3_options, version_file_path) -> VersionFile:
    """
    Download and validate the version file from S3.
    Raises:
        1. FileNotFoundError: If the version file does not exist (GRZ needs to fix).
        2. RuntimeError: For S3 access issues (network, permissions, etc.).
        3. ValueError: For outdated versions or validation issues.
    """
    try:
        s3_client = init_s3_client(s3_options)
        logger.debug(f"Attempting to fetch version file from S3: s3://{s3_options.bucket}/{version_file_path}")

        # try to get the object
        response = s3_client.get_object(Bucket=s3_options.bucket, Key=version_file_path)
        version_content = response["Body"].read().decode("utf-8").strip()
        version_data = json.loads(version_content)

        version_file = VersionFile(**version_data)
        logger.info("Successfully retrieved and parsed version file from S3.")
        return version_file

    except botocore.exceptions.ClientError as e:
        error_code = e.response.get("Error", {}).get("Code")
        if error_code == "NoSuchKey":
            msg = (
                f"Version file not found at s3://{s3_options.bucket}/{version_file_path}. "
                "Please contact GRZ to resolve this."
            )
            logger.critical(msg)
            raise VersionFileNotFoundError(msg) from e
        else:
            msg = (
                f"Unable to access s3://{s3_options.bucket}/{version_file_path} "
                f"(Error code: {error_code}). Possible permission or policy issue."
            )
            logger.error(msg, exc_info=e)
            raise VersionFileAccessError(msg) from e
    except (ValidationError, json.JSONDecodeError, KeyError) as e:
        msg = f"Invalid version file format or content: {e}"
        logger.error(msg, exc_info=e)
        raise VersionFileValidationError(msg) from e
