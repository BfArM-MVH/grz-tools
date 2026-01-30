import json
import logging

import botocore
from packaging import version
from pydantic import BaseModel, ConfigDict, Field, ValidationError, field_validator

from ..transfer import init_s3_client

logger = logging.getLogger(__name__)


class VersionFile(BaseModel):
    schema_version: int = Field(1, description="Version of this schema")
    minimal_version: version.Version = Field(..., description="Minimum supported version")
    recommended_version: version.Version = Field(..., description="Recommended version")

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @field_validator("minimal_version", "recommended_version", mode="before")
    @classmethod
    def parse_version(cls, v):
        if isinstance(v, str):
            return version.Version(v)
        return v

    @classmethod
    def from_s3(cls, s3_options, version_file_path) -> "VersionFile":
        """
        Download and validate the version file from S3.

        raises the following errors:
            VersionFileNotFoundError
            VersionFileAccessError
            VersionFileValidationError
        """
        from ..exceptions import (
            VersionFileAccessError,
            VersionFileNotFoundError,
            VersionFileValidationError,
        )

        try:
            s3_client = init_s3_client(s3_options)
            logger.debug(f"Attempting to fetch version file from S3: s3://{s3_options.bucket}/{version_file_path}")

            response = s3_client.get_object(
                Bucket=s3_options.bucket,
                Key=version_file_path,
            )
            version_content = response["Body"].read().decode("utf-8").strip()
            version_data = json.loads(version_content)

            version_file = cls(**version_data)
            logger.info("Successfully retrieved and parsed version file from S3.")
            return version_file

        except botocore.exceptions.ClientError as e:
            error_code = e.response.get("Error", {}).get("Code")
            if error_code == "NoSuchKey":
                msg = (
                    f"Version file not found at "
                    f"s3://{s3_options.bucket}/{version_file_path}. "
                    "Please contact GRZ to resolve this."
                )
                logger.critical(msg)
                raise VersionFileNotFoundError(msg) from e

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
