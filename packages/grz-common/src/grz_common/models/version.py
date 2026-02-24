import json
import logging
from datetime import UTC, datetime
from typing import Annotated, Self

import botocore
from packaging.version import Version
from pydantic import (
    BaseModel,
    BeforeValidator,
    Field,
    PlainSerializer,
    ValidationError,
    field_validator,
    model_validator,
)

from ..exceptions import (
    VersionFileAccessError,
    VersionFileNotFoundError,
    VersionFileValidationError,
)
from ..transfer import init_s3_client
from .s3 import S3Options

logger = logging.getLogger(__name__)


# Custom type converters for packaging.version.Version
def parse_version(v) -> Version:
    """Convert string or Version to Version object."""
    if v is None or isinstance(v, Version):
        return v
    return Version(str(v))


def serialize_version(v: Version) -> str:
    """Convert Version object to string."""
    return str(v)


# Annotated type that Pydantic understands
PydanticVersion = Annotated[
    Version,
    BeforeValidator(parse_version),
    PlainSerializer(serialize_version),
]


class VersionInfo(BaseModel):
    """Version constraint definition for a specific enforcement window."""

    minimal_version: PydanticVersion = Field(..., description="Minimum supported version")
    enforced_from: datetime = Field(..., description="Datetime after which minimal_version becomes compulsory")
    recommended_version: PydanticVersion | None = Field(None, description="Recommended version")
    max_version: PydanticVersion | None = Field(None, description="Highest version tested")

    @field_validator("enforced_from", mode="before")
    @classmethod
    def parse_datetime(cls, v):
        if isinstance(v, datetime):
            # Already a datetime — ensure it is timezone-aware by treating naive datetimes as UTC.
            # This guarantees safe comparison with datetime.now(UTC) regardless of whether the
            # value was constructed programmatically (e.g. in tests) or parsed from a string.
            return v if v.tzinfo is not None else v.replace(tzinfo=UTC)

        # Parse ISO 8601 strings. Strings without a timezone offset (e.g. "2026-03-01T00:00:00")
        # are treated as UTC rather than local time, ensuring consistent behaviour across
        # environments and avoiding ambiguity when comparing against datetime.now(UTC).
        dt = datetime.fromisoformat(str(v))
        return dt if dt.tzinfo is not None else dt.replace(tzinfo=UTC)

    @model_validator(mode="after")
    def check_version_order(self):
        if self.recommended_version is not None and self.recommended_version < self.minimal_version:
            raise ValueError("recommended_version must be >= minimal_version")

        if (
            self.max_version is not None
            and self.recommended_version is not None
            and self.max_version < self.recommended_version
        ):
            raise ValueError("max_version must be >= recommended_version")

        return self


class VersionFile(BaseModel):
    """Version policy definition retrieved from S3."""

    schema_version: int = Field(1, description="Version of this schema")

    # List allows staged future policies
    grzcli_version: list[VersionInfo] = Field(..., description="List of version policies for grz-cli")

    @model_validator(mode="after")
    def ensure_non_empty(self):
        if not self.grzcli_version:
            raise ValueError("At least one version policy must be defined.")
        return self

    @classmethod
    def from_s3(cls, s3_options: S3Options, version_file_key: str) -> Self:
        """Download and validate the version file from S3.

        :param s3_options: Configuration options used to access S3.
        :param version_file_key: Object key of the version file in S3.
        :returns: The validated version file.
        :raises VersionFileNotFoundError: If the version file does not exist in S3.
        :raises VersionFileAccessError: If the version file cannot be accessed.
        :raises VersionFileValidationError: If the version file contents are invalid.
        """
        try:
            s3_client = init_s3_client(s3_options)
            logger.debug(f"Fetching version file: s3://{s3_options.bucket}/{version_file_key}")

            response = s3_client.get_object(
                Bucket=s3_options.bucket,
                Key=version_file_key,
            )

            content = response["Body"].read().decode("utf-8").strip()
            data = json.loads(content)

            version_file = cls(**data)
            logger.info("Version file retrieved and validated successfully.")
            return version_file

        except botocore.exceptions.ClientError as e:
            error_code = e.response.get("Error", {}).get("Code")

            if error_code == "NoSuchKey":
                msg = f"Version file not found at s3://{s3_options.bucket}/{version_file_key}."
                logger.critical(msg)
                raise VersionFileNotFoundError(msg) from e

            msg = f"Unable to access s3://{s3_options.bucket}/{version_file_key} (Error code: {error_code})."
            logger.error(msg, exc_info=e)
            raise VersionFileAccessError(msg) from e

        except (ValidationError, json.JSONDecodeError, ValueError) as e:
            msg = f"Invalid version file format or content: {e}"
            logger.error(msg, exc_info=e)
            raise VersionFileValidationError(msg) from e