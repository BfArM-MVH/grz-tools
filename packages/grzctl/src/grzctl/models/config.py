from typing import Annotated

from grz_common.models.base import IgnoringBaseModel, IgnoringBaseSettings
from grz_common.models.keys import KeyConfigModel
from grz_common.models.s3 import S3ConfigModel, S3Options
from grz_pydantic_models.submission.metadata import GenomicDataCenterId
from pydantic import Field

from .db import DbModel
from .pruefbericht import PruefberichtModel


class ArchiveConfig(S3ConfigModel):
    pass


class DownloadConfig(S3ConfigModel):
    pass


class DecryptConfig(KeyConfigModel):
    pass


class CleanConfig(S3ConfigModel):
    pass


class InboxConfig(IgnoringBaseModel):
    """Configuration for a specific inbox."""

    pass


class ProcessS3Options(S3Options):
    """S3 options for the process command, supporting multiple inboxes."""

    bucket: str | None = None  # type: ignore[assignment]
    """Default bucket name. Required if inboxes is not used or no matching inbox is found."""

    inboxes: dict[str, dict[str, InboxConfig]] | None = None
    """Mapping of submitter IDs to a mapping of bucket names to inbox configurations."""


class ProcessKeyConfigModel(IgnoringBaseSettings):
    """Key configuration for the process command."""

    grz_private_key_path: Annotated[str, Field(min_length=1)]
    """Path to the GRZ private key for decryption."""

    consented_archive_public_key_path: Annotated[str, Field(min_length=1)]
    """Path to the public key for re-encryption of consented submissions."""

    non_consented_archive_public_key_path: Annotated[str, Field(min_length=1)]
    """Path to the public key for re-encryption of non-consented submissions."""


class DetailedQcModel(IgnoringBaseSettings):
    local_storage: Annotated[str, Field(min_length=1)]
    """Path to local storage for detailed QC staging."""

    salt: str | float | int
    """Salt to use for deterministic determination of submissions selected for detailed QC."""

    target_percentage: Annotated[float, Field(ge=0.0, le=100.0)] = 2.0
    """Target percentage of submissions selected for detailed QC per month."""


class ProcessConfig(IgnoringBaseSettings):
    """Configuration for the streaming pipeline process command."""

    s3: ProcessS3Options

    keys: ProcessKeyConfigModel
    """Key configuration for decryption and re-encryption."""

    consented_archive_s3: S3Options
    """S3 configuration for the consented archive destination."""

    non_consented_archive_s3: S3Options
    """S3 configuration for the non-consented archive destination."""

    # Prüfbericht configuration
    pruefbericht: PruefberichtModel
    """Configuration for Prüfbericht submission."""

    # Database configuration
    db: DbModel
    """Database configuration for submission tracking."""

    detailed_qc: DetailedQcModel
    """Configuration for detailed QC selection and staging."""


class PruefberichtConfig(IgnoringBaseSettings):
    pruefbericht: PruefberichtModel


class DbConfig(IgnoringBaseSettings):
    db: DbModel


class ReportIdentifiersConfigModel(IgnoringBaseModel):
    grz: GenomicDataCenterId
    """Id of the GRZ."""


class ReportConfig(DbConfig):
    identifiers: ReportIdentifiersConfigModel


class ListConfig(S3ConfigModel):
    # invalid DbModel will be silently ignored as dict
    db: Annotated[DbModel | dict | None, Field(union_mode="left_to_right")] = None
