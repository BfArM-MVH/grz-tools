from typing import Annotated

from grz_common.models.base import IgnoringBaseModel, IgnoringBaseSettings
from grz_common.models.keys import KeyConfigModel
from grz_common.models.s3 import S3ConfigModel, S3ConnectionBase, S3Options
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


class InboxConfig(S3ConnectionBase):
    """
    Configuration for a specific inbox.
    Includes connection details and the private key needed to decrypt its contents.
    """

    private_key_path: Annotated[str, Field(min_length=1)]
    """Path to the GRZ private key used to decrypt files from this inbox."""

    private_key_passphrase: Annotated[str | None, Field(default=None)] = None
    """Passphrase to the GRZ private key used to decrypt files from this inbox."""


class InboxTarget(IgnoringBaseModel):
    """
    Fully resolved source configuration.
    Encapsulates everything needed to read and decrypt from a specific inbox.
    """

    s3: S3Options
    """Fully resolved S3 options, including the bucket name."""

    private_key_path: Annotated[str, Field(min_length=1)]
    """Path to the GRZ private key used to decrypt files from this inbox."""

    private_key_passphrase: Annotated[str | None, Field(default=None)] = None
    """Passphrase to the GRZ private key used to decrypt files from this inbox."""


class ProcessS3Options(IgnoringBaseModel):
    """
    Root S3 configuration for the process command.
    Enforces that at least one LE and one inbox are defined.
    """

    inboxes: Annotated[
        dict[str, Annotated[dict[str, InboxConfig], Field(min_length=1)]],
        Field(min_length=1),
    ]
    """
    Mapping: LE-Id -> BucketName -> InboxConfig.
    """


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

    salt: str
    """Salt to use for deterministic determination of submissions selected for detailed QC."""

    target_percentage: Annotated[float, Field(ge=0.0, le=100.0)] = 2.0
    """Target percentage of submissions selected for detailed QC per month."""


class ArchiveTarget(IgnoringBaseModel):
    """Encapsulates everything needed to write to a specific archive."""

    s3: S3Options
    """S3 connection details and bucket for this archive."""

    public_key_path: Annotated[str, Field(min_length=1)]
    """Path to the public key for re-encryption of files destined for this archive."""


class ArchivesConfig(IgnoringBaseModel):
    """Configuration for consented and non-consented archives."""

    consented: ArchiveTarget
    """Target definition for consented submissions."""

    non_consented: ArchiveTarget
    """Target definition for non-consented submissions."""


class ProcessConfig(IgnoringBaseSettings):
    """Configuration for the streaming pipeline process command."""

    s3: ProcessS3Options
    """Configuration for S3 connections."""

    archives: ArchivesConfig
    """Configuration for consented and non-consented archives."""

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
