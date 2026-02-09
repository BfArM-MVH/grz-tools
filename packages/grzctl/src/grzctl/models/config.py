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


class ProcessKeyConfigModel(IgnoringBaseSettings):
    """Key configuration for the process command."""

    grz_private_key_path: str
    """Path to the GRZ private key for decryption."""

    consented_archive_public_key_path: str
    """Path to the public key for re-encryption of consented submissions."""

    non_consented_archive_public_key_path: str
    """Path to the public key for re-encryption of non-consented submissions."""


class ProcessConfig(S3ConfigModel):
    """Configuration for the streaming pipeline process command."""

    keys: ProcessKeyConfigModel
    """Key configuration for decryption and re-encryption."""

    consented_archive_s3: S3Options
    """S3 configuration for the consented archive destination."""

    non_consented_archive_s3: S3Options
    """S3 configuration for the non-consented archive destination."""

    # Prüfbericht configuration
    pruefbericht: PruefberichtModel
    """Configuration for Prüfbericht submission (optional)."""

    # Database configuration
    db: DbModel
    """Database configuration for submission tracking (optional)."""


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
