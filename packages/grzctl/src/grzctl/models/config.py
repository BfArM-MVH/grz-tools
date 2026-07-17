import logging
import sys
from pathlib import Path
from typing import Annotated, Any, ClassVar

from grz_common.models.base import IgnoringBaseModel, IgnoringBaseSettings
from grz_common.models.identifiers import IdentifiersModel
from grz_common.models.s3 import S3ConnectionBase, S3Options
from pydantic import Field, model_validator
from pydantic.fields import FieldInfo
from pydantic_settings import PydanticBaseSettingsSource

log = logging.getLogger(__name__)

from .db import DbModel
from .pruefbericht import PruefberichtModel


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


class ProcessKeyConfigModel(IgnoringBaseModel):
    """Key configuration for the process command."""

    grz_private_key_path: Annotated[str, Field(min_length=1)]
    """Path to the GRZ private key for decryption."""

    consented_archive_public_key_path: Annotated[str, Field(min_length=1)]
    """Path to the public key for re-encryption of consented submissions."""

    non_consented_archive_public_key_path: Annotated[str, Field(min_length=1)]
    """Path to the public key for re-encryption of non-consented submissions."""


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

    @model_validator(mode="after")
    def check_buckets_are_unique(self) -> "ArchivesConfig":
        if self.consented.s3.bucket == self.non_consented.s3.bucket:
            raise ValueError("consented and non-consented buckets must be distinct.")
        return self


class DictConfigSettingsSource(PydanticBaseSettingsSource):
    """A settings source that loads values from a dict (e.g. merged YAML config).

    This source has lower priority than env vars, so environment variables
    can override config file values.
    """

    def __init__(self, settings_cls: type, config_dict: dict[str, Any]):
        super().__init__(settings_cls)
        self.config_dict = config_dict

    def get_field_value(self, field: FieldInfo, field_name: str) -> tuple[Any, str, bool]:
        field_value = self.config_dict.get(field_name)
        return field_value, field_name, False

    def __call__(self) -> dict[str, Any]:
        d: dict[str, Any] = {}
        for field_name, field in self.settings_cls.model_fields.items():
            field_value, field_key, value_is_complex = self.get_field_value(field, field_name)
            if field_value is not None:
                d[field_key] = field_value
        return d


class GrzctlConfig(IgnoringBaseSettings):
    """Unified configuration for all grzctl commands."""

    _config_dict: ClassVar[dict[str, Any] | None] = None

    s3: ProcessS3Options
    """Configuration for S3 inbox connections."""

    archives: ArchivesConfig
    """Configuration for consented and non-consented archives."""

    db: DbModel
    """Database configuration for submission tracking."""

    pruefbericht: PruefberichtModel
    """Configuration for Prüfbericht submission."""

    keys: GrzctlKeyModel
    """Key configuration for encryption/decryption commands."""

    identifiers: IdentifiersModel
    """Identifiers for the GRZ and LE."""

    @classmethod
    def settings_customise_sources(
        cls,
        settings_cls: type,
        init_settings: PydanticBaseSettingsSource,
        env_settings: PydanticBaseSettingsSource,
        dotenv_settings: PydanticBaseSettingsSource,
        file_secret_settings: PydanticBaseSettingsSource,
    ) -> tuple[PydanticBaseSettingsSource, ...]:
        config_dict = getattr(cls, "_config_dict", None)
        if config_dict is not None:
            cls._config_dict = None
            return (
                init_settings,
                env_settings,
                DictConfigSettingsSource(settings_cls, config_dict),
                dotenv_settings,
                file_secret_settings,
            )
        return (init_settings, env_settings, dotenv_settings, file_secret_settings)

    @classmethod
    def from_path(cls, path: str | Path) -> "GrzctlConfig":  # type: ignore[override]
        """Load config from a single YAML file, letting env vars override file values."""
        import yaml

        with open(path) as fd:
            config_dict = yaml.safe_load(fd)
        return cls.from_configuration(config_dict)

    @classmethod
    def from_configuration(cls, configuration: dict[str, Any]) -> "GrzctlConfig":
        """Load config from a dict, letting env vars override dict values."""
        cls._config_dict = configuration
        return cls()

    def resolve_inbox_by_submission_id(self, submission_id: str, bucket: str | None = None) -> InboxTarget:
        """Resolve an inbox by extracting the LE ID from the submission ID.

        The LE ID is everything before the first ``_`` in the submission ID.
        """
        le_id = submission_id.split("_", maxsplit=1)[0]

        if le_id not in self.s3.inboxes:
            available = ", ".join(self.s3.inboxes.keys())
            log.error("Submitter ID '%s' not found in configuration. Available: %s", le_id, available)
            sys.exit(1)

        submitter_inboxes = self.s3.inboxes[le_id]

        if bucket:
            if bucket not in submitter_inboxes:
                log.error("Inbox bucket '%s' not found for '%s'.", bucket, le_id)
                sys.exit(1)
            bucket_name = bucket
            inbox_config = submitter_inboxes[bucket]
        elif len(submitter_inboxes) == 1:
            bucket_name, inbox_config = next(iter(submitter_inboxes.items()))
        else:
            available_buckets = ", ".join(submitter_inboxes.keys())
            log.error(
                "Multiple inboxes found for '%s' (%s). Please specify --inbox-bucket.",
                le_id,
                available_buckets,
            )
            sys.exit(1)

        s3_options = S3Options(bucket=bucket_name, **inbox_config.model_dump())
        return InboxTarget(
            s3=s3_options, **inbox_config.model_dump(include={"private_key_path", "private_key_passphrase"})
        )

    def resolve_inbox_by_bucket(self, bucket: str | None = None) -> S3Options:
        """Resolve S3Options by iterating all inboxes.

        If *bucket* is specified, find it across all LEs.
        Otherwise require exactly one inbox total.
        """
        if bucket:
            for _le_id, buckets in self.s3.inboxes.items():
                if bucket in buckets:
                    inbox_config = buckets[bucket]
                    return S3Options(bucket=bucket, **inbox_config.model_dump())
            log.error("Inbox bucket '%s' not found in any configured LE.", bucket)
            sys.exit(1)

        all_buckets = []
        for le_id, buckets in self.s3.inboxes.items():
            for bucket_name, inbox_config in buckets.items():
                all_buckets.append((le_id, bucket_name, inbox_config))

        if len(all_buckets) == 1:
            _, bucket_name, inbox_config = all_buckets[0]
            return S3Options(bucket=bucket_name, **inbox_config.model_dump())

        bucket_list = ", ".join(f"{le}/{b}" for le, b, _ in all_buckets)
        log.error("Multiple inboxes configured (%s). Please specify --inbox-bucket.", bucket_list)
        sys.exit(1)
