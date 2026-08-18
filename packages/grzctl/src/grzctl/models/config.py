import logging
import sys
from contextvars import ContextVar
from pathlib import Path
from typing import Annotated, Any

import yaml
from grz_common.models.base import IgnoringBaseModel, IgnoringBaseSettings
from grz_common.models.identifiers import IdentifiersModel
from grz_common.models.s3 import S3ConnectionBase, S3Options
from pydantic import Field, model_validator
from pydantic.fields import FieldInfo
from pydantic_settings import PydanticBaseSettingsSource

log = logging.getLogger(__name__)

from .db import DbModel
from .pruefbericht import PruefberichtModel

_config_ctx: ContextVar[dict[str, Any] | None] = ContextVar("_config_ctx", default=None)


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


class GrzctlKeyModel(IgnoringBaseModel):
    """Key configuration for grzctl commands."""

    grz_private_key_path: Annotated[str, Field(min_length=1)]
    """Path to the GRZ private key for decryption."""

    grz_public_key_path: Annotated[str | None, Field(default=None)] = None
    """Path to the GRZ public key (used by encrypt wrapper, derived from archives if not set)."""

    submitter_private_key_path: Annotated[str | None, Field(default=None)] = None
    """Path to the submitter's private key (used by encrypt wrapper)."""


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
            field_value, field_key, _value_is_complex = self.get_field_value(field, field_name)
            if field_value is not None:
                d[field_key] = field_value
        return d


class GrzctlConfig(IgnoringBaseSettings):
    """Unified configuration for all grzctl commands."""

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
        config_dict = _config_ctx.get()
        if config_dict is not None:
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
        with open(path) as fd:
            config_dict = yaml.safe_load(fd)
        return cls.from_configuration(config_dict)

    @classmethod
    def from_configuration(cls, configuration: dict[str, Any]) -> "GrzctlConfig":
        """Load config from a dict, letting env vars override dict values."""
        token = _config_ctx.set(configuration)
        try:
            return cls()
        finally:
            _config_ctx.reset(token)

    def resolve_inbox(
        self,
        *,
        submission_id: str | None = None,
        submitter_id: str | None = None,
        bucket: str | None = None,
    ) -> InboxTarget:
        """Resolve an InboxTarget from submission_id, submitter_id, and/or bucket name.

        - If submission_id is provided, the submitter (LE) ID is extracted automatically.
        - If submitter_id is explicitly provided, it takes precedence.
        - If neither is provided, all configured LEs are considered.
        - If bucket is provided, candidates are filtered by bucket name.
        - Auto-resolves if there is exactly one unambiguous candidate.
        """
        target_le = submitter_id
        if not target_le and submission_id:
            target_le = submission_id.split("_", maxsplit=1)[0]

        # filter candidate LEs
        if target_le:
            if target_le not in self.s3.inboxes:
                available = ", ".join(self.s3.inboxes.keys())
                log.error(f"Submitter ID '{target_le}' not found in configuration. Available: {available}")
                sys.exit(1)
            candidate_les = {target_le: self.s3.inboxes[target_le]}
        else:
            candidate_les = self.s3.inboxes

        # collect matching (le_id, bucket_name, config)
        matches = [
            (le, b_name, cfg)
            for le, buckets in candidate_les.items()
            for b_name, cfg in buckets.items()
            if bucket is None or b_name == bucket
        ]

        # handle resolution
        if len(matches) == 1:
            _le, b_name, cfg = matches[0]
            s3_options = S3Options(bucket=b_name, **cfg.model_dump())
            return InboxTarget(
                s3=s3_options,
                **cfg.model_dump(include={"private_key_path", "private_key_passphrase"}),
            )

        if not matches:
            criteria = []
            if target_le:
                criteria.append(f"submitter '{target_le}'")
            if bucket:
                criteria.append(f"bucket '{bucket}'")
            log.error(f"No inbox found matching {' and '.join(criteria) or 'any criteria'}.")
            sys.exit(1)

        # ambiguity handling
        available_choices = ", ".join(f"{le}/{b_name}" for le, b_name, _ in matches)
        if target_le is None and len({le for le, _, _ in matches}) > 1:
            log.error(f"Multiple submitters match ({available_choices}). Please specify --submitter.")
        else:
            log.error(f"Multiple inboxes match ({available_choices}). Please specify --bucket / --inbox-bucket.")
        sys.exit(1)
