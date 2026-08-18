import logging
import sys
from contextvars import ContextVar
from pathlib import Path
from typing import Annotated, Any

import yaml
from grz_common.models.base import IgnoringBaseModel, IgnoringBaseSettings
from grz_common.models.identifiers import IdentifiersModel
from grz_common.models.s3 import S3ConnectionBase, S3Options
from pydantic import Field, PrivateAttr, model_validator
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

    bucket: Annotated[str | None, Field(default=None)] = None
    """S3 bucket name. Defaults to the inbox name key if not set."""

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


class LeistungserbringerEntry(IgnoringBaseModel):
    """A single Leistungserbringer (submitter) with its inbox configurations."""

    alias: Annotated[str | None, Field(default=None)] = None
    """Human-friendly alias. Defaults to the LE id (dict key) if not set."""

    inbox_buckets: Annotated[dict[str, InboxConfig], Field(min_length=1)]
    """Mapping: InboxName -> InboxConfig."""


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

    leistungserbringer: Annotated[dict[str, LeistungserbringerEntry], Field(min_length=1)]
    """Mapping: LE-Id -> LeistungserbringerEntry."""

    _le_by_id: dict[str, LeistungserbringerEntry] = PrivateAttr(default_factory=dict)
    _le_by_alias: dict[str, LeistungserbringerEntry] = PrivateAttr(default_factory=dict)

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

    @model_validator(mode="after")
    def build_le_lookups(self) -> "GrzctlConfig":
        for le_id, entry in self.leistungserbringer.items():
            if entry.alias is None:
                entry.alias = le_id
            if le_id in self._le_by_id:
                raise ValueError(f"Duplicate LE id: '{le_id}'")
            if entry.alias in self._le_by_alias:
                raise ValueError(f"Duplicate LE alias: '{entry.alias}'")
            self._le_by_id[le_id] = entry
            self._le_by_alias[entry.alias] = entry
        return self

    @staticmethod
    def _describe_le(le_id: str, entry: LeistungserbringerEntry) -> str:
        """Return a human-readable description with both id and alias."""
        if entry.alias == le_id:
            return f"'{le_id}'"
        return f"'{le_id}' (alias '{entry.alias}')"

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

    def resolve_inbox(self, submitter_id: str, inbox_name: str) -> InboxTarget:
        """Retrieve a specific inbox target by exact submitter (LE) ID or alias and inbox name.

        No auto-guessing or fallback search.
        """
        entry = self._le_by_id.get(submitter_id) or self._le_by_alias.get(submitter_id)
        if entry is None:
            available = ", ".join(self._describe_le(le_id, e) for le_id, e in self.leistungserbringer.items())
            log.error(f"Submitter '{submitter_id}' not found. Available: {available}")
            sys.exit(1)

        le_id = next(k for k, v in self.leistungserbringer.items() if v is entry)
        if inbox_name not in entry.inbox_buckets:
            available = ", ".join(entry.inbox_buckets.keys())
            log.error(
                f"Inbox '{inbox_name}' not configured for submitter {self._describe_le(le_id, entry)}. "
                f"Available: {available}"
            )
            sys.exit(1)

        inbox_cfg = entry.inbox_buckets[inbox_name]
        bucket = inbox_cfg.bucket or inbox_name
        return InboxTarget(
            s3=S3Options(bucket=bucket, **inbox_cfg.model_dump(exclude={"bucket"})),
            **inbox_cfg.model_dump(include={"private_key_path", "private_key_passphrase"}),
        )
