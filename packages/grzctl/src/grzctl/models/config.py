import contextlib
import logging
import sys
import tempfile
from collections.abc import Iterator
from contextvars import ContextVar
from pathlib import Path
from typing import Annotated, Any

import grz_common.exceptions as grzexc
import yaml
from grz_common.models.base import Crypt4GHPublicKey, IgnoringBaseModel, IgnoringBaseSettings, get_secret_value
from grz_common.models.identifiers import IdentifiersModel
from grz_common.models.s3 import S3ConnectionBase, S3Options
from grz_common.utils.crypt import Crypt4GH
from pydantic import Field, PrivateAttr, SecretStr, model_validator
from pydantic.fields import FieldInfo
from pydantic_settings import PydanticBaseSettingsSource

log = logging.getLogger(__name__)

from .db import DbModel
from .pruefbericht import PruefberichtModel

_config_ctx: ContextVar[dict[str, Any] | None] = ContextVar("_config_ctx", default=None)


def _check_key_fields(name: str, key: object | None, key_path: str | None, *, required: bool) -> None:
    """Check that at most one of the fields ``<name>`` and ``<name>_path`` is set.

    :param name: Name of the field with the inline key.
    :param key: Value of that field.
    :param key_path: Value of the field ``<name>_path``.
    :param required: Whether one of the two fields must be set.
    :raises ValueError: If both are set, or if neither is set although one is required.
    """
    if required and key is None and key_path is None:
        raise ValueError(f"Either {name} or {name}_path must be set.")
    if key is not None and key_path is not None:
        raise ValueError(f"Only one of {name} or {name}_path must be set.")


def _load_private_key(
    location: str, private_key: SecretStr | None, private_key_path: str | None, passphrase: SecretStr | None
) -> bytes:
    """Load a crypt4gh private key in memory, from its inline text or from its file.

    The passphrase is only asked for if the key is encrypted. It is the first of: *passphrase*,
    the ``C4GH_PASSPHRASE`` environment variable, and an interactive prompt.

    :param location: Config location of the inline key. Names the key in the prompt and in errors.
    :param private_key: The private key, given inline.
    :param private_key_path: Path to the private key.
    :param passphrase: Passphrase of the private key.
    :returns: The private key.
    :raises ConfigurationError: If neither the key nor its path is set, or if the key cannot be loaded.
    """
    if private_key is not None:
        return Crypt4GH.load_private_key(
            private_key.get_secret_value(), passphrase=get_secret_value(passphrase), key_name=location
        )
    if private_key_path is not None:
        return Crypt4GH.retrieve_private_key(private_key_path, passphrase=get_secret_value(passphrase))
    raise grzexc.ConfigurationError(f"Neither {location} nor {location}_path is set.")


class InboxConfig(S3ConnectionBase):
    """
    Configuration for a specific inbox.
    Includes connection details and the private key needed to decrypt its contents.
    """

    bucket: Annotated[str | None, Field(default=None)] = None
    """S3 bucket name. Defaults to the inbox name key if not set."""

    private_key_path: Annotated[str, Field(min_length=1)]
    """Path to the GRZ private key used to decrypt files from this inbox."""

    private_key_passphrase: SecretStr | None = None
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

    private_key_passphrase: SecretStr | None = None
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
    """Path to the GRZ public key (optional; encryption targets are configured via archives instead)."""


class ArchiveTarget(IgnoringBaseModel):
    """Encapsulates everything needed to write to a specific archive."""

    s3: S3Options
    """S3 connection details and bucket for this archive."""

    public_key: Crypt4GHPublicKey | None = None
    """The crypt4gh public key for re-encryption of files destined for this archive."""

    public_key_path: Annotated[str | None, Field(default=None, min_length=1)] = None
    """Path to the crypt4gh public key for re-encryption of files destined for this archive."""

    @model_validator(mode="after")
    def validate_public_key(self) -> "ArchiveTarget":
        _check_key_fields("public_key", self.public_key, self.public_key_path, required=True)
        return self

    @contextlib.contextmanager
    def public_key_file(self) -> Iterator[str]:
        """Give a path to the crypt4gh public key, for callers that need a file.

        ``public_key_path`` is given as is. ``public_key`` is written to a temporary file first,
        which is deleted again once the caller is done with it.

        :yields: Path to a file with the crypt4gh public key.
        """
        if self.public_key_path is not None:
            yield self.public_key_path
            return
        if self.public_key is None:
            raise RuntimeError("Either public_key or public_key_path must be set.")
        with tempfile.NamedTemporaryFile("w") as public_key_file:
            public_key_file.write(self.public_key)
            public_key_file.flush()
            yield public_key_file.name


class ArchivesConfig(IgnoringBaseModel):
    """Configuration for consented and non-consented archives."""

    consented: ArchiveTarget
    """Target definition for consented submissions."""

    non_consented: ArchiveTarget
    """Target definition for non-consented submissions."""

    signing_key: SecretStr | None = None
    """The GRZ crypt4gh private key that signs the files re-encrypted for either archive."""

    signing_key_path: Annotated[str | None, Field(default=None, min_length=1)] = None
    """Path to the GRZ crypt4gh private key that signs the files re-encrypted for either archive."""

    signing_key_passphrase: SecretStr | None = None
    """Passphrase to the GRZ crypt4gh private key that signs the files re-encrypted for either archive."""

    @model_validator(mode="after")
    def check_buckets_are_unique(self) -> "ArchivesConfig":
        if self.consented.s3.bucket == self.non_consented.s3.bucket:
            raise ValueError("consented and non-consented buckets must be distinct.")
        return self

    @model_validator(mode="after")
    def validate_signing_key(self) -> "ArchivesConfig":
        _check_key_fields("signing_key", self.signing_key, self.signing_key_path, required=True)
        return self

    def load_signing_key(self) -> bytes:
        """Load the signing key in memory.

        :returns: The signing key.
        :raises ConfigurationError: If the key cannot be loaded.
        """
        return _load_private_key(
            "archives.signing_key", self.signing_key, self.signing_key_path, self.signing_key_passphrase
        )


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
        """Retrieve a specific inbox target by exact submitter (LE) ID and inbox name.

        No auto-guessing, fallback, or alias lookup.
        """
        entry = self._le_by_id.get(submitter_id)
        if entry is None:
            available = ", ".join(self._describe_le(le_id, e) for le_id, e in self.leistungserbringer.items())
            log.error(f"Submitter '{submitter_id}' not found. Available: {available}")
            sys.exit(1)

        if inbox_name not in entry.inbox_buckets:
            available = ", ".join(entry.inbox_buckets.keys())
            log.error(
                f"Inbox '{inbox_name}' not configured for submitter {self._describe_le(submitter_id, entry)}. "
                f"Available: {available}"
            )
            sys.exit(1)

        inbox_cfg = entry.inbox_buckets[inbox_name]
        bucket = inbox_cfg.bucket or inbox_name
        return InboxTarget(
            s3=S3Options(bucket=bucket, **inbox_cfg.model_dump(exclude={"bucket"})),
            **inbox_cfg.model_dump(include={"private_key_path", "private_key_passphrase"}),
        )
