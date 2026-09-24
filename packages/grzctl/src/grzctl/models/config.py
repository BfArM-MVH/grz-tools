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
from grz_common.models.base import (
    Crypt4GHPublicKey,
    FilePath,
    IgnoringBaseModel,
    IgnoringBaseSettings,
    get_secret_value,
)
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


def _check_key_fields(name: str, key: object | None, key_path: Path | None, *, required: bool) -> None:
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
    location: str, private_key: SecretStr | None, private_key_path: Path | None, passphrase: SecretStr | None
) -> bytes:
    """Load a crypt4gh private key in memory, from its inline text or from its file.

    The passphrase is only asked for if the key is encrypted. It is the first of: *passphrase*,
    the ``C4GH_PASSPHRASE`` environment variable, and an interactive prompt.

    :param location: Config location of the inline key. Errors name the key by it, or by it with a ``_path`` suffix.
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
        try:
            return Crypt4GH.retrieve_private_key(private_key_path, passphrase=get_secret_value(passphrase))
        except grzexc.ConfigurationError as e:
            raise grzexc.ConfigurationError(f"{location}_path: {e}") from e
    raise grzexc.ConfigurationError(f"Neither {location} nor {location}_path is set.")


class InboxConfig(S3ConnectionBase):
    """
    Configuration for a specific inbox.
    Includes connection details and the private key needed to decrypt its contents.
    """

    bucket: Annotated[str | None, Field(default=None)] = None
    """S3 bucket name. Defaults to the inbox name key if not set."""

    private_key: SecretStr | None = None
    """The GRZ crypt4gh private key used to decrypt files from this inbox."""

    private_key_path: FilePath | None = None
    """Path to the GRZ crypt4gh private key used to decrypt files from this inbox."""

    private_key_passphrase: SecretStr | None = None
    """Passphrase to the GRZ private key used to decrypt files from this inbox."""

    @model_validator(mode="after")
    def validate_private_key(self) -> "InboxConfig":
        _check_key_fields("private_key", self.private_key, self.private_key_path, required=True)
        return self


class InboxTarget(IgnoringBaseModel):
    """
    Fully resolved source configuration.
    Encapsulates everything needed to read and decrypt from a specific inbox.
    """

    s3: S3Options
    """Fully resolved S3 options, including the bucket name."""

    private_key: SecretStr | None = None
    """The GRZ crypt4gh private key used to decrypt files from this inbox."""

    private_key_path: FilePath | None = None
    """Path to the GRZ crypt4gh private key used to decrypt files from this inbox."""

    private_key_passphrase: SecretStr | None = None
    """Passphrase to the GRZ private key used to decrypt files from this inbox."""

    @model_validator(mode="after")
    def validate_private_key(self) -> "InboxTarget":
        _check_key_fields("private_key", self.private_key, self.private_key_path, required=True)
        return self


class LeistungserbringerEntry(IgnoringBaseModel):
    """A single Leistungserbringer (submitter) with its inbox configurations."""

    alias: Annotated[str | None, Field(default=None)] = None
    """Human-friendly alias. Defaults to the LE id (dict key) if not set."""

    inbox_buckets: Annotated[dict[str, InboxConfig], Field(min_length=1)]
    """Mapping: InboxName -> InboxConfig."""


class ArchiveTarget(IgnoringBaseModel):
    """Encapsulates everything needed to write to a specific archive, and to decrypt submissions from it."""

    s3: S3Options
    """S3 connection details and bucket for this archive."""

    public_key: Crypt4GHPublicKey | None = None
    """The crypt4gh public key for re-encryption of files destined for this archive."""

    public_key_path: FilePath | None = None
    """Path to the crypt4gh public key for re-encryption of files destined for this archive."""

    private_key: SecretStr | None = None
    """The crypt4gh private key of this archive, to decrypt archived submissions (optional)."""

    private_key_path: FilePath | None = None
    """Path to the crypt4gh private key of this archive, to decrypt archived submissions (optional)."""

    private_key_passphrase: SecretStr | None = None
    """Passphrase to the crypt4gh private key of this archive."""

    @model_validator(mode="after")
    def validate_public_key(self) -> "ArchiveTarget":
        _check_key_fields("public_key", self.public_key, self.public_key_path, required=True)
        return self

    @model_validator(mode="after")
    def validate_private_key(self) -> "ArchiveTarget":
        _check_key_fields("private_key", self.private_key, self.private_key_path, required=False)
        return self

    @contextlib.contextmanager
    def public_key_file(self) -> Iterator[Path]:
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
            yield Path(public_key_file.name)


class ArchivesConfig(IgnoringBaseModel):
    """Configuration for consented and non-consented archives."""

    consented: ArchiveTarget
    """Target definition for consented submissions."""

    non_consented: ArchiveTarget
    """Target definition for non-consented submissions."""

    signing_key: SecretStr | None = None
    """The GRZ crypt4gh private key that signs the files re-encrypted for either archive."""

    signing_key_path: FilePath | None = None
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
            **inbox_cfg.model_dump(include={"private_key", "private_key_path", "private_key_passphrase"}),
        )

    def iter_decryption_keys(self, submitter_id: str) -> Iterator[tuple[str, bytes]]:
        """Load the private keys that may decrypt a submission of a submitter, one at a time.

        First come the keys of all inboxes of the submitter (LE), then the private keys of the consented
        and the non-consented archive, if set. A submitter missing from the config has no inbox keys.
        Each key is only loaded, and its passphrase only asked for, when the caller asks for it.

        Locations that share a key, for example through a YAML anchor, give one entry, so their key is
        loaded once. Two locations share a key if they name the same file or hold the same inline text.
        The entry takes its place from the first of them, and its passphrase from the first of them that sets one.

        :param submitter_id: Submitter (LE) ID, as in the submission's metadata.
        :yields: Pairs of the config locations of a key, joined by ``", "``, and the key.
        :raises ConfigurationError: If a key cannot be loaded.
        """
        candidates: list[tuple[str, InboxConfig | ArchiveTarget]] = []
        entry = self._le_by_id.get(submitter_id)
        if entry is None:
            log.warning(f"Submitter '{submitter_id}' is not in the config, so no inbox private key is tried.")
        else:
            for inbox_name, inbox in entry.inbox_buckets.items():
                candidates.append((f"leistungserbringer.{submitter_id}.inbox_buckets.{inbox_name}", inbox))
        for archive_name, target in (
            ("consented", self.archives.consented),
            ("non_consented", self.archives.non_consented),
        ):
            if target.private_key is not None or target.private_key_path is not None:
                candidates.append((f"archives.{archive_name}", target))

        # Group the locations by where their key comes from: the inline text, or the resolved path of the file.
        # dicts keep their insertion order, so the groups keep the order of their first location.
        groups: dict[tuple[str, str], list[tuple[str, InboxConfig | ArchiveTarget]]] = {}
        for prefix, holder in candidates:
            if holder.private_key is not None:
                source = ("inline", holder.private_key.get_secret_value())
            elif holder.private_key_path is not None:
                source = ("path", str(holder.private_key_path.resolve()))
            else:
                # validation leaves no inbox without a key, and the archives without one are not candidates
                continue
            groups.setdefault(source, []).append((prefix, holder))

        for group in groups.values():
            locations = ", ".join(
                f"{prefix}.private_key" if holder.private_key is not None else f"{prefix}.private_key_path"
                for prefix, holder in group
            )
            prefix, holder = group[0]
            passphrase = next(
                (h.private_key_passphrase for _, h in group if h.private_key_passphrase is not None), None
            )
            yield (
                locations,
                _load_private_key(f"{prefix}.private_key", holder.private_key, holder.private_key_path, passphrase),
            )
