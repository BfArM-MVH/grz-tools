from os import PathLike
from pathlib import Path
from typing import Annotated, Any, Self

import yaml
from pydantic import (
    AfterValidator,
    BaseModel,
    ConfigDict,
    SecretStr,
    SerializationInfo,
    SerializerFunctionWrapHandler,
    field_serializer,
)
from pydantic.types import PathType
from pydantic_settings import BaseSettings, SettingsConfigDict

FilePath = Annotated[Path, AfterValidator(lambda v: v.expanduser()), PathType("file")]


def _check_crypt4gh_public_key(public_key: str) -> str:
    """Check that a crypt4gh public key given as text is framed like one.

    Only the markers are checked, not the key itself. The check fails only if both the
    ``BEGIN CRYPT4GH PUBLIC KEY`` and the ``END CRYPT4GH PUBLIC KEY`` marker are missing.

    :param public_key: The public key.
    :returns: *public_key*, unchanged.
    :raises ValueError: If both markers are missing.
    """
    if "BEGIN CRYPT4GH PUBLIC KEY" not in public_key and "END CRYPT4GH PUBLIC KEY" not in public_key:
        raise ValueError("Invalid public key format")
    return public_key


Crypt4GHPublicKey = Annotated[str, AfterValidator(_check_crypt4gh_public_key)]
"""A crypt4gh public key as text. Its validation fails only if both its BEGIN and its END marker are missing."""


def get_secret_value(value: SecretStr | None) -> str | None:
    """Extract the plain-text value from a ``SecretStr``.

    This is the canonical way to retrieve secret values from Pydantic models that
    use ``SecretStr`` fields.  Use it at every consumption site (boto3 clients,
    cryptographic key loaders, etc.) rather than accessing the raw value.
    """
    if value is None:
        return None
    return value.get_secret_value()


class _RevealableSecrets:
    """In JSON mode, writes ``SecretStr`` fields in plain text when the caller asks for it.

    By default pydantic writes them as ``"**********"``. With ``context={"reveal_secrets": True}``,
    ``model_dump(mode="json")`` and ``model_dump_json()`` return output that loads back into an equal model.
    """

    @field_serializer("*", mode="wrap", when_used="json")
    def _serialize_secret(self, value: Any, handler: SerializerFunctionWrapHandler, info: SerializationInfo) -> Any:
        if isinstance(value, SecretStr) and (info.context or {}).get("reveal_secrets"):
            return value.get_secret_value()
        return handler(value)


class IgnoringBaseModel(_RevealableSecrets, BaseModel):
    model_config = ConfigDict(
        extra="ignore",
        validate_assignment=True,
        use_enum_values=True,
    )

    def to_yaml(self, fd):
        """Writes the configuration as YAML, with secrets in plain text, so that ``from_path`` loads it back."""
        data = self.model_dump(
            mode="json", exclude_none=True, exclude_unset=True, exclude_defaults=True, context={"reveal_secrets": True}
        )
        yaml.dump(data, fd)

    @classmethod
    def from_path(cls, path: str | PathLike) -> Self:
        """Reads the configuration file and validates it against the schema."""
        with open(path, encoding="utf-8") as f:
            config = cls(**yaml.safe_load(f))

        return config


class IgnoringBaseSettings(_RevealableSecrets, BaseSettings):
    model_config = SettingsConfigDict(
        extra="ignore",
        validate_assignment=True,
        use_enum_values=True,
        env_nested_delimiter="__",
        env_prefix="grz_",
        # errors would show the raw input, and SecretStr does not mask passphrases and keys there
        hide_input_in_errors=True,
    )

    def to_yaml(self, fd):
        """Writes the configuration as YAML, with secrets in plain text, so that ``from_path`` loads it back."""
        data = self.model_dump(
            mode="json", exclude_none=True, exclude_unset=True, exclude_defaults=True, context={"reveal_secrets": True}
        )
        yaml.dump(data, fd)

    @classmethod
    def from_path(cls, path: str | PathLike) -> Self:
        """Reads the configuration file and validates it against the schema."""
        with open(path, encoding="utf-8") as f:
            config = cls(**yaml.safe_load(f))

        return config
