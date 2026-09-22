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
