from collections.abc import Sequence
from os import PathLike
from pathlib import Path
from typing import Annotated, Any, Self

import yaml
from grz_common.utils.config import read_and_merge_config_files
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
        """Reads the configuration file and validates it against the schema."""
        yaml.dump(self.model_dump(mode="json", exclude_none=True, exclude_unset=True, exclude_defaults=True), fd)

    @classmethod
    def from_path(cls, path: str | PathLike | Sequence[str | PathLike]) -> Self:
        """Reads one or more configuration files, merges them in order, and validates the result against the schema."""
        paths = [path] if isinstance(path, str | PathLike) else list(path)
        return cls.model_validate(read_and_merge_config_files([Path(p) for p in paths]))


class IgnoringBaseSettings(_RevealableSecrets, BaseSettings):
    model_config = SettingsConfigDict(
        extra="ignore",
        validate_assignment=True,
        use_enum_values=True,
        env_nested_delimiter="__",
        env_prefix="grz_",
    )

    def to_yaml(self, fd):
        """Reads the configuration file and validates it against the schema."""
        yaml.dump(self.model_dump(mode="json", exclude_none=True, exclude_unset=True, exclude_defaults=True), fd)

    @classmethod
    def from_path(cls, path: str | PathLike | Sequence[str | PathLike]) -> Self:
        """Reads one or more configuration files, merges them in order, and validates the result against the schema."""
        paths = [path] if isinstance(path, str | PathLike) else list(path)
        return cls.model_validate(read_and_merge_config_files([Path(p) for p in paths]))
