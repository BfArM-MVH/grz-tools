from os import PathLike
from pathlib import Path
from typing import Annotated, Self

import yaml
from grz_common.utils.config import read_and_merge_config_files
from pydantic import AfterValidator, BaseModel, ConfigDict, SecretStr
from pydantic.types import PathType
from pydantic_settings import BaseSettings, SettingsConfigDict

FilePath = Annotated[Path, AfterValidator(lambda v: v.expanduser()), PathType("file")]


def get_secret_value(value: SecretStr | str | None) -> str | None:
    """Extract the plain-text value from a ``SecretStr`` (or pass through a plain ``str``).

    This is the canonical way to retrieve secret values from Pydantic models that
    use ``SecretStr`` fields.  Use it at every consumption site (boto3 clients,
    cryptographic key loaders, etc.) rather than accessing the raw value.
    """
    if value is None:
        return None
    if isinstance(value, SecretStr):
        return value.get_secret_value()
    return str(value)


def _mask_secrets(data: dict) -> dict:
    """Recursively replace ``SecretStr`` instances in a dict with ``"**********"``."""
    masked: dict = {}
    for key, value in data.items():
        if isinstance(value, SecretStr):
            masked[key] = "**********"
        elif isinstance(value, dict):
            masked[key] = _mask_secrets(value)
        elif isinstance(value, list):
            masked[key] = [
                _mask_secrets(v) if isinstance(v, dict) else ("**********" if isinstance(v, SecretStr) else v)
                for v in value
            ]
        else:
            masked[key] = value
    return masked


class IgnoringBaseModel(BaseModel):
    model_config = ConfigDict(
        extra="ignore",
        validate_assignment=True,
        use_enum_values=True,
    )

    def model_dump(
        self,
        *,
        mode: str = "python",
        exclude_none: bool = False,
        exclude_unset: bool = False,
        exclude_defaults: bool = False,
        **kwargs,
    ) -> dict:
        """Serialize the model.

        In ``mode="json"`` secret fields are masked with ``"**********"`` so
        that serialized output never leaks credentials.  Use ``mode="python"``
        (the default) when you need the actual ``SecretStr`` objects so that
        downstream code can call ``.get_secret_value()``.
        """
        data = super().model_dump(
            mode=mode,
            exclude_none=exclude_none,
            exclude_unset=exclude_unset,
            exclude_defaults=exclude_defaults,
            **kwargs,
        )
        if mode == "json":
            return _mask_secrets(data)
        return data

    def to_yaml(self, fd):
        """Reads the configuration file and validates it against the schema."""
        yaml.dump(self.model_dump(mode="json", exclude_none=True, exclude_unset=True, exclude_defaults=True), fd)

    @classmethod
    def from_path(cls, path: str | PathLike) -> Self:
        """Reads the configuration file and validates it against the schema."""
        with open(path, encoding="utf-8") as f:
            config = cls(**yaml.safe_load(f))

        return config


class IgnoringBaseSettings(BaseSettings):
    model_config = SettingsConfigDict(
        extra="ignore",
        validate_assignment=True,
        use_enum_values=True,
        env_nested_delimiter="__",
        env_prefix="grz_",
        env_file=".env",
    )

    def model_dump(
        self,
        *,
        mode: str = "python",
        exclude_none: bool = False,
        exclude_unset: bool = False,
        exclude_defaults: bool = False,
        **kwargs,
    ) -> dict:
        """Serialize the model.

        In ``mode="json"`` secret fields are masked with ``"**********"`` so
        that serialized output never leaks credentials.  Use ``mode="python"``
        (the default) when you need the actual ``SecretStr`` objects so that
        downstream code can call ``.get_secret_value()``.
        """
        data = super().model_dump(
            mode=mode,
            exclude_none=exclude_none,
            exclude_unset=exclude_unset,
            exclude_defaults=exclude_defaults,
            **kwargs,
        )
        if mode == "json":
            return _mask_secrets(data)
        return data

    def to_yaml(self, fd):
        """Reads the configuration file and validates it against the schema."""
        yaml.dump(self.model_dump(mode="json", exclude_none=True, exclude_unset=True, exclude_defaults=True), fd)

    @classmethod
    def from_path(cls, path: str | PathLike | list[str | PathLike]) -> Self:
        """Reads the configuration file and validates it against the schema."""
        if isinstance(path, tuple):
            path = list(path)

        if not isinstance(path, list):
            path = [path]

        paths = [Path(p) for p in path]
        config = read_and_merge_config_files(paths)
        return cls.model_validate(config)
