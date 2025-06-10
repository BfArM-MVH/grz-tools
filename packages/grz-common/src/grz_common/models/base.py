from os import PathLike
from pathlib import Path
from typing import Annotated, Self

from pydantic import AfterValidator, BaseModel, ConfigDict
from pydantic.types import PathType
from pydantic_settings import BaseSettings, SettingsConfigDict

FilePath = Annotated[Path, AfterValidator(lambda v: v.expanduser()), PathType("file")]


class IgnoringBaseModel(BaseModel):
    model_config = ConfigDict(
        extra="ignore",
        validate_assignment=True,
        use_enum_values=True,
    )

    def to_yaml(self, fd):
        """Reads the configuration file and validates it against the schema."""
        import yaml

        yaml.dump(self.model_dump(mode="json", exclude_none=True, exclude_unset=True, exclude_defaults=True))

    @classmethod
    def from_path(cls, path: str | PathLike) -> Self:
        """Reads the configuration file and validates it against the schema."""
        import yaml

        with open(path, encoding="utf-8") as f:
            config = cls(**yaml.safe_load(f))  # noqa:

        return config


class IgnoringBaseSettings(BaseSettings):
    model_config = SettingsConfigDict(
        extra="ignore",
        validate_assignment=True,
        use_enum_values=True,
        env_nested_delimiter="__",
        env_prefix="grz_",
    )

    def to_yaml(self, fd):
        """Reads the configuration file and validates it against the schema."""
        import yaml

        yaml.dump(self.model_dump(mode="json", exclude_none=True, exclude_unset=True, exclude_defaults=True), fd)

    @classmethod
    def from_path(cls, path: str | PathLike) -> Self:
        """Reads the configuration file and validates it against the schema."""
        import yaml

        with open(path, encoding="utf-8") as f:
            config = cls(**yaml.safe_load(f))  # noqa:

        return config
