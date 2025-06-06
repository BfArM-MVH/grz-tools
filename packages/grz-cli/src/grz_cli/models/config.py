import logging
from os import PathLike
from typing import Self, override

from grz_common.models.keys import KeyConfigModel, KeyModel
from grz_common.models.s3 import S3ConfigModel

log = logging.getLogger(__name__)


class UploadConfig(S3ConfigModel):
    @classmethod
    @override
    def from_path(cls, path: str | PathLike) -> Self:
        """
        Reads the configuration file and validates it against the schema.
        Overrides IgnoringBaseSettings.from_path to handle deprecated schemas,
        where the "s3" key is named "s3_options".
        """
        import yaml

        with open(path, encoding="utf-8") as f:
            contents = yaml.safe_load(f)
            if not contents.get("s3"):
                s3_options = contents.get("s3_options") or {}
                log.warning(
                    "The key 's3_options' has been deprecated in favour of 's3'. Support will be removed in a future release. Please rename 's3_options' to 's3'."
                )
                contents["s3"] = s3_options
            config = cls(**contents)  # noqa:

        return config


class EncryptConfig(KeyConfigModel):
    @classmethod
    @override
    def from_path(cls, path: str | PathLike) -> Self:
        """
        Reads the configuration file and validates it against the schema.
        Overrides IgnoringBaseSettings.from_path to handle deprecated schemas,
        where public/private keys aren't grouped under a "keys" parent.
        """
        import yaml

        with open(path, encoding="utf-8") as f:
            contents = yaml.safe_load(f)
            if not contents.get("keys"):
                key_fields = KeyModel.model_fields.keys()
                log.warning(
                    f"Bare keys {', '.join(key_fields)} are deprecated. Support will be removed in a future release. Please group them under 'keys'."
                )
                keys = {field: contents.get(field) for field in key_fields}
                contents["keys"] = keys
            config = cls(**contents)  # noqa:

        return config
