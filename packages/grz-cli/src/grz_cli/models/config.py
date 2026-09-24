import logging
from typing import Self

from grz_common.models.base import Crypt4GHPublicKey, FilePath, IgnoringBaseModel
from grz_common.models.identifiers import IdentifiersConfigModel
from grz_common.models.s3 import S3ConfigModel
from pydantic import model_validator

log = logging.getLogger(__name__)


class KeyModel(IgnoringBaseModel):
    grz_public_key: Crypt4GHPublicKey | None = None
    """
    The public key of the recipient (the associated GRZ).
    """

    grz_public_key_path: FilePath | None = None
    """
    Path to the crypt4gh public key of the recipient (the associated GRZ).
    """

    submitter_private_key_path: FilePath | None = None
    """
    Path to the submitter's private key (optional).
    """

    @model_validator(mode="after")
    def validate_grz_public_key(self) -> Self:
        if self.grz_public_key is None and self.grz_public_key_path is None:
            raise ValueError("Either grz_public_key or grz_public_key_path must be set.")
        if self.grz_public_key is not None and self.grz_public_key_path is not None:
            raise ValueError("Only one of grz_public_key or grz_public_key_path must be set.")
        return self


class KeyConfigModel(IgnoringBaseModel):
    keys: KeyModel


class UploadConfig(S3ConfigModel):
    pass


class EncryptConfig(KeyConfigModel):
    pass


class ValidateConfig(IdentifiersConfigModel):
    pass
