from typing import Self

from grz_common.models.base import FilePath, IgnoringBaseSettings
from pydantic import field_validator, model_validator


class KeyModel(IgnoringBaseSettings):
    grz_public_key: str | None = None
    """
    The public key of the recipient (the associated GRZ).
    """

    grz_public_key_path: FilePath | None = None
    """
    Path to the crypt4gh public key of the recipient (the associated GRZ).
    """

    grz_private_key_path: FilePath | None = None
    """
    Path to the crypt4gh private key of the recipient (optional).
    """

    submitter_private_key_path: FilePath | None = None
    """
    Path to the submitter's private key (optional).
    """

    @field_validator("grz_public_key")
    @classmethod
    def check_grz_public_key(cls, v):
        if v is not None and "BEGIN CRYPT4GH PUBLIC KEY" not in v and "END CRYPT4GH PUBLIC KEY" not in v:
            raise ValueError("Invalid public key format")
        return v

    @model_validator(mode="after")
    def validate_grz_public_key(self) -> Self:
        if self.grz_public_key is None and self.grz_public_key_path is None:
            raise ValueError("Either grz_public_key or grz_public_key_path must be set.")
        if self.grz_public_key is not None and self.grz_public_key_path is not None:
            raise ValueError("Only one of grz_public_key or grz_public_key_path must be set.")
        return self


class KeyConfigModel(IgnoringBaseSettings):
    keys: KeyModel
