from typing import Annotated, Self

from grz_common.models.base import IgnoringBaseSettings
from pydantic import Field, FilePath, model_validator

# No whitespace (\s)
# No control characters (\x00-\x1f and \x7f)
AuthorNameStr = Annotated[
    str,
    Field(
        pattern=r"^[^\s\x00-\x1f\x7f]+$",
        min_length=1,
        description="A username without whitespace and control characters",
    ),
]


class Author(IgnoringBaseSettings):
    name: AuthorNameStr
    """Name of the author"""

    private_key: str | None = None
    """Author's private key (needed to sign DB modifications)."""

    private_key_path: FilePath | None = None
    """Path to the author's private key (needed to sign DB modifications)."""

    private_key_passphrase: str | None = None
    """Passphrase to author's private key (should almost always be provided in an environment variable)"""

    @model_validator(mode="after")
    def validate_private_key(self) -> Self:
        if self.private_key is not None and self.private_key_path is not None:
            raise ValueError("Only one of private_key or private_key_path must be set.")
        return self


class DbModel(IgnoringBaseSettings):
    """Submission database related configuration."""

    database_url: Annotated[str, Field(examples=["sqlite:///submission.sqlite"])]
    """URL to a database."""

    author: Author
    """Author information for submission database."""

    known_public_keys: FilePath | str = "~/.config/grz-cli/known_public_keys"
    """
    File listing public keys. Used for DB verification.

    Format: key_format public_key_base64 author_name"""
