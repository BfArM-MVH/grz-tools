from typing import Annotated, Any, Self

from grz_common.models.base import FilePath, IgnoringBaseSettings
from pydantic import Field, SecretStr, field_validator, model_validator

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

    private_key_passphrase: SecretStr | None = None
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

    known_public_keys: list[str] | None = None
    """Public keys that verify the signatures in the DB, one ``<format> <key> <author name>``
    entry per key, the same format as a line of ``known_public_keys_file``."""

    known_public_keys_file: FilePath | None = None
    """File with one ``<format> <key> <author name>`` line per key. Blank lines and lines
    starting with ``#`` are skipped."""

    @field_validator("known_public_keys", mode="before")
    @classmethod
    def reject_a_path_for_known_public_keys(cls, value: Any) -> Any:
        if isinstance(value, str):
            raise ValueError(
                "known_public_keys must be a list of keys, not a path. "
                "Put a path to a file with these lines into known_public_keys_file."
            )
        return value

    @model_validator(mode="after")
    def validate_known_public_keys(self) -> Self:
        if self.known_public_keys is not None and self.known_public_keys_file is not None:
            raise ValueError("Only one of known_public_keys or known_public_keys_file must be set.")
        return self
