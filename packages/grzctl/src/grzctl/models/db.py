from collections.abc import Iterable
from functools import cached_property
from pathlib import Path
from typing import Annotated, Any, Self

import platformdirs
from cryptography.exceptions import UnsupportedAlgorithm
from cryptography.hazmat.primitives.serialization import SSHPublicKeyTypes, load_ssh_public_key
from grz_common.models.base import FilePath, IgnoringBaseSettings, get_secret_value
from grz_db.errors import DatabaseConfigurationError
from grz_db.models.author import Author as SigningAuthor
from pydantic import Field, SecretStr, field_validator, model_validator

DEFAULT_KNOWN_PUBLIC_KEYS_FILE = Path(platformdirs.user_config_dir("grzctl")) / "known_public_keys"
"""The known public keys file that grzctl reads if the config sets neither ``known_public_keys`` nor
``known_public_keys_file``. It sits next to the default ``config.yaml``."""

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

    private_key: SecretStr | None = None
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


def _parse_known_public_keys(entries: Iterable[tuple[str, str]]) -> dict[str, list[SSHPublicKeyTypes]]:
    """Parse ``<format> <key> <comment>`` entries, one per key.

    The comment names the key's owner, and signature checks look keys up by it. So the comment
    is required, and it may contain spaces. Several keys may share a comment, for example across
    a key rotation, and all of them are kept.

    :param entries: Pairs of the entry's location for error messages, such as ``path:line``, and the entry.
    :returns: The public keys, grouped by their comment in the order of *entries*.
    :raises DatabaseConfigurationError: for an entry without a comment, or with a key that does not load.
    """
    public_keys: dict[str, list[SSHPublicKeyTypes]] = {}
    for location, entry in entries:
        parts = entry.split(maxsplit=2)
        if len(parts) < 3:
            raise DatabaseConfigurationError(
                f"{location}: expected '<format> <key> <comment>', where the comment names the key's owner."
            )
        try:
            public_keys.setdefault(parts[2], []).append(load_ssh_public_key(entry.encode()))
        except (ValueError, UnsupportedAlgorithm) as e:
            raise DatabaseConfigurationError(f"{location}: cannot load the public key: {e}") from e
    return public_keys


def _read_known_public_keys_file(path: Path) -> dict[str, list[SSHPublicKeyTypes]]:
    """Read a file with one ``<format> <key> <comment>`` line per key.

    Blank lines and lines starting with ``#`` are skipped.

    :param path: Path to the known public keys file.
    :returns: The public keys, grouped by their comment in file order.
    :raises DatabaseConfigurationError: if the file cannot be read, or for a line without a comment,
        or with a key that does not load.
    """
    try:
        with open(path) as f:
            lines = [(f"{path}:{line_number}", line.strip()) for line_number, line in enumerate(f, start=1)]
    except OSError as e:
        raise DatabaseConfigurationError(
            f"Cannot read the known public keys file {path}: {e}. "
            "Set db.known_public_keys or db.known_public_keys_file in the config."
        ) from e
    return _parse_known_public_keys((location, line) for location, line in lines if line and not line.startswith("#"))


class DbModel(IgnoringBaseSettings):
    """Submission database related configuration."""

    database_url: Annotated[str, Field(examples=["sqlite:///submission.sqlite"])]
    """URL to a database."""

    author: Author
    """Author information for submission database."""

    known_public_keys: list[str] | None = None
    """Public keys that verify the signatures in the DB, one ``<format> <key> <author name>``
    entry per key, the same format as a line of ``known_public_keys_file``. Every entry must be
    a key. Use YAML comments, not entries starting with ``#``."""

    known_public_keys_file: FilePath | None = None
    """File with one ``<format> <key> <author name>`` line per key. Blank lines and lines
    starting with ``#`` are skipped. If neither this nor ``known_public_keys`` is set, grzctl
    reads :data:`DEFAULT_KNOWN_PUBLIC_KEYS_FILE`."""

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

    @cached_property
    def public_keys_by_owner(self) -> dict[str, list[SSHPublicKeyTypes]]:
        """The known public keys, grouped by the owner that their comment names.

        Reads the file on first access, so only commands that verify signatures read it.

        :returns: The public keys, grouped by their comment in the order of the entries.
        :raises DatabaseConfigurationError: if the file cannot be read, or for an entry without a comment
            or with a key that does not load.
        """
        if self.known_public_keys is not None:
            return _parse_known_public_keys(
                (f"db.known_public_keys[{index}]", entry.strip()) for index, entry in enumerate(self.known_public_keys)
            )
        return _read_known_public_keys_file(self.known_public_keys_file or DEFAULT_KNOWN_PUBLIC_KEYS_FILE)

    @cached_property
    def signing_author(self) -> SigningAuthor:
        """The author signing this run's DB writes, holding the private key behind them.

        Cached on the configuration, so a command that opens several ``DbContext`` unlocks the
        key once and, without a configured passphrase, asks for it once.

        :returns: The author to hand to :class:`~grz_db.models.submission.SubmissionDb`.
        :raises ValueError: If neither ``private_key`` nor ``private_key_path`` is configured.
        """
        if self.author.private_key_path is not None:
            private_key_bytes = Path(self.author.private_key_path).read_bytes()
        elif self.author.private_key is not None:
            private_key_bytes = self.author.private_key.get_secret_value().encode("utf-8")
        else:
            raise ValueError("Either private_key or private_key_path must be provided.")

        return SigningAuthor(
            name=self.author.name,
            private_key_bytes=private_key_bytes,
            private_key_passphrase=get_secret_value(self.author.private_key_passphrase),
        )
