"""Utilities for handling crypt4gh keys, encryption and decryption."""

import io
import logging
import os
from base64 import b64decode
from contextlib import nullcontext
from functools import partial
from getpass import getpass
from os import PathLike
from pathlib import Path

import crypt4gh.keys
import crypt4gh.keys.c4gh
import crypt4gh.keys.ssh
import grz_common.exceptions as grzexc
from cryptography.hazmat.primitives.asymmetric.x25519 import X25519PrivateKey, X25519PublicKey
from grz_common.pipeline.components import ReadStream, Tee, TqdmObserver
from tqdm.auto import tqdm

from ..constants import TQDM_DEFAULTS

log = logging.getLogger(__name__)


class Crypt4GH:
    """Crypt4GH encryption/decryption utility class using the streaming pipeline."""

    Key = tuple[int, bytes, bytes]

    VERSION = 1
    SEGMENT_SIZE = 65536
    FILE_EXTENSION = ".c4gh"

    @staticmethod
    def prepare_c4gh_keys(
        recipient_public_key: X25519PublicKey,
        sender_private_key: X25519PrivateKey | None = None,
    ) -> tuple[Key]:
        """
        Prepare the key format that Crypt4GH needs. While it can contain multiple
         keys for multiple recipients, in our use case there is only a single recipient.

        :param recipient_public_key: the public key of the recipient
        :param sender_private_key: the private key of the sender. If ``None``, a random one is generated.
        """
        if sender_private_key is None:
            sender_private_key = X25519PrivateKey.generate()
        # crypt4gh works with the raw 32 bytes of each key
        keys = ((0, sender_private_key.private_bytes_raw(), recipient_public_key.public_bytes_raw()),)
        return keys

    @staticmethod
    def encrypt_file(
        input_path: str | PathLike,
        output_path: str | PathLike,
        public_keys: tuple[Key],
        show_progress: bool = True,
    ):
        """
        Encrypt a file using the Crypt4GH streaming pipeline.

        :param input_path: Path to the input file
        :param output_path: Path to the output encrypted file
        :param public_keys: Prepared Crypt4GH keys for encryption
        :param show_progress: Whether to show progress bar
        """
        from ..pipeline.components.crypt4gh import Crypt4GHEncryptor  # noqa: PLC0415

        input_path = Path(input_path)
        output_path = Path(output_path)

        # extract public key and signing key from prepared keys tuple
        _, signing_key, public_key = public_keys[0]

        with (
            open(input_path, "rb") as in_fd,
            open(output_path, "wb") as out_fd,
            (
                tqdm(  # type: ignore[call-overload]
                    total=input_path.stat().st_size,
                    desc="ENCRYPT ",
                    postfix={"file": input_path.name},
                    **TQDM_DEFAULTS,
                )
                if show_progress
                else nullcontext()
            ) as pbar,
        ):
            pipeline = ReadStream(in_fd)
            if show_progress and pbar:
                pipeline = pipeline | Tee(TqdmObserver(pbar))
            pipeline = pipeline | Crypt4GHEncryptor(recipient_pubkey=public_key, sender_privkey=signing_key)
            pipeline >> out_fd

    @staticmethod
    def retrieve_public_key(pubkey_path: str | PathLike) -> X25519PublicKey:
        """
        Read Crypt4GH public key from specified path.

        :param pubkey_path: Path to the public key
        :returns: The public key
        :raises ConfigurationError: If the key is missing or cannot be read.
        """
        try:
            public_key = Path(pubkey_path).expanduser().read_bytes()
        except OSError as e:
            raise grzexc.ConfigurationError(f"Public key {pubkey_path} cannot be read: {e}") from e
        return Crypt4GH.load_public_key(public_key, key_name=str(pubkey_path))

    @staticmethod
    def load_public_key(public_key: str | bytes, key_name: str = "(inline)") -> X25519PublicKey:
        """
        Load a Crypt4GH public key from its text, in memory.

        Supports the same formats as ``crypt4gh.keys.get_public_key``, which only reads from a file:
        a Crypt4GH public key in PEM format, and an OpenSSH ``ssh-ed25519`` public key line.

        :param public_key: The public key, as the content of a public key file.
        :param key_name: Names the key in errors.
        :returns: The public key
        :raises ConfigurationError: If the key is in no supported format, or is malformed.
        """
        if isinstance(public_key, str):
            public_key = public_key.encode("utf-8")

        # Reads the lines like crypt4gh.keys.get_public_key does from a file, but checks both PEM markers
        lines = [line.strip() for line in public_key.splitlines() if line.strip()]
        if (
            lines
            and lines[0].startswith(b"-----BEGIN CRYPT4GH PUBLIC KEY")
            and lines[-1].startswith(b"-----END CRYPT4GH PUBLIC KEY")
        ):
            try:
                return X25519PublicKey.from_public_bytes(b64decode(b"".join(lines[1:-1])))
            except ValueError as e:
                # b64decode and from_public_bytes raise ValueError for a payload that is no base64 or not 32 bytes long
                raise grzexc.ConfigurationError(f"Public key {key_name} cannot be read: {e}") from e

        if lines and lines[0].startswith(b"ssh-ed25519 "):
            try:
                return X25519PublicKey.from_public_bytes(crypt4gh.keys.ssh.get_public_key(lines[0]))
            except (AssertionError, RuntimeError, ValueError) as e:
                # crypt4gh asserts the key type, and raises RuntimeError for a key that is no ed25519 point
                raise grzexc.ConfigurationError(
                    f"Public key {key_name} cannot be read: it is no valid ssh-ed25519 key"
                ) from e

        raise grzexc.ConfigurationError(
            f"Public key {key_name} cannot be read: it is neither a Crypt4GH nor an OpenSSH ssh-ed25519 public key"
        )

    @staticmethod
    def retrieve_private_key(seckey_path: str | PathLike, passphrase: str | None = None) -> X25519PrivateKey:
        """
        Read Crypt4GH private key from specified path.

        :param seckey_path: Path to the private key
        :param passphrase: Passphrase for the private key. If None, will check C4GH_PASSPHRASE envvar, if that is also undefined, will prompt for user input.
        :returns: The private key
        :raises ConfigurationError: If the key is missing, or cannot be read with the passphrase.
        """
        seckeypath = Path(seckey_path).expanduser()
        if not seckeypath.exists():
            raise grzexc.ConfigurationError(f"Secret key not found: {seckey_path}")

        try:
            private_key = seckeypath.read_bytes()
        except OSError as e:
            raise grzexc.ConfigurationError(f"Secret key {seckey_path} cannot be read: {e}") from e
        return Crypt4GH.load_private_key(private_key, passphrase=passphrase, key_name=str(seckey_path))

    @staticmethod
    def load_private_key(
        private_key: str | bytes, passphrase: str | None = None, key_name: str = "(inline)"
    ) -> X25519PrivateKey:
        """
        Load a Crypt4GH private key from its text, in memory.

        Supports the same formats as ``crypt4gh.keys.get_private_key``, which only reads from a file:
        a Crypt4GH private key and an OpenSSH private key, both in PEM format.
        The passphrase is only asked for if the key is encrypted. It is the first of: *passphrase*,
        the ``C4GH_PASSPHRASE`` environment variable, and an interactive prompt.

        :param private_key: The private key, as the content of a private key file.
        :param passphrase: Passphrase for the private key.
        :param key_name: Names the key in the passphrase prompt and in errors, which never show the key itself.
        :returns: The private key
        :raises ConfigurationError: If the key is in no supported format, or cannot be read with the passphrase.
        """
        if isinstance(private_key, str):
            private_key = private_key.encode("utf-8")

        # Mirrors crypt4gh.keys.load_from_pem, which reads the same lines from a file
        lines = [line.strip() for line in private_key.splitlines() if line.strip()]
        if not lines or not lines[0].startswith(b"-----BEGIN ") or not lines[-1].startswith(b"-----END "):
            raise grzexc.ConfigurationError(f"Secret key {key_name} cannot be read: it is not in PEM format")
        try:
            stream = io.BytesIO(b64decode(b"".join(lines[1:-1])))
        except ValueError as e:
            raise grzexc.ConfigurationError(f"Secret key {key_name} cannot be read: {e}") from e

        if passphrase:
            passphrase_callback = lambda: passphrase
        elif global_passphrase := os.getenv("C4GH_PASSPHRASE"):
            passphrase_callback = lambda: global_passphrase
        else:
            passphrase_callback = partial(getpass, prompt=f"Passphrase for {key_name}: ")

        # Mirrors crypt4gh.keys.get_private_key
        magic_word = stream.read(len(crypt4gh.keys.c4gh.MAGIC_WORD))
        try:
            if magic_word == crypt4gh.keys.c4gh.MAGIC_WORD:
                return X25519PrivateKey.from_private_bytes(
                    crypt4gh.keys.c4gh.parse_private_key(stream, passphrase_callback)
                )

            magic_word += stream.read(len(crypt4gh.keys.ssh.MAGIC_WORD) - len(crypt4gh.keys.c4gh.MAGIC_WORD))
            if magic_word == crypt4gh.keys.ssh.MAGIC_WORD:
                return X25519PrivateKey.from_private_bytes(
                    crypt4gh.keys.ssh.parse_private_key(stream, passphrase_callback)[0]
                )
        except SystemExit as e:
            # crypt4gh exits the process for a key or a passphrase that it cannot use
            raise grzexc.ConfigurationError(f"Secret key {key_name} cannot be read with the given passphrase") from e

        raise grzexc.ConfigurationError(
            f"Secret key {key_name} cannot be read: it is neither a Crypt4GH nor an OpenSSH private key"
        )

    @staticmethod
    def decrypt_file(
        input_path: str | PathLike,
        output_path: str | PathLike,
        private_key: X25519PrivateKey,
        show_progress: bool = True,
    ):
        """
        Decrypt a file using the Crypt4GH streaming pipeline.

        :param input_path: Path to the encrypted file
        :param output_path: Path to the decrypted file
        :param private_key: The private key
        :param show_progress: Whether to show progress bar
        :raises DecryptionError: If the private key does not open the header,
            or if the header or a segment of the file cannot be decrypted.
        """
        from ..pipeline.components.crypt4gh import Crypt4GHDecryptor  # noqa: PLC0415

        input_path = Path(input_path)
        output_path = Path(output_path)

        with (
            open(input_path, "rb") as in_fd,
            open(output_path, "wb") as out_fd,
            (
                tqdm(  # type: ignore[call-overload]
                    total=input_path.stat().st_size,
                    desc="DECRYPT ",
                    postfix={"file": input_path.name},
                    **TQDM_DEFAULTS,
                )
                if show_progress
                else nullcontext()
            ) as pbar,
        ):
            pipeline = ReadStream(in_fd)
            if show_progress and pbar:
                pipeline = pipeline | Tee(TqdmObserver(pbar))
            # the decryptor works with the raw 32 bytes of the key
            pipeline = pipeline | Crypt4GHDecryptor(private_key=private_key.private_bytes_raw())
            try:
                pipeline >> out_fd
            except grzexc.DecryptionError as e:
                raise grzexc.DecryptionError(f"Cannot decrypt {input_path}: {e}") from e
