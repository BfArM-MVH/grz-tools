"""Utilities for handling crypt4gh keys, encryption and decryption"""

import io
import logging
import os
import typing
from base64 import b64decode
from functools import partial
from getpass import getpass
from os import PathLike
from os.path import getsize
from pathlib import Path

import crypt4gh.header
import crypt4gh.keys
import crypt4gh.keys.c4gh
import crypt4gh.keys.ssh
import crypt4gh.lib
import grz_common.exceptions as grzexc
from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives.asymmetric.x25519 import X25519PrivateKey
from tqdm.auto import tqdm

from ..constants import TQDM_DEFAULTS
from .io import TqdmIOWrapper

log = logging.getLogger(__name__)


class Crypt4GH:
    """Crypt4GH encryption/decryption utility class"""

    Key = tuple[int, bytes, bytes]

    VERSION = 1
    SEGMENT_SIZE = 65536
    FILE_EXTENSION = ".c4gh"

    @staticmethod
    def prepare_c4gh_keys(
        recipient_key_file_path: str | PathLike,
        sender_private_key_file_path: str | PathLike | None = None,
        *,
        sender_private_key_bytes: bytes | None = None,
    ) -> tuple[Key]:
        """
        Prepare the key format that Crypt4GH needs. While it can contain multiple
         keys for multiple recipients, in our use case there is only a single recipient.

        If neither sender key is given, a random one is generated.

        :param recipient_key_file_path: path to the public key file of the recipient
        :param sender_private_key_file_path: path to the private key file of the sender.
        :param sender_private_key_bytes: the private key of the sender, as returned by
            :meth:`load_private_key`. Mutually exclusive with ``sender_private_key_file_path``.
        :raises ValueError: If both sender keys are given.
        """
        if sender_private_key_file_path is not None and sender_private_key_bytes is not None:
            raise ValueError("Only one of sender_private_key_file_path or sender_private_key_bytes must be given.")
        if sender_private_key_bytes is not None:
            sk = sender_private_key_bytes
        elif sender_private_key_file_path is not None:
            sk = Crypt4GH.retrieve_private_key(sender_private_key_file_path)
        else:
            sk = X25519PrivateKey.generate().private_bytes(
                encoding=serialization.Encoding.Raw,
                format=serialization.PrivateFormat.Raw,
                encryption_algorithm=serialization.NoEncryption(),
            )
        keys = ((0, sk, Crypt4GH.retrieve_public_key(recipient_key_file_path)),)
        return keys

    @staticmethod
    def encrypt_file(
        input_path: str | PathLike,
        output_path: str | PathLike,
        public_keys: tuple[Key],
    ):
        """
        Encrypt the file, properly handling the Crypt4GH header.

        :param public_keys:
        :param output_path:
        :param input_path:
        :return: tuple with md5 values for original file, encrypted file
        """
        # TODO: Progress bar?
        # TODO: store header in separate file?
        input_path = Path(input_path)
        output_path = Path(output_path)

        total_size = getsize(input_path)
        with (
            open(input_path, "rb") as in_fd,
            open(output_path, "wb") as out_fd,
            TqdmIOWrapper(
                typing.cast(io.RawIOBase, in_fd),
                tqdm(total=total_size, desc="ENCRYPT ", postfix=f"{input_path.name}", **TQDM_DEFAULTS),  # type: ignore[call-overload]
            ) as pbar_in_fd,
        ):
            crypt4gh.lib.encrypt(
                keys=public_keys,
                infile=pbar_in_fd,
                outfile=out_fd,
            )

    @staticmethod
    def retrieve_public_key(pubkey_path: str | PathLike) -> bytes:
        """
        Read Crypt4GH public key from specified path.

        :param pubkey_path: Path to the public key
        :returns: Public key bytes
        :raises ConfigurationError: If the key is missing or cannot be read.
        """
        try:
            return crypt4gh.keys.get_public_key(Path(pubkey_path).expanduser())
        except (OSError, ValueError, NotImplementedError) as e:
            # crypt4gh raises NotImplementedError for a file in no key format it knows
            raise grzexc.ConfigurationError(f"Public key {pubkey_path} cannot be read: {e}") from e

    @staticmethod
    def retrieve_private_key(seckey_path: str | PathLike, passphrase: str | None = None) -> bytes:
        """
        Read Crypt4GH private key from specified path.

        :param seckey_path: Path to the private key
        :param passphrase: Passphrase for the private key. If None, will check C4GH_PASSPHRASE envvar, if that is also undefined, will prompt for user input.
        :returns: Private key bytes
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
    def load_private_key(private_key: str | bytes, passphrase: str | None = None, key_name: str = "(inline)") -> bytes:
        """
        Load a Crypt4GH private key from its text, in memory.

        Supports the same formats as ``crypt4gh.keys.get_private_key``, which only reads from a file:
        a Crypt4GH private key and an OpenSSH private key, both in PEM format.
        The passphrase is only asked for if the key is encrypted. It is the first of: *passphrase*,
        the ``C4GH_PASSPHRASE`` environment variable, and an interactive prompt.

        :param private_key: The private key, as the content of a private key file.
        :param passphrase: Passphrase for the private key.
        :param key_name: Names the key in the passphrase prompt and in errors, which never show the key itself.
        :returns: Private key bytes
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
                return crypt4gh.keys.c4gh.parse_private_key(stream, passphrase_callback)

            magic_word += stream.read(len(crypt4gh.keys.ssh.MAGIC_WORD) - len(crypt4gh.keys.c4gh.MAGIC_WORD))
            if magic_word == crypt4gh.keys.ssh.MAGIC_WORD:
                return crypt4gh.keys.ssh.parse_private_key(stream, passphrase_callback)[0]
        except SystemExit as e:
            # crypt4gh exits the process for a key or a passphrase that it cannot use
            raise grzexc.ConfigurationError(f"Secret key {key_name} cannot be read with the given passphrase") from e

        raise grzexc.ConfigurationError(
            f"Secret key {key_name} cannot be read: it is neither a Crypt4GH nor an OpenSSH private key"
        )

    @staticmethod
    def key_opens_header(input_path: str | PathLike, private_key: bytes) -> bool:
        """
        Check whether a private key opens the Crypt4GH header of a file, without decrypting its body.

        :param input_path: Path to the encrypted file
        :param private_key: The private key, as returned by :meth:`load_private_key`
        :returns: ``True`` if the key decrypts at least one header packet.
        :raises DecryptionError: If the file has no valid Crypt4GH header.
        """
        with open(input_path, "rb") as in_fd:
            try:
                header_packets = list(crypt4gh.header.parse(in_fd))
            except ValueError as e:
                # crypt4gh raises ValueError for a header that the file gets wrong
                raise grzexc.DecryptionError(f"Cannot read the Crypt4GH header of {input_path}: {e}") from e

        # crypt4gh.header.decrypt logs every key that does not fit as an error, so try the key directly.
        # Like crypt4gh.header.decrypt_packet, this supports X25519 with ChaCha20-Poly1305 (method 0) only.
        for packet in header_packets:
            if int.from_bytes(packet[:4], byteorder="little") != 0:
                continue
            try:
                crypt4gh.header.decrypt_X25519_Chacha20_Poly1305(packet[4:], private_key)
            except Exception:  # noqa: S112
                continue
            return True
        return False

    @staticmethod
    def decrypt_file(input_path: Path, output_path: Path, private_key: bytes):
        """
        Decrypt a file using the provided private key
        :param input_path: Path to the encrypted file
        :param output_path: Path to the decrypted file
        :param private_key: The private key
        :raises DecryptionError: If the private key does not open the header,
            or if the header or a segment of the file cannot be decrypted.
        """
        total_size = getsize(input_path)
        file_name = input_path.name
        with (
            open(input_path, "rb") as in_fd,
            open(output_path, "wb") as out_fd,
            TqdmIOWrapper(
                typing.cast(io.RawIOBase, in_fd),
                tqdm(total=total_size, desc="DECRYPT ", postfix=f"{file_name}", **TQDM_DEFAULTS),  # type: ignore[call-overload]
            ) as pbar_in_fd,
        ):
            try:
                crypt4gh.lib.decrypt(
                    keys=[(0, private_key, None)],  # list of (method, privkey, recipient_pubkey=None),
                    infile=pbar_in_fd,
                    outfile=out_fd,
                )
            except ValueError as e:
                if str(e) == "No supported encryption method":
                    # crypt4gh raises this if the key opens no packet of the header
                    raise grzexc.DecryptionError(
                        f"Cannot decrypt {input_path}: the private key does not open its Crypt4GH header"
                    ) from e
                # crypt4gh raises ValueError for a header or a segment that the file gets wrong
                raise grzexc.DecryptionError(f"Cannot decrypt {input_path}: {e}") from e
