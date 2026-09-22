"""Utilities for handling crypt4gh keys, encryption and decryption."""

import logging
import os
from contextlib import nullcontext
from functools import partial
from getpass import getpass
from os import PathLike
from pathlib import Path

import crypt4gh.keys
from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives.asymmetric.x25519 import X25519PrivateKey
from grz_common.pipeline.components import ReadStream, Tee, TqdmObserver
from tqdm.auto import tqdm

from ..constants import TQDM_DEFAULTS
from ..exceptions import ConfigurationError

log = logging.getLogger(__name__)


class Crypt4GH:
    """Crypt4GH encryption/decryption utility class using the streaming pipeline."""

    Key = tuple[int, bytes, bytes]

    VERSION = 1
    SEGMENT_SIZE = 65536
    FILE_EXTENSION = ".c4gh"

    @staticmethod
    def prepare_c4gh_keys(
        recipient_key_file_path: str | PathLike,
        sender_private_key: str | PathLike | None = None,
    ) -> tuple[Key]:
        """
        Prepare the key format that Crypt4GH needs.

        :param recipient_key_file_path: Path to the public key file of the recipient
        :param sender_private_key: Path to the private key file of the sender.
            If None, will be generated randomly.
        """
        if sender_private_key is not None:
            sk = Crypt4GH.retrieve_private_key(sender_private_key)
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
    def retrieve_public_key(pubkey_path: str | PathLike) -> bytes:
        """
        Read Crypt4GH public key from specified path.

        :param pubkey_path: Path to the public key
        :returns: Public key bytes
        :raises ConfigurationError: If the key is missing or cannot be read.
        """
        try:
            return crypt4gh.keys.get_public_key(os.path.expanduser(str(pubkey_path)))
        except (OSError, ValueError, NotImplementedError) as e:
            # crypt4gh raises NotImplementedError for a file in no key format it knows
            raise ConfigurationError(f"Public key {pubkey_path} cannot be read: {e}") from e

    @staticmethod
    def retrieve_private_key(seckey_path: str | PathLike, passphrase: str | None = None) -> bytes:
        """
        Read Crypt4GH private key from specified path.

        :param seckey_path: Path to the private key
        :param passphrase: Passphrase for the private key. If None, will check C4GH_PASSPHRASE envvar, if that is also undefined, will prompt for user input.
        :returns: Private key bytes
        :raises ConfigurationError: If the key is missing, or cannot be read with the passphrase.
        """
        seckeypath = os.path.expanduser(str(seckey_path))
        if not os.path.exists(seckeypath):
            raise ConfigurationError(f"Secret key not found: {seckey_path}")

        if passphrase:
            passphrase_callback = lambda: passphrase
        elif global_passphrase := os.getenv("C4GH_PASSPHRASE"):
            passphrase_callback = lambda: global_passphrase
        else:
            passphrase_callback = partial(getpass, prompt=f"Passphrase for {seckey_path}: ")

        try:
            return crypt4gh.keys.get_private_key(seckeypath, passphrase_callback)
        except SystemExit as e:
            # crypt4gh exits the process for a key or a passphrase that it cannot use
            raise ConfigurationError(f"Secret key {seckey_path} cannot be read with the given passphrase") from e
        except (OSError, ValueError, NotImplementedError) as e:
            raise ConfigurationError(f"Secret key {seckey_path} cannot be read: {e}") from e

    @staticmethod
    def decrypt_file(
        input_path: str | PathLike,
        output_path: str | PathLike,
        private_key: bytes,
        show_progress: bool = True,
    ):
        """
        Decrypt a file using the Crypt4GH streaming pipeline.

        :param input_path: Path to the encrypted file
        :param output_path: Path to the decrypted file
        :param private_key: The private key bytes for decryption
        :param show_progress: Whether to show progress bar
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
            pipeline = pipeline | Crypt4GHDecryptor(private_key=private_key)
            pipeline >> out_fd
