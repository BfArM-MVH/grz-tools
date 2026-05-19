"""Utilities for handling crypt4gh keys, encryption and decryption"""

import io
import logging
import os
from functools import partial
from getpass import getpass
from os import PathLike
from pathlib import Path
from typing import Any, cast

import crypt4gh.header
import crypt4gh.keys
import crypt4gh.lib
from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives.asymmetric.x25519 import X25519PrivateKey
from grz_common.utils.io import TqdmIOWrapper

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
        sender_private_key: str | PathLike | None = None,
    ) -> tuple[Key]:
        """
        Prepare the key format that Crypt4GH needs. While it can contain multiple
         keys for multiple recipients, in our use case there is only a single recipient.

        :param recipient_key_file_path: path to the public key file of the recipient
        :param sender_private_key: path to the private key file of the sender.
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
        keys = ((0, sk, crypt4gh.keys.get_public_key(recipient_key_file_path)),)
        return keys

    @staticmethod
    def encrypt_file(
        input_path: str | PathLike,
        output_path: str | PathLike,
        public_keys: tuple[Key],
        progress_bar: Any | None = None,
    ):
        """
        Encrypt the file, properly handling the Crypt4GH header.

        :param public_keys:
        :param output_path:
        :param input_path:
        :param progress_bar: Anything that has an update(n) method, e.g., tqdm.tqdm
        """
        # TODO: store header in separate file?
        input_path = Path(input_path)
        output_path = Path(output_path)

        with (
            open(input_path, "rb") as in_fd,
            open(output_path, "wb") as out_fd,
        ):
            wrapped_in_fd = TqdmIOWrapper(cast(io.RawIOBase, in_fd), progress_bar) if progress_bar else in_fd

            crypt4gh.lib.encrypt(
                keys=public_keys,
                infile=wrapped_in_fd,
                outfile=out_fd,
            )

            if progress_bar:
                wrapped_in_fd.flush()

    @staticmethod
    def retrieve_private_key(seckey_path) -> bytes:
        """
        Read Crypt4GH private key from specified path.
        :param seckey_path: path to the private key
        :return:
        """
        seckeypath = os.path.expanduser(seckey_path)
        if not os.path.exists(seckeypath):
            raise ValueError("Secret key not found")

        passphrase = os.getenv("C4GH_PASSPHRASE")
        if passphrase:
            passphrase_callback = lambda: passphrase
        else:
            passphrase_callback = partial(getpass, prompt=f"Passphrase for {seckey_path}: ")

        return crypt4gh.keys.get_private_key(seckeypath, passphrase_callback)

    @staticmethod
    def decrypt_file(
        input_path: Path,
        output_path: Path,
        private_key: bytes,
        progress_bar: Any | None = None,
    ):
        """
        Decrypt a file using the provided private key
        :param input_path: Path to the encrypted file
        :param output_path: Path to the decrypted file
        :param private_key: The private key
        :param progress_bar: Anything that has an update(n) method, e.g., tqdm.tqdm
        """
        with (
            open(input_path, "rb") as in_fd,
            open(output_path, "wb") as out_fd,
        ):
            wrapped_in_fd = TqdmIOWrapper(cast(io.RawIOBase, in_fd), progress_bar) if progress_bar else in_fd
            crypt4gh.lib.decrypt(
                keys=[(0, private_key, None)],
                infile=wrapped_in_fd,
                outfile=out_fd,
            )

            if progress_bar:
                wrapped_in_fd.flush()
