"""Utilities for handling crypt4gh keys, encryption and decryption."""

import logging
import os
import shutil
from functools import partial
from getpass import getpass
from os import PathLike
from pathlib import Path

import crypt4gh.keys
import tqdm.auto as tqdm
from grz_common.pipeline.components import TqdmObserver
from nacl.public import PrivateKey

log = logging.getLogger(__name__)


class Crypt4GH:
    """Crypt4GH encryption/decryption utility class."""

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
            sk = bytes(PrivateKey.generate())
        keys = ((0, sk, crypt4gh.keys.get_public_key(recipient_key_file_path)),)
        return keys

    @staticmethod
    def encrypt_file(
        input_path: str | PathLike,
        output_path: str | PathLike,
        public_keys: tuple[Key],
        show_progress: bool = True,
    ):
        """
        Encrypt a file using Crypt4GH.

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
                tqdm.tqdm(
                    total=input_path.stat().st_size,
                    desc="ENCRYPT",
                    postfix=input_path.name,
                    unit="B",
                    unit_scale=True,
                )
                if show_progress
                else None
            ) as pbar,
        ):
            if pbar:
                in_fd = TqdmObserver(in_fd, pbar=pbar)
            decrypted_stream = Crypt4GHEncryptor(in_fd, public_key, signing_key)
            shutil.copyfileobj(decrypted_stream, out_fd)

    @staticmethod
    def retrieve_private_key(seckey_path: str | PathLike) -> bytes:
        """
        Read Crypt4GH private key from specified path.

        :param seckey_path: Path to the private key
        :returns: Private key bytes
        """
        seckeypath = os.path.expanduser(str(seckey_path))
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
        input_path: str | Path,
        output_path: str | Path,
        private_key: bytes,
        show_progress: bool = True,
    ):
        """
        Decrypt a file using the provided private key.

        :param input_path: Path to the encrypted file
        :param output_path: Path to the decrypted file
        :param private_key: The private key bytes
        :param show_progress: Whether to show progress bar
        """
        from ..pipeline.components.crypt4gh import Crypt4GHDecryptor  # noqa: PLC0415

        input_path = Path(input_path)
        output_path = Path(output_path)
        
        with (
            open(input_path, "rb") as in_fd,
            open(output_path, "wb") as out_fd,
            (
                tqdm.tqdm(
                    total=input_path.stat().st_size,
                    desc="DECRYPT",
                    postfix=input_path.name,
                    unit="B",
                    unit_scale=True,
                )
                if show_progress
                else None
            ) as pbar,
        ):
            if pbar:
                in_fd = TqdmObserver(in_fd, pbar=pbar)

            decrypted_stream = Crypt4GHDecryptor(in_fd, private_key)
            shutil.copyfileobj(decrypted_stream, out_fd)
