"""Utilities for handling crypt4gh keys, encryption and decryption"""

import io
import logging
import os
import typing
from functools import partial
from getpass import getpass
from os import PathLike
from os.path import getsize
from pathlib import Path

import crypt4gh.header
import crypt4gh.keys
import crypt4gh.lib
from nacl.public import PrivateKey
from tqdm.auto import tqdm

from ..constants import TQDM_SMOOTHING
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
            sk = bytes(PrivateKey.generate())
        keys = ((0, sk, crypt4gh.keys.get_public_key(recipient_key_file_path)),)
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
                tqdm(
                    total=total_size,
                    desc=f"Encrypting: '{input_path.name}'",
                    unit="B",
                    unit_scale=True,
                    # unit_divisor=1024,  # make use of standard units e.g. KB, MB, etc.
                    miniters=1,
                    smoothing=TQDM_SMOOTHING,
                ),
            ) as pbar_in_fd,
        ):
            crypt4gh.lib.encrypt(
                keys=public_keys,
                infile=pbar_in_fd,
                outfile=out_fd,
            )

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
    def decrypt_file(input_path: Path, output_path: Path, private_key: bytes):
        """
        Decrypt a file using the provided private key
        :param input_path: Path to the encrypted file
        :param output_path: Path to the decrypted file
        :param private_key: The private key
        """
        total_size = getsize(input_path)
        with (
            open(input_path, "rb") as in_fd,
            open(output_path, "wb") as out_fd,
            TqdmIOWrapper(
                typing.cast(io.RawIOBase, in_fd),
                tqdm(
                    total=total_size,
                    desc=f"Decrypting: '{input_path.name}'",
                    unit="B",
                    unit_scale=True,
                    # unit_divisor=1024,  # make use of standard units e.g. KB, MB, etc.
                    miniters=1,
                ),
            ) as pbar_in_fd,
        ):
            crypt4gh.lib.decrypt(
                keys=[(0, private_key, None)],  # list of (method, privkey, recipient_pubkey=None),
                infile=pbar_in_fd,
                outfile=out_fd,
            )
