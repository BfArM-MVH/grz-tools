"""
Utility module for basic file operations.
Includes functions for calculating checksums, streaming json objects, file system checks and Crypt4GH encryption/decryption.
"""

from __future__ import annotations

import hashlib
import io
import json
import logging
import os
import typing
from functools import partial
from getpass import getpass
from os import PathLike
from os.path import getsize
from pathlib import Path
from typing import TextIO

import crypt4gh.header  # type: ignore[import-untyped]
import crypt4gh.keys  # type: ignore[import-untyped]
import crypt4gh.lib  # type: ignore[import-untyped]
from nacl.public import PrivateKey
from tqdm.auto import tqdm

log = logging.getLogger(__name__)


def calculate_sha256(file_path: str | PathLike, chunk_size=2**16, progress=True) -> str:
    """
    Calculate the sha256 value of a file in chunks

    :param file_path: path to the file
    :param chunk_size: Chunk size in bytes
    :param progress: Print progress
    :return: calculated sha256 value of file_path
    """
    file_path = Path(file_path)
    total_size = getsize(file_path)
    sha256_hash = hashlib.sha256()
    with open(file_path, "rb") as f:
        if progress and (total_size > chunk_size):
            with tqdm(
                total=total_size,
                unit="B",
                unit_scale=True,
                desc=f"Calculating SHA256 {file_path.name}",
            ) as pbar:
                while chunk := f.read(chunk_size):
                    sha256_hash.update(chunk)
                    pbar.update(len(chunk))
        else:
            while chunk := f.read(chunk_size):
                sha256_hash.update(chunk)
    return sha256_hash.hexdigest()


def calculate_md5(file_path, chunk_size=2**16, progress=True) -> str:
    """
    Calculate the md5 value of a file in chunks

    :param file_path: path to the file
    :param chunk_size: Chunk size in bytes
    :param progress: Print progress
    :return: calculated md5 value of file_path
    """
    total_size = getsize(file_path)
    md5_hash = hashlib.md5()  # noqa: S324
    with open(file_path, "rb") as f:
        if progress and (total_size > chunk_size):
            with tqdm(total=total_size, unit="B", unit_scale=True, desc="Calculating MD5") as pbar:
                while chunk := f.read(chunk_size):
                    md5_hash.update(chunk)
                    pbar.update(len(chunk))
        else:
            while chunk := f.read(chunk_size):
                md5_hash.update(chunk)
    return md5_hash.hexdigest()


def read_multiple_json(input: TextIO, buffer_size=65536, max_buffer_size=134217728):
    """
    Read multiple JSON objects from a text stream.
    :param input:
    :param buffer_size:
    :param max_buffer_size:
    """
    decoder = json.JSONDecoder()
    buffer = io.StringIO()  # Use StringIO as the buffer

    while chunk := input.read(buffer_size):
        if len(buffer.getvalue()) + len(chunk) > max_buffer_size:
            raise MemoryError("Reached maximum buffer size while reading input")

        # append chunk to buffer
        buffer.write(chunk)

        while True:
            try:
                data = buffer.getvalue().lstrip()
                # Attempt to decode a JSON object from the current buffer content
                obj, idx = decoder.raw_decode(data)
                yield obj  # Process the decoded object

                # Reset the buffer with the unprocessed content
                buffer = io.StringIO(data[idx:])
            except json.JSONDecodeError:
                # If a JSONDecodeError occurs, we need more data, so break out to read the next chunk
                break

    # If there is any remaining content after the loop, try to process it
    remaining_data = buffer.getvalue().strip()
    if remaining_data != "":
        raise ValueError("Remaining data is not empty. Is there invalid JSON?")


def is_relative_subdirectory(relative_path: str | PathLike, root_directory: str | PathLike) -> bool:
    """
    Check if the target path is a subdirectory of the root path
    using os.path.commonpath() without checking the file system.

    :param relative_path: The target path.
    :param root_directory: The root directory.
    :return: True if relative_path is a subdirectory of root_directory, otherwise False.
    """
    # Convert both paths to absolute paths without resolving symlinks
    root_directory = os.path.abspath(root_directory)
    relative_path = os.path.abspath(relative_path)

    common_path = os.path.commonpath([root_directory, relative_path])

    # Check if the common path is equal to the root path
    return common_path == root_directory


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
        public_keys: tuple[Crypt4GH.Key],
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
    def decrypt_file(input_path, output_path, private_key: bytes):
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


class TqdmIOWrapper(io.RawIOBase):
    """
    Wrapper to record reads and writes in a tqdm progress bar.

    Example:
        very_long_input = io.StringIO("0123456789abcdef" * 100000)
        total_size = len(very_long_input.getvalue())

        batch_size = 10 ** 4  # 10kb
        with open("/dev/null", "w") as fd:
            with TqdmIOWrapper(fd, tqdm(
                    total=total_size,
                    desc="Printing: ",
                    unit="B",
                    unit_scale=True,
                    # unit_divisor=1024,  # make use of standard units e.g. KB, MB, etc.
                    miniters=1,
            )) as pbar_fd:
                while chunk := very_long_input.read(batch_size):
                    pbar_fd.write(chunk)

                    time.sleep(0.1)
    """

    def __init__(self, io_buf: io.RawIOBase, progress_bar):
        """

        :param io_buf: the buffer to wrap
        :param progress_bar: tqdm progress bar
        """
        self.io_buf = io_buf
        self.callback = progress_bar.update

    def write(self, data):
        """Write data to the buffer and update the progress bar"""
        nbytes_written = self.io_buf.write(data)
        if nbytes_written:
            self.callback(nbytes_written)
        return nbytes_written

    def read(self, size=-1) -> bytes | None:
        """Read data from the buffer and update the progress bar"""
        data = self.io_buf.read(size)
        if data:
            self.callback(len(data))

        return data

    def readinto(self, buffer, /):
        """Read data into a buffer and update the progress bar"""
        nbytes_written = self.io_buf.readinto(buffer)
        if nbytes_written:
            self.callback(nbytes_written)

        return nbytes_written

    def flush(self):
        """Flush the buffer"""
        # Ensure all data is flushed to the underlying binary IO object
        self.io_buf.flush()

    def close(self):
        """Close the buffer"""
        # Close the underlying binary IO object
        self.io_buf.close()


# class HashLoggingWriter(io.RawIOBase):
#     def __init__(self, binary_io: io.RawIOBase, hash_algorithm: str = 'sha256'):
#         self._binary_io = binary_io  # The underlying binary I/O object
#         self.hash_algorithm = hash_algorithm  # Hashing algorithm (e.g., 'md5', 'sha256', etc.)
#
#         # Create the hash object based on the specified algorithm
#         if hash_algorithm not in hashlib.algorithms_available:
#             raise ValueError(f"Unsupported hash algorithm: {hash_algorithm}")
#         self.hash_obj = hashlib.new(hash_algorithm)
#
#     def write(self, b: bytes) -> int:
#         # Update the hash with the bytes being written
#         self.hash_obj.update(b)
#         # Write the bytes to the underlying binary IO object
#         return self._binary_io.write(b)
#
#     def read(self, size=-1) -> bytes:
#         raise NotImplementedError("This class is write-only")
#
#     def readinto(self, b: bytearray) -> int:
#         raise NotImplementedError("This class is write-only")
#
#     def flush(self):
#         # Ensure all data is flushed to the underlying binary IO object
#         self._binary_io.flush()
#
#     def close(self):
#         # Close the underlying binary IO object
#         self._binary_io.close()
#
#     def get_hash(self):
#         # Return the hex digest of the hash
#         return self.hash_obj
