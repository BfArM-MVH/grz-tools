"""Hash calculation utilities."""

import hashlib
import logging
import traceback
from collections.abc import Callable
from os import PathLike
from os.path import getsize
from pathlib import Path
from typing import BinaryIO

from tqdm.auto import tqdm

from ..constants import TQDM_DEFAULTS

log = logging.getLogger(__name__)


class ChecksumTrackingReader:
    """
    A wrapper around a binary stream that calculates SHA256 checksum while reading.

    This allows calculating checksums on-the-fly during streaming operations
    without needing to read the data twice or store it on disk.

    Example usage:
        with open("file.txt", "rb") as f:
            reader = ChecksumTrackingReader(f)
            while chunk := reader.read(8192):
                process(chunk)
            checksum = reader.hexdigest()
    """

    def __init__(self, stream: BinaryIO, callback: Callable[[int], None] | None = None):
        """
        Initialize the checksum tracking reader.

        :param stream: The underlying binary stream to read from
        :param callback: Optional callback(bytes_read) for progress tracking
        """
        self._stream = stream
        self._hasher = hashlib.sha256()
        self._bytes_read = 0
        self._callback = callback

    def read(self, size: int = -1) -> bytes:
        """Read data from the stream and update the checksum."""
        data = self._stream.read(size)
        if data:
            self._hasher.update(data)
            self._bytes_read += len(data)
            if self._callback:
                self._callback(len(data))
        return data

    def hexdigest(self) -> str:
        """Return the hex digest of all data read so far."""
        return self._hasher.hexdigest()

    @property
    def bytes_read(self) -> int:
        """Return the total number of bytes read."""
        return self._bytes_read

    def close(self):
        """Close the underlying stream."""
        self._stream.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


class ChecksumTrackingWriter:
    """
    A writer that calculates SHA256 checksum while writing data.

    This allows calculating checksums on-the-fly during streaming operations
    while also writing to an underlying destination.

    Example usage:
        output = io.BytesIO()
        writer = ChecksumTrackingWriter(output)
        writer.write(data)
        checksum = writer.hexdigest()
    """

    def __init__(self, stream, callback: Callable[[int], None] | None = None):
        """
        Initialize the checksum tracking writer.

        :param stream: The underlying stream to write to
        :param callback: Optional callback(bytes_written) for progress tracking
        """
        self._stream = stream
        self._hasher = hashlib.sha256()
        self._bytes_written = 0
        self._callback = callback

    def write(self, data: bytes) -> int:
        """Write data to the stream and update the checksum."""
        self._stream.write(data)
        self._hasher.update(data)
        self._bytes_written += len(data)
        if self._callback:
            self._callback(len(data))
        return len(data)

    def hexdigest(self) -> str:
        """Return the hex digest of all data written so far."""
        return self._hasher.hexdigest()

    @property
    def bytes_written(self) -> int:
        """Return the total number of bytes written."""
        return self._bytes_written

    def close(self):
        """Close the underlying stream."""
        if hasattr(self._stream, "close"):
            self._stream.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


class TeeWriter:
    """
    A writer that writes data to multiple destinations simultaneously.

    This enables true streaming where data flows through multiple processing
    steps (e.g., validation + encryption) without buffering the entire content.

    Example usage:
        from grz_common.pipeline import FastqValidator, PipelineContext

        context = PipelineContext()
        validator = FastqValidator()
        validator.initialize(context)
        encryptor_input = io.BytesIO()
        tee = TeeWriter([validator, encryptor_input])

        # Write decrypted data - goes to both validator and encryptor
        tee.write(decrypted_chunk)
    """

    def __init__(self, writers: list):
        """
        Initialize the tee writer.

        :param writers: List of file-like objects with write() method
        """
        self._writers = writers
        self._bytes_written = 0

    def write(self, data: bytes) -> int:
        """Write data to all destinations."""
        for writer in self._writers:
            writer.write(data)
        self._bytes_written += len(data)
        return len(data)

    @property
    def bytes_written(self) -> int:
        """Return total bytes written."""
        return self._bytes_written

    def close(self):
        """Close all writers that support close()."""  # noqa: D402
        for writer in self._writers:
            if hasattr(writer, "close"):
                try:
                    writer.close()
                except Exception:
                    log.exception("Failed to close writer")
                    traceback.print_exc()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


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
            with tqdm(total=total_size, desc="SHA256  ", postfix=f"{file_path.name}", **TQDM_DEFAULTS) as pbar:  # type: ignore[call-overload]
                while chunk := f.read(chunk_size):
                    sha256_hash.update(chunk)
                    pbar.update(len(chunk))
        else:
            while chunk := f.read(chunk_size):
                sha256_hash.update(chunk)
    return sha256_hash.hexdigest()
