"""Unified file operations for individual stage commands.

These operations provide file-based wrappers around the streaming pipeline
stages. They are used by the stage commands (download, decrypt, encrypt,
upload, validate) to provide consistent behavior with progress bars.

The streaming pipeline (GrzctlProcessPipeline) uses the low-level stages
directly for the TEE pattern required by parallel validation + encryption.

Each operation class:
- Uses streaming internally for memory efficiency
- Supports optional progress bar display
- Accepts an optional PipelineContext for state sharing
"""

from __future__ import annotations

import contextlib
import logging
from collections.abc import Iterator
from os import PathLike
from pathlib import Path
from typing import Any

import crypt4gh.keys

from ..utils.crypt import Crypt4GH
from .base import PipelineContext
from .compressors import GzipDecompressor
from .crypt4gh import Crypt4GHDecryptor, Crypt4GHEncryptor
from .s3 import S3Downloader, S3MultipartUploader
from .validators import BamValidator, FastqValidator, RawChecksumValidator

log = logging.getLogger(__name__)

# default chunk size for streaming operations (matches crypt4gh segment size)
STREAMING_CHUNK_SIZE = 64 * 1024


def _get_or_create_context(context: PipelineContext | None) -> PipelineContext:
    """Return provided context or create a new one."""
    return context if context is not None else PipelineContext()


@contextlib.contextmanager
def _progress_bar(
    total: int,
    desc: str,
    postfix: str,
    enabled: bool = True,
):
    """
    Context manager for optional tqdm progress bar.

    :param total: Total size for progress bar
    :param desc: Description label (8 chars recommended)
    :param postfix: Postfix text (usually filename)
    :param enabled: Whether to show progress bar
    :yields: Update callback function (no-op if disabled)
    """
    if not enabled or total <= 0:
        yield lambda n: None
        return

    from tqdm.auto import tqdm  # noqa: PLC0415

    from ..constants import TQDM_DEFAULTS  # noqa: PLC0415

    with tqdm(
        total=total,
        desc=desc,
        postfix=postfix,
        **TQDM_DEFAULTS,
    ) as pbar:  # type: ignore[call-overload]
        yield pbar.update


class DownloadOperation:
    """Download operation for S3 objects."""

    def __init__(
        self,
        s3_client: Any,
        bucket: str,
        chunk_size: int = STREAMING_CHUNK_SIZE,
    ):
        """
        Initialize download operation.

        :param s3_client: Boto3 S3 client
        :param bucket: S3 bucket name
        :param chunk_size: Chunk size for streaming
        """
        self._client = s3_client
        self._bucket = bucket
        self._chunk_size = chunk_size
        self._log = log.getChild("DownloadOperation")

    def _stream(self, key: str, context: PipelineContext) -> Iterator[bytes]:
        """Stream download from S3."""
        downloader = S3Downloader(self._client, self._bucket, key)
        downloader.initialize(context)

        try:
            yield from downloader.iter_chunks(self._chunk_size)
        finally:
            downloader.finalize()

    def to_file(
        self,
        key: str,
        output_path: Path,
        show_progress: bool = True,
        file_size: int | None = None,
        context: PipelineContext | None = None,
    ) -> int:
        """
        Download from S3 to a local file.

        :param key: S3 object key
        :param output_path: Local file path
        :param show_progress: Whether to show progress bar
        :param file_size: Optional file size for progress bar (fetched if not provided)
        :param context: Optional pipeline context for state sharing
        :returns: Number of bytes written
        """
        ctx = _get_or_create_context(context)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # get file size for progress bar if not provided
        if show_progress and file_size is None:
            try:
                head = self._client.head_object(Bucket=self._bucket, Key=key)
                file_size = head["ContentLength"]
            except Exception:
                file_size = 0

        bytes_written = 0
        with (
            output_path.open("wb") as f,
            _progress_bar(file_size or 0, "DOWNLOAD", output_path.name, show_progress) as update,
        ):
            for chunk in self._stream(key, ctx):
                f.write(chunk)
                bytes_written += len(chunk)
                update(len(chunk))

        return bytes_written


class DecryptOperation:
    """Crypt4GH decryption operation."""

    def __init__(
        self,
        private_key: bytes,
        chunk_size: int = STREAMING_CHUNK_SIZE,
    ):
        """
        Initialize decryption operation.

        :param private_key: Crypt4GH private key bytes
        :param chunk_size: Chunk size for streaming
        """
        self._private_key = private_key
        self._chunk_size = chunk_size
        self._log = log.getChild("DecryptOperation")

    @classmethod
    def from_key_file(
        cls,
        key_path: str | PathLike,
        chunk_size: int = STREAMING_CHUNK_SIZE,
    ) -> DecryptOperation:
        """
        Create DecryptOperation from a key file.

        :param key_path: Path to the private key file
        :param chunk_size: Chunk size for streaming
        """
        private_key = Crypt4GH.retrieve_private_key(key_path)
        return cls(private_key, chunk_size)

    def _stream(self, input_stream: Iterator[bytes], context: PipelineContext) -> Iterator[bytes]:
        """Decrypt a stream of encrypted data."""
        decryptor = Crypt4GHDecryptor(self._private_key)
        decryptor.initialize(context)

        try:
            for chunk in input_stream:
                decrypted = decryptor.process(chunk)
                if decrypted:
                    yield decrypted

            final = decryptor.flush()
            if final:
                yield final
        finally:
            decryptor.finalize()

    def _read_file_chunks(self, input_path: Path) -> Iterator[bytes]:
        """Read file in chunks."""
        with input_path.open("rb") as f:
            while chunk := f.read(self._chunk_size):
                yield chunk

    def file_to_file(
        self,
        input_path: Path,
        output_path: Path,
        show_progress: bool = True,
        context: PipelineContext | None = None,
    ) -> int:
        """
        Decrypt a file to another file.

        :param input_path: Path to encrypted file
        :param output_path: Path to output decrypted file
        :param show_progress: Whether to show progress bar
        :param context: Optional pipeline context for state sharing
        :returns: Number of bytes written
        """
        ctx = _get_or_create_context(context)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        file_size = input_path.stat().st_size

        bytes_written = 0
        bytes_read = 0
        with (
            output_path.open("wb") as out,
            _progress_bar(file_size, "DECRYPT ", input_path.name, show_progress) as update,
        ):
            for chunk in self._stream(self._read_file_chunks(input_path), ctx):
                out.write(chunk)
                bytes_written += len(chunk)
                # estimate read progress from write progress
                new_read = min(int(bytes_written / max(bytes_written, 1) * file_size), file_size)
                update(new_read - bytes_read)
                bytes_read = new_read
            update(file_size - bytes_read)

        return bytes_written


class EncryptOperation:
    """Crypt4GH encryption operation."""

    def __init__(
        self,
        recipient_public_key: bytes,
        signing_key: bytes | None = None,
        chunk_size: int = STREAMING_CHUNK_SIZE,
    ):
        """
        Initialize encryption operation.

        :param recipient_public_key: Recipient's Crypt4GH public key
        :param signing_key: Optional private key for signing
        :param chunk_size: Chunk size for streaming
        """
        self._public_key = recipient_public_key
        self._signing_key = signing_key
        self._chunk_size = chunk_size
        self._log = log.getChild("EncryptOperation")

    @classmethod
    def from_key_files(
        cls,
        public_key_path: str | PathLike,
        signing_key_path: str | PathLike | None = None,
        chunk_size: int = STREAMING_CHUNK_SIZE,
    ) -> EncryptOperation:
        """
        Create EncryptOperation from key files.

        :param public_key_path: Path to the recipient's public key
        :param signing_key_path: Optional path to the signing private key
        :param chunk_size: Chunk size for streaming
        """
        public_key = crypt4gh.keys.get_public_key(str(public_key_path))
        signing_key = None
        if signing_key_path:
            signing_key = Crypt4GH.retrieve_private_key(signing_key_path)
        return cls(public_key, signing_key, chunk_size)

    def _stream(self, input_stream: Iterator[bytes]) -> Iterator[bytes]:
        """Encrypt a stream of plaintext data."""
        encryptor = Crypt4GHEncryptor(self._public_key, self._signing_key)
        context = PipelineContext()
        encryptor.initialize(context)

        try:
            for chunk in input_stream:
                encrypted = encryptor.process(chunk)
                if encrypted:
                    yield encrypted

            final = encryptor.flush()
            if final:
                yield final
        finally:
            encryptor.finalize()

    def _read_file_chunks(self, input_path: Path) -> Iterator[bytes]:
        """Read file in chunks."""
        with input_path.open("rb") as f:
            while chunk := f.read(self._chunk_size):
                yield chunk

    def file_to_file(
        self,
        input_path: Path,
        output_path: Path,
        show_progress: bool = True,
    ) -> int:
        """
        Encrypt a file to another file.

        :param input_path: Path to plaintext file
        :param output_path: Path to output encrypted file
        :param show_progress: Whether to show progress bar
        :returns: Number of bytes written
        """
        output_path.parent.mkdir(parents=True, exist_ok=True)
        file_size = input_path.stat().st_size

        bytes_written = 0
        bytes_read = 0
        with (
            output_path.open("wb") as out,
            _progress_bar(file_size, "ENCRYPT ", input_path.name, show_progress) as update,
        ):
            for chunk in self._stream(self._read_file_chunks(input_path)):
                out.write(chunk)
                bytes_written += len(chunk)
                new_read = min(bytes_read + self._chunk_size, file_size)
                update(new_read - bytes_read)
                bytes_read = new_read
            update(file_size - bytes_read)

        return bytes_written


class UploadOperation:
    """Multipart upload operation for S3."""

    def __init__(
        self,
        s3_client: Any,
        bucket: str,
        max_concurrent_uploads: int = 4,
        chunk_size: int = STREAMING_CHUNK_SIZE,
    ):
        """
        Initialize upload operation.

        :param s3_client: Boto3 S3 client
        :param bucket: S3 bucket name
        :param max_concurrent_uploads: Maximum concurrent part uploads
        :param chunk_size: Chunk size for streaming
        """
        self._client = s3_client
        self._bucket = bucket
        self._max_concurrent_uploads = max_concurrent_uploads
        self._chunk_size = chunk_size
        self._log = log.getChild("UploadOperation")

    def _from_stream(
        self,
        key: str,
        input_stream: Iterator[bytes],
        expected_size: int | None = None,
    ) -> int:
        """Upload data from a stream to S3."""
        uploader = S3MultipartUploader(
            self._client,
            self._bucket,
            key,
            max_concurrent_uploads=self._max_concurrent_uploads,
            expected_size=expected_size,
        )
        context = PipelineContext()
        uploader.initialize(context)

        bytes_uploaded = 0
        try:
            for chunk in input_stream:
                uploader.write(chunk)
                bytes_uploaded += len(chunk)
            uploader.finalize()
        except Exception:
            uploader.abort()
            raise

        return bytes_uploaded

    def from_file(
        self,
        key: str,
        input_path: Path,
        show_progress: bool = True,
    ) -> int:
        """
        Upload a file to S3.

        :param key: S3 object key
        :param input_path: Path to the file to upload
        :param show_progress: Whether to show progress bar
        :returns: Number of bytes uploaded
        """
        file_size = input_path.stat().st_size

        uploader = S3MultipartUploader(
            self._client,
            self._bucket,
            key,
            max_concurrent_uploads=self._max_concurrent_uploads,
            expected_size=file_size,
        )
        context = PipelineContext()
        uploader.initialize(context)

        bytes_uploaded = 0
        try:
            with (
                input_path.open("rb") as f,
                _progress_bar(file_size, "UPLOAD  ", input_path.name, show_progress) as update,
            ):
                while chunk := f.read(self._chunk_size):
                    uploader.write(chunk)
                    bytes_uploaded += len(chunk)
                    update(len(chunk))
            uploader.finalize()
        except Exception:
            uploader.abort()
            raise

        return bytes_uploaded


class ValidateOperation:
    """Validation operations for submission files."""

    def __init__(
        self,
        chunk_size: int = STREAMING_CHUNK_SIZE,
    ):
        """
        Initialize validation operation.

        :param chunk_size: Chunk size for streaming
        """
        self._chunk_size = chunk_size
        self._log = log.getChild("ValidateOperation")

    def _read_file_chunks(self, file_path: Path) -> Iterator[bytes]:
        """Read file in chunks."""
        with file_path.open("rb") as f:
            while chunk := f.read(self._chunk_size):
                yield chunk

    def validate_checksum_stream(
        self,
        input_stream: Iterator[bytes],
        expected_checksum: str,
        expected_size: int,
    ) -> tuple[bool, list[str]]:
        """
        Validate checksum and size of a data stream.

        :param input_stream: Iterator of data chunks
        :param expected_checksum: Expected SHA256 checksum
        :param expected_size: Expected size in bytes
        :returns: Tuple of (passed, errors)
        """
        validator = RawChecksumValidator(
            expected_checksum=expected_checksum,
            expected_size=expected_size,
        )
        context = PipelineContext()
        validator.initialize(context)

        for chunk in input_stream:
            validator.observe(chunk)

        validator.finalize()
        return not context.has_errors(), context.errors

    def validate_checksum_file(
        self,
        file_path: Path,
        expected_checksum: str,
        expected_size: int,
    ) -> tuple[bool, list[str]]:
        """
        Validate checksum and size of a file.

        :param file_path: Path to the file
        :param expected_checksum: Expected SHA256 checksum
        :param expected_size: Expected size in bytes
        :returns: Tuple of (passed, errors)
        """
        return self.validate_checksum_stream(self._read_file_chunks(file_path), expected_checksum, expected_size)

    def validate_fastq_stream(
        self,
        input_stream: Iterator[bytes],
        is_gzipped: bool = False,
        mean_read_length_threshold: int = 0,
    ) -> tuple[bool, list[str]]:
        """
        Validate FASTQ format from a data stream.

        :param input_stream: Iterator of data chunks (compressed or not)
        :param is_gzipped: Whether the data is gzip compressed
        :param mean_read_length_threshold: Minimum mean read length
        :returns: Tuple of (passed, errors)
        """
        context = PipelineContext()
        decompressor = GzipDecompressor() if is_gzipped else None
        validator = FastqValidator(mean_read_length_threshold=mean_read_length_threshold)

        if decompressor:
            decompressor.initialize(context)
        validator.initialize(context)

        try:
            for chunk in input_stream:
                if decompressor:
                    decompressed = decompressor.process(chunk)
                    if decompressed:
                        validator.observe(decompressed)
                else:
                    validator.observe(chunk)

            if decompressor:
                final = decompressor.flush()
                if final:
                    validator.observe(final)
        finally:
            if decompressor:
                decompressor.finalize()
            validator.finalize()

        return not context.has_errors(), context.errors

    def validate_fastq_file(
        self,
        file_path: Path,
        mean_read_length_threshold: int = 0,
    ) -> tuple[bool, list[str]]:
        """
        Validate FASTQ format from a file.

        :param file_path: Path to the FASTQ file
        :param mean_read_length_threshold: Minimum mean read length
        :returns: Tuple of (passed, errors)
        """
        is_gzipped = str(file_path).endswith(".gz")
        return self.validate_fastq_stream(self._read_file_chunks(file_path), is_gzipped, mean_read_length_threshold)

    def validate_bam_file(self, file_path: Path) -> tuple[bool, list[str]]:
        """
        Validate BAM format from a file.

        Note: BAM validation requires a complete file due to index parsing.

        :param file_path: Path to the BAM file
        :returns: Tuple of (passed, errors)
        """
        context = PipelineContext()
        validator = BamValidator()
        validator.initialize(context)

        for chunk in self._read_file_chunks(file_path):
            validator.observe(chunk)

        validator.finalize()
        return not context.has_errors(), context.errors
