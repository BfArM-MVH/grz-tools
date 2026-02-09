"""Processing pipeline that orchestrates the streaming workflow."""

from __future__ import annotations

import contextlib
import logging
import queue
import threading
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

from grz_common.constants import TQDM_DEFAULTS
from tqdm.auto import tqdm

from .base import PipelineContext, StreamStage
from .compressors import GzipDecompressor
from .crypt4gh import Crypt4GHDecryptor, Crypt4GHEncryptor
from .s3 import S3Downloader, S3MultipartUploader
from .utils import abort_all_stages, drain_queue, finalize_stages_in_order, safe_join_thread
from .validators import BamValidator, FastqValidator, RawChecksumValidator

if TYPE_CHECKING:
    from grz_pydantic_models.submission.metadata.v1 import File as SubmissionFileMetadata

log = logging.getLogger(__name__)


@dataclass
class FileProcessingResult:
    """Result of processing a single file."""

    file_path: str
    success: bool
    errors: list[str] = field(default_factory=list)
    bytes_read: int = 0
    bytes_written: int = 0
    bytes_decrypted: int = 0
    checksum: str | None = None


class GrzctlProcessPipeline:
    """
    A specialized pipeline for the GRZ streaming workflow.

    Implements:
        S3 Download -> Decrypt -> TEE -> [Decompress -> Validate]
                                   |
                                   +---> Encrypt -> S3 Upload

    The decrypted stream is tee'd (duplicated):
    - One copy goes to validation (after decompression for format validators)
    - One copy goes to encryption and upload

    This allows validation and encryption to run on the same data without
    requiring disk storage. On validation failure, the upload is aborted.

    Key design principles:
    - Checksum validation happens on decrypted (but still compressed) data
    - Format validation (FASTQ/BAM) happens on decompressed data
    - Re-encryption happens on decrypted data (preserves compression)
    - Decompressor is a separate pluggable stage, not built into validators
    """

    def __init__(  # noqa: PLR0913
        self,
        source_client: Any,
        target_client: Any,
        source_bucket: str,
        target_bucket: str,
        decryption_key: bytes,
        encryption_public_key: bytes,
        signing_key: bytes | None = None,
        chunk_size: int = 65536,
        skip_validation: bool = False,
        max_concurrent_uploads: int = 1,
    ):
        """
        Initialize the pipeline.

        :param source_client: Boto3 S3 client for source bucket
        :param target_client: Boto3 S3 client for target bucket
        :param source_bucket: Source S3 bucket name
        :param target_bucket: Target S3 bucket name
        :param decryption_key: Private key for Crypt4GH decryption
        :param encryption_public_key: Public key for Crypt4GH encryption
        :param signing_key: Optional private key for signing
        :param chunk_size: Size of chunks for streaming
        :param skip_validation: Skip format validation
        :param max_concurrent_uploads: Maximum concurrent part uploads (1 = sequential)
        """
        self._source_client = source_client
        self._target_client = target_client
        self._source_bucket = source_bucket
        self._target_bucket = target_bucket
        self._decryption_key = decryption_key
        self._encryption_public_key = encryption_public_key
        self._signing_key = signing_key
        self._chunk_size = chunk_size
        self._skip_validation = skip_validation
        self._max_concurrent_uploads = max_concurrent_uploads
        self._log = log.getChild("DecryptValidateEncryptPipeline")

    def process_file(  # noqa: C901
        self,
        source_key: str,
        target_key: str,
        file_metadata: SubmissionFileMetadata,
        progress_callback: Callable[[int], None] | None = None,
    ) -> FileProcessingResult:
        """
        Process a single file through the complete pipeline.

        The pipeline tees the decrypted stream:
        - Validation branch: decompress (if needed) -> validate
        - Encryption branch: encrypt -> upload

        :param source_key: Source S3 object key
        :param target_key: Target S3 object key
        :param file_metadata: Metadata for the file
        :param progress_callback: Optional callback called with bytes processed
        :returns: Processing result with success/error info
        """
        file_path = str(file_metadata.file_path)
        context = PipelineContext()

        # Create pipeline stages
        downloader = S3Downloader(
            self._source_client,
            self._source_bucket,
            source_key,
        )
        decryptor = Crypt4GHDecryptor(self._decryption_key)
        encryptor = Crypt4GHEncryptor(
            self._encryption_public_key,
            self._signing_key,
        )
        uploader = S3MultipartUploader(
            self._target_client,
            self._target_bucket,
            target_key,
            max_concurrent_uploads=self._max_concurrent_uploads,
            # Use file size for dynamic part size calculation
            # (encrypted size is ~similar to plaintext + header overhead)
            expected_size=file_metadata.file_size_in_bytes,
        )

        # Create validation stages
        checksum_validator = RawChecksumValidator(
            expected_checksum=file_metadata.file_checksum,
            expected_size=file_metadata.file_size_in_bytes,
            name="DecryptedChecksumValidator",
        )

        # Format validator and decompressor (if needed)
        format_validator: FastqValidator | BamValidator | None = None
        decompressor: GzipDecompressor | None = None

        if not self._skip_validation:
            file_type = file_metadata.file_type
            is_gzipped = str(file_metadata.file_path).endswith(".gz")

            if file_type == "fastq":
                format_validator = FastqValidator(mean_read_length_threshold=0)
                if is_gzipped:
                    decompressor = GzipDecompressor()
            elif file_type == "bam":
                format_validator = BamValidator()
                # BAM files handle their own compression internally

        # Initialize all stages
        stages: list[StreamStage] = [downloader, decryptor, encryptor, uploader, checksum_validator]
        if format_validator:
            stages.append(format_validator)
        if decompressor:
            stages.append(decompressor)

        for stage in stages:
            stage.initialize(context)

        try:
            # Run the tee'd pipeline
            self._run_teed_pipeline(
                downloader=downloader,
                decryptor=decryptor,
                encryptor=encryptor,
                uploader=uploader,
                checksum_validator=checksum_validator,
                format_validator=format_validator,
                decompressor=decompressor,
                file_metadata=file_metadata,
                progress_callback=progress_callback,
            )

            # Finalize validation stages FIRST to collect any errors
            # Order matters: decompressor -> format_validator -> checksum_validator
            validation_stages: list[StreamStage] = []
            if decompressor:
                validation_stages.append(decompressor)
            if format_validator:
                validation_stages.append(format_validator)
            validation_stages.append(checksum_validator)

            validation_success = finalize_stages_in_order(validation_stages, context, self._log)

            # Check for validation errors BEFORE completing the upload
            if not validation_success or context.has_errors():
                self._log.error(f"Validation failed for {file_path}: {context.errors}")
                # Abort remaining stages
                abort_all_stages([uploader, decryptor, encryptor, downloader], self._log)
                return FileProcessingResult(
                    file_path=file_path,
                    success=False,
                    errors=context.errors,
                    bytes_read=context.bytes_read,
                    bytes_decrypted=context.bytes_decrypted,
                )

            # Validation passed - finalize remaining stages
            for stage in [decryptor, encryptor, uploader, downloader]:
                stage.finalize()

            return FileProcessingResult(
                file_path=file_path,
                success=True,
                errors=[],
                bytes_read=context.bytes_read,
                bytes_written=context.bytes_written,
                bytes_decrypted=context.bytes_decrypted,
                checksum=context.checksums.get("sha256"),
            )

        except Exception as e:
            self._log.error(f"Pipeline failed for {file_path}: {e}")
            # Abort all stages
            abort_all_stages(stages, self._log)
            return FileProcessingResult(
                file_path=file_path,
                success=False,
                errors=[str(e)],
            )

    def _run_teed_pipeline(  # noqa: C901, PLR0913, PLR0915, PLR0912
        self,
        downloader: S3Downloader,
        decryptor: Crypt4GHDecryptor,
        encryptor: Crypt4GHEncryptor,
        uploader: S3MultipartUploader,
        checksum_validator: RawChecksumValidator,
        format_validator: FastqValidator | BamValidator | None,
        decompressor: GzipDecompressor | None,
        file_metadata: SubmissionFileMetadata,
        progress_callback: Callable[[int], None] | None = None,
    ) -> None:
        """
        Execute the tee'd pipeline with parallel validation and upload.

        Data flow (parallel branches after decrypt):
            Download -> Decrypt -> TEE -> [Checksum + Decompress -> Validate] (Thread 1)
                                    |
                                    +---> [Encrypt -> Upload] (Thread 2)

        Validation and encryption/upload run in parallel threads, allowing the
        CPU-bound validation (especially decompression) to proceed while the
        network-bound upload happens concurrently.

        On any error (validation, encryption, upload), the pipeline aborts immediately
        and signals all threads to stop.

        :param progress_callback: Optional callback called with bytes processed per chunk
        """
        file_name = Path(file_metadata.file_path).name
        content_length = downloader.content_length or 0

        # sentinel to signal end of stream
        sentinel = object()

        # queues for parallel processing
        # using bounded queues to avoid memory issues with large files
        validation_queue: queue.Queue[bytes | object] = queue.Queue(maxsize=32)
        upload_queue: queue.Queue[bytes | object] = queue.Queue(maxsize=32)

        # error tracking and abort signaling for threads
        validation_error: list[Exception] = []
        upload_error: list[Exception] = []
        abort_event = threading.Event()

        def should_abort() -> bool:
            """Check if pipeline should abort due to internal error."""
            return abort_event.is_set()

        def validation_worker() -> None:  # noqa: C901
            """Worker thread for validation branch."""
            try:
                while not should_abort():
                    try:
                        chunk = validation_queue.get(timeout=0.1)
                    except queue.Empty:
                        continue

                    if chunk is sentinel:
                        break
                    if not isinstance(chunk, bytes) or not chunk:
                        continue

                    # checksum validation (on compressed data)
                    checksum_validator.observe(chunk)

                    # format validation (on decompressed data if needed)
                    if format_validator is not None:
                        if decompressor is not None:
                            decompressed = decompressor.process(chunk)
                            if decompressed:
                                format_validator.observe(decompressed)
                        else:
                            format_validator.observe(chunk)

                # flush decompressor (unless aborted)
                if not should_abort() and decompressor is not None and format_validator is not None:
                    final_decompressed = decompressor.flush()
                    if final_decompressed:
                        format_validator.observe(final_decompressed)

            except Exception as e:
                validation_error.append(e)
                abort_event.set()  # signal other threads to stop
                self._log.error(f"Validation worker error: {e}")

        def upload_worker() -> None:
            """Worker thread for encryption and upload branch."""
            try:
                while not should_abort():
                    try:
                        chunk = upload_queue.get(timeout=0.1)
                    except queue.Empty:
                        continue

                    if chunk is sentinel:
                        break
                    if not isinstance(chunk, bytes) or not chunk:
                        continue

                    # encrypt and upload
                    encrypted_out = encryptor.process(chunk)
                    if encrypted_out:
                        uploader.write(encrypted_out)

                # flush encryptor (unless aborted)
                if not should_abort():
                    final_encrypted = encryptor.flush()
                    if final_encrypted:
                        uploader.write(final_encrypted)

            except Exception as e:
                upload_error.append(e)
                abort_event.set()  # signal other threads to stop
                self._log.error(f"Upload worker error: {e}")

        # start worker threads
        validation_thread = threading.Thread(target=validation_worker, name="validation-worker")
        upload_thread = threading.Thread(target=upload_worker, name="upload-worker")
        validation_thread.start()
        upload_thread.start()

        try:
            with tqdm(
                total=content_length,
                desc="PROCESS ",
                postfix={"file": file_name},
                leave=False,
                **TQDM_DEFAULTS,
            ) as pbar:  # type: ignore[call-overload]
                # process chunks - download and decrypt in main thread, tee to workers
                for encrypted_chunk in downloader.iter_chunks(self._chunk_size):
                    # check for abort before processing (worker error or external interrupt)
                    if should_abort():
                        self._log.warning("Aborting pipeline due to error or interrupt")
                        break

                    chunk_size = len(encrypted_chunk)
                    decrypted = decryptor.process(encrypted_chunk)

                    if decrypted:
                        # tee to both queues (validation and upload run in parallel)
                        # use timeout to avoid blocking if abort happens
                        try:
                            validation_queue.put(decrypted, timeout=1.0)
                            upload_queue.put(decrypted, timeout=1.0)
                        except queue.Full:
                            if should_abort():
                                break
                            raise

                    pbar.update(chunk_size)
                    if progress_callback is not None:
                        progress_callback(chunk_size)

                # only flush if not aborted
                if not should_abort():
                    final_decrypted = decryptor.flush()
                    if final_decrypted:
                        try:
                            validation_queue.put(final_decrypted, timeout=1.0)
                            upload_queue.put(final_decrypted, timeout=1.0)
                        except queue.Full:
                            pass  # abort in progress

            # normal completion: send sentinel and wait for workers
            if not should_abort():
                validation_queue.put(sentinel)
                upload_queue.put(sentinel)
                safe_join_thread(validation_thread, timeout=30, logger=self._log)
                safe_join_thread(upload_thread, timeout=30, logger=self._log)
            else:
                # abort case: signal abort, drain queues, send sentinel
                abort_event.set()
                drain_queue(validation_queue)
                drain_queue(upload_queue)
                with contextlib.suppress(queue.Full):
                    validation_queue.put(sentinel)
                with contextlib.suppress(queue.Full):
                    upload_queue.put(sentinel)
                safe_join_thread(validation_thread, timeout=5, logger=self._log)
                safe_join_thread(upload_thread, timeout=5, logger=self._log)

        except Exception:
            # on exception, ensure workers are stopped
            abort_event.set()
            drain_queue(validation_queue)
            drain_queue(upload_queue)
            with contextlib.suppress(queue.Full):
                validation_queue.put(sentinel)
            with contextlib.suppress(queue.Full):
                upload_queue.put(sentinel)
            safe_join_thread(validation_thread, timeout=5, logger=self._log)
            safe_join_thread(upload_thread, timeout=5, logger=self._log)
            raise

        # raise any errors from workers
        if validation_error:
            raise validation_error[0]
        if upload_error:
            raise upload_error[0]
