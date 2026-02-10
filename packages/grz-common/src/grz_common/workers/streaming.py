"""Streaming pipeline worker for processing submissions."""

from __future__ import annotations

import json
import logging
import tempfile
import threading
from concurrent.futures import Future, ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import date
from os import PathLike
from pathlib import Path
from typing import TYPE_CHECKING, Any

from grz_pydantic_models.submission.metadata.v1 import File as SubmissionFileMetadata
from tqdm.auto import tqdm

from ..constants import TQDM_DEFAULTS, TqdmLoggingHandler
from ..models.s3 import S3Options
from ..pipeline import FileProcessingResult, GrzctlProcessPipeline, S3KeyBuilder, SignalManager
from ..progress import FileProgressLogger
from ..progress.states import ProcessingState
from ..transfer import init_s3_client, init_s3_resource
from ..utils.crypt import Crypt4GH
from ..utils.redaction import redact_file_patterns

if TYPE_CHECKING:
    from .submission import EncryptedSubmission

log = logging.getLogger(__name__)


class StreamingPipelineError(Exception):
    """Exception raised when a streaming pipeline operation fails."""

    pass


class PipelineInterruptedError(Exception):
    """Exception raised when the pipeline is interrupted by user (Ctrl+C)."""

    pass


@dataclass
class StreamingPipelineConfig:
    """Configuration for the streaming pipeline."""

    # source S3 (inbox)
    source_s3: S3Options

    # target S3 for consented submissions
    consented_archive_s3: S3Options
    # target S3 for non-consented submissions
    non_consented_archive_s3: S3Options

    # path to the GRZ private key for decryption
    grz_private_key_path: str | PathLike
    # path to the archive public key for re-encryption of consented submissions
    consented_archive_public_key_path: str | PathLike
    # path to the archive public key for re-encryption of non-consented submissions
    non_consented_archive_public_key_path: str | PathLike

    # number of parallel file processing threads
    threads: int = 1
    # whether to skip validation (use with caution)
    skip_validation: bool = False
    # chunk size for streaming operations
    chunk_size: int = 64 * 1024  # 64KB, matches Crypt4GH segment size
    # maximum concurrent part uploads per file (for multipart S3 uploads)
    max_concurrent_uploads: int = 4
    # patterns to redact from logs before uploading: list of (pattern, replacement) tuples
    redact_patterns: list[tuple[str, str]] | None = None


class StreamingPipelineWorker:
    """
    Worker that processes submissions in a streaming fashion.

    The workflow is: download -> decrypt -> [validate, re-encrypt -> archive]

    This minimizes disk I/O by processing data as streams wherever possible.
    The target archive is selected based on whether the submission is consented for research.

    Uses the modular pipeline components from grz_common.pipeline for all
    streaming operations, ensuring proper decompression for format validation.
    """

    __log = log.getChild("StreamingPipelineWorker")

    def __init__(
        self,
        config: StreamingPipelineConfig,
        status_file_path: str | PathLike,
        log_dir: str | PathLike,
    ):
        """
        Initialize the streaming pipeline worker.

        :param config: Pipeline configuration
        :param status_file_path: Path to the status/progress file
        :param log_dir: Path to the log directory
        """
        self._config = config
        self._status_file_path = Path(status_file_path)
        self._log_dir = Path(log_dir)

        # initialize S3 clients for source
        self._source_s3_client = init_s3_client(config.source_s3)

        # initialize S3 clients for both archive destinations
        self._consented_s3_client = init_s3_client(config.consented_archive_s3)
        self._consented_s3_resource = init_s3_resource(config.consented_archive_s3)
        self._non_consented_s3_client = init_s3_client(config.non_consented_archive_s3)
        self._non_consented_s3_resource = init_s3_resource(config.non_consented_archive_s3)

        # prepare decryption key
        self._private_key = Crypt4GH.retrieve_private_key(config.grz_private_key_path)

        # prepare encryption keys for both archives
        self._consented_keys = Crypt4GH.prepare_c4gh_keys(
            config.consented_archive_public_key_path, config.grz_private_key_path
        )
        self._non_consented_keys = Crypt4GH.prepare_c4gh_keys(
            config.non_consented_archive_public_key_path, config.grz_private_key_path
        )

        # track uploaded keys for cleanup on failure
        self._uploaded_keys: list[str] = []
        self._target_s3_client: Any = None  # will be set based on consent
        self._target_s3_options: S3Options | None = None  # will be set based on consent
        self._archive_keys: tuple[tuple[int, bytes, bytes]] | None = None  # will be set based on consent

        # interrupt handling
        self._interrupted = threading.Event()
        self._files_in_progress: set[str] = set()  # track files currently being processed
        self._files_in_progress_lock = threading.Lock()

        # tqdm-compatible logging handler
        self._tqdm_handler: TqdmLoggingHandler | None = None
        self._original_handlers: list[logging.Handler] = []

    def _setup_tqdm_logging(self) -> None:
        """Replace root logger handlers with tqdm-compatible handler during progress bar display."""
        root_logger = logging.getLogger()
        self._original_handlers = root_logger.handlers.copy()

        # Create tqdm handler with same formatter as existing handlers
        self._tqdm_handler = TqdmLoggingHandler()
        if self._original_handlers:
            self._tqdm_handler.setFormatter(self._original_handlers[0].formatter)

        # Replace handlers
        for handler in self._original_handlers:
            root_logger.removeHandler(handler)
        root_logger.addHandler(self._tqdm_handler)

    def _restore_logging(self) -> None:
        """Restore original logging handlers."""
        if self._tqdm_handler is None:
            return

        root_logger = logging.getLogger()
        root_logger.removeHandler(self._tqdm_handler)
        for handler in self._original_handlers:
            root_logger.addHandler(handler)
        self._tqdm_handler = None
        self._original_handlers = []

    def _check_interrupted(self) -> None:
        """Check if interrupted and raise if so."""
        if self._interrupted.is_set():
            raise PipelineInterruptedError("Pipeline interrupted by user")

    def _select_archive_for_submission(self, encrypted_submission: EncryptedSubmission) -> bool:
        """
        Determine consent status and select appropriate archive destination.

        :param encrypted_submission: The submission to check
        :returns: True if consented, False otherwise
        """
        is_consented = encrypted_submission.metadata.content.consents_to_research(date.today())
        self.__log.info(f"Submission consent status: {'consented' if is_consented else 'non-consented'}")

        if is_consented:
            self._target_s3_client = self._consented_s3_client
            self._target_s3_options = self._config.consented_archive_s3
            self._archive_keys = self._consented_keys
            self.__log.info("Using consented archive destination")
        else:
            self._target_s3_client = self._non_consented_s3_client
            self._target_s3_options = self._config.non_consented_archive_s3
            self._archive_keys = self._non_consented_keys
            self.__log.info("Using non-consented archive destination")

        return is_consented

    def _process_files_parallel(
        self,
        submission_id: str,
        files_to_process: list[tuple[Path, SubmissionFileMetadata]],
        progress_logger: FileProgressLogger[ProcessingState],
        global_pbar: tqdm,
    ) -> tuple[int, bool]:
        """
        Process files in parallel using ThreadPoolExecutor.

        :param submission_id: Submission ID
        :param files_to_process: List of (local_path, metadata) tuples
        :param progress_logger: Progress logger
        :param global_pbar: Global progress bar
        :returns: Tuple of (files_processed, interrupted)
        """
        files_processed = 0
        interrupted = False

        with ThreadPoolExecutor(max_workers=self._config.threads) as executor:
            futures: dict[Future[None], tuple[Path, SubmissionFileMetadata, str]] = {}

            # submit all files for processing
            for local_file_path, file_metadata in files_to_process:
                if self._interrupted.is_set():
                    self.__log.info("Interrupt received - not starting new files")
                    break

                source_key = self._get_source_key(submission_id, file_metadata)
                target_key = self._get_target_key(submission_id, file_metadata)

                future = executor.submit(
                    self._process_file_streaming,
                    local_file_path,
                    file_metadata,
                    source_key,
                    target_key,
                    progress_logger,
                    global_pbar,
                )
                futures[future] = (local_file_path, file_metadata, source_key)

            # wait for all submitted files to complete
            for future in as_completed(futures):
                local_file_path, file_metadata, source_key = futures[future]
                try:
                    future.result()
                    files_processed += 1
                except PipelineInterruptedError:
                    self.__log.info(f"File processing interrupted: {source_key}")
                    interrupted = True
                except Exception as e:
                    self.__log.error(f"File processing failed: {source_key}: {e}")
                    raise

            # check if we were interrupted
            if self._interrupted.is_set():
                interrupted = True

        return files_processed, interrupted

    def _process_files_sequential(
        self,
        submission_id: str,
        files_to_process: list[tuple[Path, SubmissionFileMetadata]],
        progress_logger: FileProgressLogger[ProcessingState],
        global_pbar: tqdm,
    ) -> tuple[int, bool]:
        """
        Process files sequentially (single-threaded).

        :param submission_id: Submission ID
        :param files_to_process: List of (local_path, metadata) tuples
        :param progress_logger: Progress logger
        :param global_pbar: Global progress bar
        :returns: Tuple of (files_processed, interrupted)
        """
        files_processed = 0
        interrupted = False

        for local_file_path, file_metadata in files_to_process:
            if self._interrupted.is_set():
                self.__log.info("Interrupt received - stopping file processing")
                interrupted = True
                break

            source_key = self._get_source_key(submission_id, file_metadata)
            target_key = self._get_target_key(submission_id, file_metadata)

            try:
                self._process_file_streaming(
                    local_file_path,
                    file_metadata,
                    source_key,
                    target_key,
                    progress_logger,
                    global_pbar,
                )
                files_processed += 1
            except PipelineInterruptedError:
                self.__log.info(f"File processing interrupted: {source_key}")
                interrupted = True
                break

        return files_processed, interrupted

    def _process_all_files(
        self,
        submission_id: str,
        encrypted_submission: EncryptedSubmission,
        progress_logger: FileProgressLogger[ProcessingState],
    ) -> int:
        """
        Process all files in the submission, either parallel or sequential.

        :param submission_id: Submission ID
        :param encrypted_submission: The submission to process
        :param progress_logger: Progress logger
        :returns: Number of files successfully processed
        :raises PipelineInterruptedError: If processing was interrupted
        """
        files_to_process = list(encrypted_submission.encrypted_files.items())
        total_files = len(files_to_process)
        total_bytes = sum(fm.file_size_in_bytes or 0 for _, fm in files_to_process)

        # create global progress bar
        with tqdm(
            total=total_bytes,
            desc="TOTAL   ",
            postfix=f"{total_files} files",
            position=0,
            leave=True,
            **TQDM_DEFAULTS,
        ) as global_pbar:  # type: ignore[call-overload]
            # process files based on thread configuration
            if self._config.threads > 1:
                files_processed, interrupted = self._process_files_parallel(
                    submission_id, files_to_process, progress_logger, global_pbar
                )
            else:
                files_processed, interrupted = self._process_files_sequential(
                    submission_id, files_to_process, progress_logger, global_pbar
                )

            # handle interruption
            if interrupted:
                self.__log.warning(
                    f"Pipeline interrupted! Processed {files_processed}/{total_files} files. "
                    "Completed files have been kept in the archive."
                )
                raise PipelineInterruptedError(
                    f"Pipeline interrupted by user. {files_processed}/{total_files} files completed."
                )

        return files_processed

    def process_submission(self, encrypted_submission: EncryptedSubmission) -> None:
        """
        Process an entire submission through the streaming pipeline.

        The target archive is determined by the submission's consent status.
        If any operation fails, all uploaded files are cleaned up from the archive.

        Handles Ctrl+C gracefully:
        - First Ctrl+C: Stop accepting new files, wait for in-progress files to complete
        - Second Ctrl+C: Force exit immediately

        :param encrypted_submission: The encrypted submission to process
        """
        submission_id = encrypted_submission.submission_id
        self.__log.info(f"Starting streaming pipeline for submission: {submission_id}")

        # initialize for processing
        self._interrupted.clear()
        self._files_in_progress.clear()
        self._uploaded_keys = []

        # determine consent status and select archive
        is_consented = self._select_archive_for_submission(encrypted_submission)

        # set up progress logging and signal handling
        progress_logger = FileProgressLogger[ProcessingState](self._status_file_path)
        self._setup_tqdm_logging()

        try:
            with SignalManager(self._interrupted, logger=self.__log):
                # process all files
                self._process_all_files(submission_id, encrypted_submission, progress_logger)

                # upload metadata and logs (metadata last as completion marker)
                self._upload_metadata_and_logs(encrypted_submission, is_consented)

                self.__log.info(f"Streaming pipeline completed for submission: {submission_id}")

        except PipelineInterruptedError:
            # re-raise interrupt errors (completed files are kept)
            raise
        except Exception as e:
            self.__log.error(f"Pipeline failed for submission {submission_id}: {e}")
            self.__log.info("Cleaning up incomplete submission from archive...")
            self._cleanup_failed_submission(submission_id)
            raise
        finally:
            self._restore_logging()

    def _cleanup_failed_submission(self, submission_id: str) -> None:
        """
        Remove all uploaded files from the archive after a failure.

        This ensures no incomplete or invalid submissions remain in the archive.

        :param submission_id: The submission ID for logging purposes
        """
        if not self._uploaded_keys:
            self.__log.debug("No files to clean up")
            return

        if self._target_s3_client is None or self._target_s3_options is None:
            self.__log.warning("Cannot cleanup: target S3 client not initialized")
            return

        self.__log.info(f"Removing {len(self._uploaded_keys)} uploaded files from archive")

        for key in self._uploaded_keys:
            try:
                self._target_s3_client.delete_object(
                    Bucket=self._target_s3_options.bucket,
                    Key=key,
                )
                self.__log.debug(f"Deleted: {key}")
            except Exception as delete_error:
                self.__log.warning(f"Failed to delete {key}: {delete_error}")

        self._uploaded_keys = []
        self.__log.info(f"Cleanup completed for submission: {submission_id}")

    def _get_source_key(self, submission_id: str, file_metadata: SubmissionFileMetadata) -> str:
        """Get the source S3 key for a file."""
        return S3KeyBuilder.source_key(submission_id, file_metadata)

    def _get_target_key(self, submission_id: str, file_metadata: SubmissionFileMetadata) -> str:
        """Get the target S3 key for a file."""
        return S3KeyBuilder.target_key(submission_id, file_metadata)

    def _process_file_streaming(  # noqa: PLR0913, C901
        self,
        local_file_path: Path,
        file_metadata: SubmissionFileMetadata,
        source_key: str,
        target_key: str,
        progress_logger: FileProgressLogger[ProcessingState],
        global_pbar: tqdm | None = None,
    ) -> None:
        """
        Process a single file through the streaming pipeline.

        Uses the modular DecryptValidateEncryptPipeline which properly
        handles decompression for format validation.

        :param local_file_path: Local file path reference
        :param file_metadata: Metadata for the file
        :param source_key: Source S3 object key
        :param target_key: Target S3 object key
        :param progress_logger: Progress logger
        :param global_pbar: Optional global progress bar to update
        :raises PipelineInterruptedError: If interrupted before processing starts
        """
        # check for interrupt before starting
        if self._interrupted.is_set():
            raise PipelineInterruptedError("Pipeline interrupted before file processing started")

        # check if already processed
        logged_state = progress_logger.get_state(local_file_path, file_metadata)
        if logged_state and logged_state.get("processing_successful"):
            self.__log.info(f"File '{source_key}' already processed, skipping.")
            # update global progress bar for skipped file
            if global_pbar is not None:
                global_pbar.update(file_metadata.file_size_in_bytes or 0)
            return

        # track this file as in-progress
        with self._files_in_progress_lock:
            self._files_in_progress.add(source_key)

        self.__log.info(f"Processing file: {source_key}")

        if self._target_s3_options is None:
            raise RuntimeError("Target S3 options not initialized")
        if self._target_s3_client is None:
            raise RuntimeError("Target S3 client not initialized")

        try:
            # get the encryption public key
            if not self._archive_keys:
                raise ValueError(
                    "No encryption public key provided, did you forget to run `self._select_archive_for_submission()`?"
                )
            else:
                encryption_public_key = self._archive_keys[0][2]

            # create the modular pipeline
            # note: files complete once started; interrupt is checked between files
            pipeline = GrzctlProcessPipeline(
                source_client=self._source_s3_client,
                target_client=self._target_s3_client,
                source_bucket=self._config.source_s3.bucket,
                target_bucket=self._target_s3_options.bucket,
                decryption_key=self._private_key,
                encryption_public_key=encryption_public_key,
                signing_key=self._private_key,
                chunk_size=self._config.chunk_size,
                skip_validation=self._config.skip_validation,
                max_concurrent_uploads=self._config.max_concurrent_uploads,
            )

            # create progress callback for global progress bar
            def progress_callback(bytes_processed: int) -> None:
                if global_pbar is not None:
                    global_pbar.update(bytes_processed)

            # process the file with progress callback
            result: FileProcessingResult = pipeline.process_file(
                source_key, target_key, file_metadata, progress_callback=progress_callback
            )

            if not result.success:
                error_msg = f"Processing failed for {source_key}: {result.errors}"
                self.__log.error(error_msg)
                progress_logger.set_state(
                    local_file_path,
                    file_metadata,
                    state=ProcessingState(
                        processing_successful=False,
                        validation_errors=result.errors,
                    ),
                )
                raise StreamingPipelineError(error_msg)

            # success
            self.__log.info(
                f"Successfully processed: {source_key} "
                f"({result.bytes_read}B -> {result.bytes_decrypted}B decrypted -> {result.bytes_written}B)"
            )
            progress_logger.set_state(
                local_file_path,
                file_metadata,
                state=ProcessingState(
                    processing_successful=True,
                    validation_errors=[],
                ),
            )

            # track uploaded file for potential cleanup
            self._uploaded_keys.append(target_key)

        except Exception as e:
            self.__log.error(f"Failed to process {source_key}: {e}")
            progress_logger.set_state(
                local_file_path,
                file_metadata,
                state=ProcessingState(
                    processing_successful=False,
                    validation_errors=[],
                    errors=[str(e)],
                ),
            )
            raise
        finally:
            # remove from in-progress tracking
            with self._files_in_progress_lock:
                self._files_in_progress.discard(source_key)

    def _upload_metadata_and_logs(self, encrypted_submission: EncryptedSubmission, is_consented: bool) -> None:
        """
        Upload log files and metadata to the archive.

        Metadata is uploaded LAST as a completion marker - this matches the
        non-streaming archive flow and allows detection of incomplete submissions.
        """
        # ensure target S3 options are set
        if self._target_s3_options is None:
            raise RuntimeError("Target S3 options not initialized")
        if self._target_s3_client is None:
            raise RuntimeError("Target S3 client not initialized")

        submission_id = encrypted_submission.submission_id
        target_metadata_key = S3KeyBuilder.metadata_key(submission_id)

        # upload log files FIRST (with redaction if patterns are configured)
        for log_path, log_key in encrypted_submission.get_log_files_and_object_id().items():
            if log_path.exists():
                # apply redaction patterns if configured
                if self._config.redact_patterns:
                    try:
                        self._redact_log_file(log_path)
                    except Exception as e:
                        self.__log.warning(f"Failed to redact log file {log_path}: {e}")

                self._target_s3_client.upload_file(
                    str(log_path),
                    self._target_s3_options.bucket,
                    log_key,
                )
                self._uploaded_keys.append(log_key)
                self.__log.info(f"Uploaded log: {log_key}")

        # create a redacted copy of metadata using the model's own method
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as temp_meta:
            temp_meta_path = temp_meta.name
            redacted_metadata = encrypted_submission.metadata.content.to_redacted_dict(is_consented)
            json.dump(redacted_metadata, temp_meta, indent=2)

        try:
            # upload metadata LAST as completion marker
            self._target_s3_client.upload_file(
                temp_meta_path,
                self._target_s3_options.bucket,
                target_metadata_key,
            )
            self._uploaded_keys.append(target_metadata_key)
            self.__log.info(f"Uploaded metadata to: {target_metadata_key} (completion marker)")
        finally:
            Path(temp_meta_path).unlink()

    def _redact_log_file(self, log_path: Path) -> None:
        """
        Redact sensitive information from a log file using configured patterns.

        :param log_path: Path to the log file to redact
        """
        if not self._config.redact_patterns:
            return

        modified = redact_file_patterns(log_path, self._config.redact_patterns)
        if modified:
            self.__log.debug(f"Redacted sensitive information from: {log_path.name}")
