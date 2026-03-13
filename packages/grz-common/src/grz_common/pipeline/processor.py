import logging
from concurrent.futures import Future, ThreadPoolExecutor
from contextlib import ExitStack
from datetime import date
from io import BytesIO
from pathlib import Path
from typing import TYPE_CHECKING, Any

from crypt4gh.keys import get_public_key
from moto.kinesis.models import Stream

from grz_common.constants import TQDM_DEFAULTS
from grz_common.workers.submission import SubmissionMetadata
from grz_db.models.submission import SubmissionDb, SubmissionStateEnum
from grz_pydantic_models.submission.metadata import File, FileType
from grz_pydantic_models.submission.thresholds import Thresholds
from grzctl.commands.clean import _clean_submission_from_bucket
from grzctl.dbcontext import DbContext
from grzctl.models.config import CleanConfig, InboxTarget, ProcessConfig
from pydantic import AnyHttpUrl
from tqdm.auto import tqdm

from ..models.s3 import S3Options
from ..progress import FileProgressLogger, ProcessingState
from ..transfer import init_s3_client
from ..utils.crypt import Crypt4GH
from .components import (
    ObserverWithMetrics,
    Tee,
    TqdmObserver,
)
from .components.crypt4gh import Crypt4GHDecryptor, Crypt4GHEncryptor
from .components.perf import StreamMetricsRegistry
from .components.s3 import S3Downloader, S3MultipartUploader, calculate_s3_part_size
from .components.validation import BamValidator, ChecksumValidator, FastqValidator
from .context import ConsistencyValidator, SubmissionContext

log = logging.getLogger(__name__)

if TYPE_CHECKING:
    from types_boto3_s3.client import S3Client


class SubmissionProcessor:
    """
    Orchestrates the streaming submission pipeline:
    Inbox -> Decrypt -> Validate -> Re-Encrypt -> Archive
    """

    def __init__(  # noqa: PLR0913
        self,
        configuration: ProcessConfig,
        inbox: InboxTarget,
        status_file_path: Path,
        redact_patterns: list[tuple[str, str]] | None = None,
        clean_inbox: bool = True,
        max_concurrent_uploads: int = 1,
        threads: int = 1,
        enable_metrics: bool = True,
        background_tee: bool = False,
        update_db: bool = True,
    ):
        """
        Initialize the SubmissionProcessor with necessary configuration.

        This sets up the execution environment, including:
        - Loading the Crypt4GH private key for the inbox and public keys for archives.
        - Initializing the S3 connection pool based on thread concurrency settings.
        - Setting up the context for paired-end FASTQ consistency checks.
        - Initializing the progress logger for state persistence.

        :param configuration: Global processing configuration (DB, Archives, etc.).
        :param inbox: Specific inbox target configuration (S3 credentials, keys).
        :param status_file_path: Path to the local file used for tracking processing state.
        :param redact_patterns: Optional list of regex patterns to redact from metadata. (TODO actually use these)
        :param clean_inbox: Whether to remove files from the inbox after successful processing.
        :param max_concurrent_uploads: Number of threads used for S3 multipart uploads _per file_.
        :param threads: Number of files to process concurrently.
        :param enable_metrics: Whether to collect and log throughput/latency metrics.
        :param background_tee: Whether to use background threads for stream branching.
        :param update_db: Whether to update the central database state during processing.
        """
        self.config = configuration
        self._source_s3_options = inbox.s3

        log.debug("Loading crypt4gh keys...")
        self._private_key = Crypt4GH.retrieve_private_key(
            inbox.private_key_path, passphrase=inbox.private_key_passphrase
        )
        self._consented_pub_key = get_public_key(configuration.archives.consented.public_key_path)
        self._non_consented_pub_key = get_public_key(configuration.archives.non_consented.public_key_path)

        self._redact_patterns = redact_patterns or []
        self._enable_metrics = enable_metrics
        self._background_tee = background_tee
        self._update_db = update_db

        self.context = SubmissionContext()
        self._consistency = ConsistencyValidator(self.context)
        self.progress_logger = FileProgressLogger[ProcessingState](status_file_path)

        self._threads = threads
        self._max_concurrent_uploads = max_concurrent_uploads

        self._pool_size = max(10, threads * (1 + max_concurrent_uploads) + 1)
        log.debug(f"Configuring S3 client pool size: {self._pool_size}")

        self._source_s3 = init_s3_client(s3_options=self._source_s3_options, max_pool_connections=self._pool_size)
        self._target_s3_map: dict[tuple[AnyHttpUrl, str], S3Client] = {}

        self._submission_id: str | None = None
        self._should_qc: bool = False
        self._target_s3: Any = None
        self._target_bucket: str | None = None
        self._target_public_key: bytes | None = None
        self._partner_map: dict[str, str] = {}
        self._clean_inbox = clean_inbox

    def _get_target_s3(self, options: S3Options):
        """Get or create an S3 client for the target options."""
        key = (options.endpoint_url, options.bucket)
        if key not in self._target_s3_map:
            self._target_s3_map[key] = init_s3_client(s3_options=options, max_pool_connections=self._pool_size)
        return self._target_s3_map[key]

    def run(self, submission_metadata: SubmissionMetadata):
        """
        Execute the processing pipeline for a single submission.

        High-level view:
        1. Determine the target archive based on consent status (at the time of execution!).
        2. Check if submission should undergo detailed QC, in which case the decrypted files are written to local storage.
        3. Spawn threads to process files (Download -> Decrypt -> Validate -> Encrypt -> Archive).
        4. Archive redacted metadata.
        5. Optionally clean the inbox.

        :param submission_metadata: The parsed metadata object containing donor and file information.
        :raises RuntimeError: If consistency checks fail or any file fails validation.
        """
        self._submission_id = submission_metadata.content.submission_id
        is_consented = submission_metadata.content.consents_to_research(date.today())
        db = SubmissionDb(self.config.db.database_url, self.config.db.author)  # type: ignore[arg-type]
        target_percentage = self.config.detailed_qc.target_percentage
        self._should_qc = (
            db.should_qc(self._submission_id, target_percentage, self.config.detailed_qc.salt)
            if target_percentage > 0.0
            else False
        )
        if self._should_qc:
            # TODO: set selected_for_qc to TRUE for this submission in database
            log.info(f"Submission {self._submission_id} should be QCed.")

        target_archive = self.config.archives.consented if is_consented else self.config.archives.non_consented
        self._target_bucket = target_archive.s3.bucket
        self._target_s3 = self._get_target_s3(target_archive.s3)
        self._target_public_key = self._consented_pub_key if is_consented else self._non_consented_pub_key

        log.info(f"Consent Status: {'Consented' if is_consented else 'Non-Consented'}")
        log.info(f"Target Archive: {self._target_bucket}")

        files_map, thresholds, total_bytes = self._prepare_fastq_validation(submission_metadata)

        log.info(f"Processing {len(files_map)} files ({total_bytes / (1024**3):.2f} GB)...")

        with (
            tqdm(total=total_bytes, desc="Total     ", position=0, **TQDM_DEFAULTS) as pbar_global,  # type: ignore[call-overload]
            ThreadPoolExecutor(max_workers=self._threads) as pool,
        ):
            futures: list[Future] = [
                pool.submit(
                    self._process_one_file,
                    file_meta=file_meta,
                    threshold=thresholds.get(file_meta.file_path),
                    pbar_global=pbar_global,
                )
                for file_meta in files_map.values()
            ]
            for future in futures:
                future.result()

        if self.context.has_errors:
            raise RuntimeError("Submission failed consistency checks or validation.")

        self._upload_final_metadata(submission_metadata)
        log.info(f"Submission {self._submission_id} processed successfully.")

        if self._clean_inbox:
            with DbContext(
                self.config.model_dump(by_alias=True),
                self._submission_id,
                start_state=SubmissionStateEnum.CLEANING,
                end_state=SubmissionStateEnum.CLEANED,
                enabled=self._update_db,
            ):
                _clean_submission_from_bucket(
                    self._source_s3_options.bucket,
                    CleanConfig.model_validate({"s3": self._source_s3_options.model_dump(by_alias=True)}),
                    self._submission_id,
                )

        self._should_qc = False
        self._submission_id = None

    def _prepare_fastq_validation(
        self, submission_metadata: SubmissionMetadata
    ) -> tuple[dict[Path, File], dict[str, Thresholds], int]:
        files_map = submission_metadata.files
        total_bytes = sum(f.file_size_in_bytes for f in files_map.values())

        self._partner_map.clear()
        for _donor, _lab_datum, pairs, _thresholds in submission_metadata.iter_paired_end_fastqs():
            for fq1, fq2 in pairs:
                key1 = fq1.file_path
                key2 = fq2.file_path
                self._partner_map[key1] = key2
                self._partner_map[key2] = key1

        thresholds: dict[str, Thresholds] = {}
        for _, _, files, t in submission_metadata.iter_single_end_fastqs():
            for f in files:
                thresholds[f.file_path] = t
        for _, _, pairs, t in submission_metadata.iter_paired_end_fastqs():
            for fq1, fq2 in pairs:
                thresholds[fq1.file_path] = t
                thresholds[fq2.file_path] = t
        return files_map, thresholds, total_bytes

    def _process_one_file(self, file_meta: File, threshold: Thresholds | None, pbar_global: Any):
        src_key = f"{self._submission_id}/files/{file_meta.file_path}.c4gh"
        file_path_str = str(file_meta.file_path)
        file_name = file_path_str.rsplit("/", maxsplit=1)[-1]

        try:
            head = self._source_s3.head_object(Bucket=self._source_s3_options.bucket, Key=src_key)
            s3_size = head["ContentLength"]
            s3_mtime = head["LastModified"].timestamp()
        except Exception as e:
            log.error(f"Failed to access source file {src_key}: {e}")
            self.context.add_error(f"Source access failed: {src_key}")
            self.progress_logger.set_state(
                file_path_str,
                file_meta,
                {"processing_successful": False, "errors": [str(e)]},
                size=-1,
                mtime=-1.0,
            )
            return

        state = self.progress_logger.get_state(file_path_str, file_meta, size=s3_size, mtime=s3_mtime)
        if state and state.get("processing_successful"):
            log.info(f"Skipping {file_meta.file_path}, already processed.")
            pbar_global.update(file_meta.file_size_in_bytes)
            return

        if self.context.has_errors:
            return

        try:
            self._run_pipeline(file_meta, threshold, pbar_global, file_name)

            partner = self._partner_map.get(file_meta.file_path)
            if partner and not self._consistency.check_pair(file_meta.file_path, partner):
                raise RuntimeError(f"Consistency Check Failed: {file_meta.file_path}")

            self.progress_logger.set_state(
                file_path_str,
                file_meta,
                {"processing_successful": True, "errors": []},
                size=s3_size,
                mtime=s3_mtime,
            )

        except Exception as e:
            log.error(f"Failed processing {file_meta.file_path}: {e}")
            self.context.add_error(str(e))
            self.progress_logger.set_state(
                file_path_str,
                file_meta,
                {"processing_successful": False, "errors": [str(e)]},
                size=s3_size,
                mtime=s3_mtime,
            )

    def _run_pipeline(self, file_meta: File, threshold: Thresholds | None, pbar_global: Any, file_name: str):
        if not self._target_public_key:
            raise RuntimeError("Target public key not set.")
        if not self._target_bucket:
            raise RuntimeError("Target bucket not set.")
        if not self._submission_id:
            raise RuntimeError("Submission id not set.")

        src_key = f"{self._submission_id}/files/{file_meta.file_path}.c4gh"
        dest_key = f"{self._submission_id}/files/{file_meta.encrypted_file_path()}"

        part_size = calculate_s3_part_size(file_meta.file_size_in_bytes, None)
        metrics = StreamMetricsRegistry(enabled=self._enable_metrics)

        checksum_validator = ChecksumValidator(expected_checksum=file_meta.file_checksum)
        validation_chain = checksum_validator | metrics.measure("3a_Checksum")

        format_validator: ObserverWithMetrics | None = None
        if file_meta.file_type == FileType.fastq:
            format_validator = FastqValidator(
                mean_read_length_threshold=threshold.mean_read_length if threshold else None
            )
        elif file_meta.file_type == FileType.bam:
            format_validator = BamValidator()

        if format_validator:
            validation_chain |= format_validator | metrics.measure("3b_Format")

        with (
            tqdm(  # type: ignore[call-overload]
                total=file_meta.file_size_in_bytes,
                desc="Processing",
                postfix={"file": file_name},
                leave=False,
                **TQDM_DEFAULTS,
            ) as pbar_local,
            ExitStack() as stack,
        ):
            pipeline = (
                S3Downloader(self._source_s3, self._source_s3_options.bucket, src_key)
                | metrics.measure("1_Source")
                | Crypt4GHDecryptor(private_key=self._private_key)
                | metrics.measure("2_Decrypt")
            )

            if self._should_qc:
                path = Path(self.config.detailed_qc.local_storage) / self._submission_id / "files" / file_meta.file_path
                path.parent.mkdir(parents=True, exist_ok=True)

                writer = stack.enter_context(open(path, "wb"))
                writer = metrics.measure("2b_Write", writer)

                pipeline |= Tee(writer, threaded=self._background_tee)

            pipeline |= Tee(validation_chain, threaded=self._background_tee)

            pipeline = (
                pipeline | Crypt4GHEncryptor(recipient_pubkey=self._target_public_key) | metrics.measure("4_Encrypt")
            )

            pipeline |= Tee(TqdmObserver([pbar_global, pbar_local]), threaded=self._background_tee)

            uploader = S3MultipartUploader(
                self._target_s3,
                self._target_bucket,
                dest_key,
                part_size=part_size,
                max_threads=self._max_concurrent_uploads,
            )

            pipeline >> uploader

        # collect statistics
        stats = checksum_validator.metrics
        if format_validator:
            stats.update(format_validator.metrics)

        if metrics.enabled:
            log.info(f"Performance for {file_meta.file_path}: {metrics.report()}")

        self.context.record_stats(file_meta.file_path, stats)

    def _upload_final_metadata(self, submission_metadata: SubmissionMetadata):
        """Redacts and uploads the final metadata."""
        dest_key = f"{self._submission_id}/metadata/metadata.json"
        redacted_metadata = submission_metadata.content.to_redacted_dict()

        self._target_s3.put_object(
            Bucket=self._target_bucket, Key=dest_key, Body=BytesIO(str(redacted_metadata).encode("utf-8"))
        )
