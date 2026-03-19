import json
import logging
import tempfile
from concurrent.futures import Future, ThreadPoolExecutor
from contextlib import ExitStack
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import TYPE_CHECKING, Any

from crypt4gh.keys import get_public_key
from grz_common.constants import TQDM_DEFAULTS
from grz_common.pipeline.components import ObserverWithMetrics
from grz_common.workers.submission import SubmissionMetadata
from grz_db.models.submission import SubmissionDb, SubmissionStateEnum
from grz_pydantic_models.submission.metadata import File, FileType
from grz_pydantic_models.submission.thresholds import Thresholds
from grzctl.commands.clean import _clean_submission_from_bucket
from grzctl.dbcontext import DbContext
from grzctl.models.config import CleanConfig, InboxTarget, ProcessConfig
from tqdm.auto import tqdm

from ..models.s3 import S3Options
from ..progress import FileProgressLogger, ProcessingState
from ..transfer import init_s3_client
from ..utils.crypt import Crypt4GH
from ..utils.redaction import redact_file
from .components import Tee, TqdmObserver
from .components.crypt4gh import Crypt4GHDecryptor, Crypt4GHEncryptor
from .components.perf import StreamMetricsRegistry
from .components.s3 import S3Downloader, S3MultipartUploader, calculate_s3_part_size
from .components.validation import BamValidator, ChecksumValidator, FastqValidator
from .context import ReadPairConsistencyValidator, SubmissionContext

log = logging.getLogger(__name__)

if TYPE_CHECKING:
    from types_boto3_s3.client import S3Client
else:
    S3Client = Any


@dataclass
class SubmissionRunState:
    submission_metadata: SubmissionMetadata
    target_s3: S3Client
    target_bucket: str
    target_public_key: bytes
    should_qc: bool
    context: SubmissionContext = field(default_factory=SubmissionContext)
    consistency_validator: ReadPairConsistencyValidator = field(init=False)

    def __post_init__(self) -> None:
        self.consistency_validator = ReadPairConsistencyValidator.from_submission_metadata(
            self.context, self.submission_metadata
        )

    @property
    def submission_id(self) -> str:
        return self.submission_metadata.content.submission_id


class S3ClientCache:
    def __init__(self, pool_size: int):
        self._pool_size = pool_size
        self._target_s3_map: dict[tuple[str, str], S3Client] = {}

    def get(self, options: S3Options) -> S3Client:
        key = (str(options.endpoint_url), options.bucket)
        if key not in self._target_s3_map:
            self._target_s3_map[key] = init_s3_client(s3_options=options, max_pool_connections=self._pool_size)
        return self._target_s3_map[key]


class RunSetupCoordinator:
    def __init__(
        self,
        config: ProcessConfig,
        consented_pub_key: bytes,
        non_consented_pub_key: bytes,
        s3_client_cache: S3ClientCache,
    ):
        self._config = config
        self._consented_pub_key = consented_pub_key
        self._non_consented_pub_key = non_consented_pub_key
        self._s3_client_cache = s3_client_cache

    def _determine_qc_flag(self, submission_id: str) -> bool:
        db = SubmissionDb(self._config.db.database_url, self._config.db.author)  # type: ignore[arg-type]
        target_percentage = self._config.detailed_qc.target_percentage
        should_qc = (
            db.should_qc(submission_id, target_percentage, self._config.detailed_qc.salt)
            if target_percentage > 0.0
            else False
        )
        if should_qc:
            # TODO: set selected_for_qc to TRUE for this submission in database
            log.info(f"Submission {submission_id} should be QCed.")

        return should_qc

    def new_submission_run(self, submission_metadata: SubmissionMetadata) -> SubmissionRunState:
        is_research_consented = submission_metadata.content.consents_to_research(date.today())
        target_archive = (
            self._config.archives.consented if is_research_consented else self._config.archives.non_consented
        )

        run_state = SubmissionRunState(
            submission_metadata=submission_metadata,
            target_s3=self._s3_client_cache.get(target_archive.s3),
            target_bucket=target_archive.s3.bucket,
            target_public_key=self._consented_pub_key if is_research_consented else self._non_consented_pub_key,
            should_qc=self._determine_qc_flag(submission_metadata.content.submission_id),
        )

        log.info(f"Consent Status: {'Consented' if is_research_consented else 'Non-Consented'}")
        log.info(f"Target Archive: {run_state.target_bucket}")

        return run_state


class FilePipelineExecutor:
    def __init__(  # noqa: PLR0913
        self,
        source_s3: S3Client,
        source_bucket: str,
        private_key: bytes,
        progress_logger: FileProgressLogger[ProcessingState],
        threads: int,
        max_concurrent_uploads: int,
        enable_metrics: bool,
        background_tee: bool,
        qc_local_storage: str,
    ):
        self._source_s3 = source_s3
        self._source_bucket = source_bucket
        self._private_key = private_key
        self._progress_logger = progress_logger
        self._threads = threads
        self._max_concurrent_uploads = max_concurrent_uploads
        self._enable_metrics = enable_metrics
        self._background_tee = background_tee
        self._qc_local_storage = qc_local_storage

    @staticmethod
    def get_thresholds(submission_metadata: SubmissionMetadata) -> dict[str, Thresholds]:
        thresholds: dict[str, Thresholds] = {}
        for _, _, files, t in submission_metadata.iter_single_end_fastqs():
            for f in files:
                thresholds[f.file_path] = t
        for _, _, pairs, t in submission_metadata.iter_paired_end_fastqs():
            for fq1, fq2 in pairs:
                thresholds[fq1.file_path] = t
                thresholds[fq2.file_path] = t

        return thresholds

    def process_submission_files(
        self,
        run_state: SubmissionRunState,
    ) -> None:
        files_map = run_state.submission_metadata.files
        total_bytes = sum(f.file_size_in_bytes for f in files_map.values())
        thresholds = self.get_thresholds(run_state.submission_metadata)

        log.info(f"Processing {len(files_map)} files ({total_bytes / (1024**3):.2f} GB)...")

        with (
            tqdm(total=total_bytes, desc="Total     ", position=0, **TQDM_DEFAULTS) as pbar_global,  # type: ignore[call-overload]
            ThreadPoolExecutor(max_workers=self._threads) as pool,
        ):
            futures: list[Future] = [
                pool.submit(
                    self._process_one_file,
                    run_state=run_state,
                    file_meta=file_meta,
                    threshold=thresholds.get(file_meta.file_path),
                    pbar_global=pbar_global,
                )
                for file_meta in files_map.values()
            ]
            for future in futures:
                future.result()

    def _process_one_file(
        self,
        run_state: SubmissionRunState,
        file_meta: File,
        threshold: Thresholds | None,
        pbar_global: Any,
    ) -> None:
        src_key = f"{run_state.submission_id}/files/{file_meta.file_path}.c4gh"
        file_path_str = str(file_meta.file_path)
        file_name = file_path_str.rsplit("/", maxsplit=1)[-1]

        try:
            head = self._source_s3.head_object(Bucket=self._source_bucket, Key=src_key)
            s3_size = head["ContentLength"]
            s3_mtime = head["LastModified"].timestamp()
        except Exception as e:
            log.exception(
                "Failed to access source file",
                extra={
                    "submission_id": run_state.submission_id,
                    "file_path": file_meta.file_path,
                    "src_key": src_key,
                },
            )
            run_state.context.add_error(f"Source access failed: {src_key}")
            self._progress_logger.set_state(
                file_path_str,
                file_meta,
                {"processing_successful": False, "errors": [str(e)]},
                size=-1,
                mtime=-1.0,
            )
            return

        state = self._progress_logger.get_state(file_path_str, file_meta, size=s3_size, mtime=s3_mtime)
        if state and state.get("processing_successful"):
            log.info(f"Skipping {file_meta.file_path}, already processed.")
            pbar_global.update(file_meta.file_size_in_bytes)
            return

        if run_state.context.has_errors:
            return

        try:
            self._run_pipeline(run_state, file_meta, threshold, pbar_global, file_name)

            if not run_state.consistency_validator.check(file_meta.file_path):
                raise RuntimeError(f"Consistency Check Failed: {file_meta.file_path}")

            self._progress_logger.set_state(
                file_path_str,
                file_meta,
                {"processing_successful": True, "errors": []},
                size=s3_size,
                mtime=s3_mtime,
            )

        except Exception as e:
            log.exception(
                "Failed processing file",
                extra={
                    "submission_id": run_state.submission_id,
                    "file_path": file_meta.file_path,
                    "src_key": src_key,
                },
            )
            run_state.context.add_error(str(e))
            self._progress_logger.set_state(
                file_path_str,
                file_meta,
                {"processing_successful": False, "errors": [str(e)]},
                size=s3_size,
                mtime=s3_mtime,
            )

    @staticmethod
    def build_format_validator(
        file_meta: File,
        threshold: Thresholds | None,
    ) -> ObserverWithMetrics | None:
        format_validator: ObserverWithMetrics | None = None
        if file_meta.file_type == FileType.fastq:
            format_validator = FastqValidator(
                mean_read_length_threshold=threshold.mean_read_length if threshold else None
            )
        elif file_meta.file_type == FileType.bam:
            format_validator = BamValidator()

        return format_validator

    def _run_pipeline(
        self,
        run_state: SubmissionRunState,
        file_meta: File,
        threshold: Thresholds | None,
        pbar_global: Any,
        file_name: str,
    ) -> None:
        if not run_state.target_public_key:
            raise RuntimeError("Target public key not set.")
        if not run_state.target_bucket:
            raise RuntimeError("Target bucket not set.")
        if run_state.target_s3 is None:
            raise RuntimeError("Target S3 client not set.")

        src_key = f"{run_state.submission_id}/files/{file_meta.file_path}.c4gh"
        dest_key = f"{run_state.submission_id}/files/{file_meta.encrypted_file_path()}"

        metrics = StreamMetricsRegistry(enabled=self._enable_metrics)

        checksum_validator = ChecksumValidator(expected_checksum=file_meta.file_checksum)
        format_validator = self.build_format_validator(
            file_meta=file_meta,
            threshold=threshold,
        )

        validation_chain = checksum_validator | metrics.measure("3a_Checksum")
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
            # download and decrypt
            pipeline = (
                S3Downloader(self._source_s3, self._source_bucket, src_key)
                | metrics.measure("1_Source")
                | Crypt4GHDecryptor(private_key=self._private_key)
                | metrics.measure("2_Decrypt")
            )

            # tee to local folder if QC is needed
            if run_state.should_qc:
                path = Path(self._qc_local_storage) / run_state.submission_id / "files" / file_meta.file_path
                path.parent.mkdir(parents=True, exist_ok=True)

                writer = stack.enter_context(open(path, "wb"))
                writer_observer: Any = metrics.measure("2b_Write", writer)
                pipeline |= Tee(writer_observer, threaded=self._background_tee)

            # add validation chain
            pipeline |= Tee(validation_chain, threaded=self._background_tee)

            # re-encrypt
            pipeline = (
                pipeline
                | Crypt4GHEncryptor(recipient_pubkey=run_state.target_public_key)
                | metrics.measure("4_Encrypt")
            )

            # progress bar
            pipeline |= Tee(TqdmObserver([pbar_global, pbar_local]), threaded=self._background_tee)

            # upload to archive bucket
            part_size = calculate_s3_part_size(file_meta.file_size_in_bytes, None)
            uploader = S3MultipartUploader(
                run_state.target_s3,
                run_state.target_bucket,
                dest_key,
                part_size=part_size,
                max_threads=self._max_concurrent_uploads,
            )

            # run the whole pipeline
            pipeline >> uploader

        stats = checksum_validator.metrics
        if format_validator:
            stats.update(format_validator.metrics)

        if metrics.enabled:
            log.info(f"Performance for {file_meta.file_path}: {metrics.report()}")

        run_state.context.record_stats(file_meta.file_path, stats)


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
        :param clean_inbox: Whether to remove files from the inbox after successful processing.
        :param max_concurrent_uploads: Number of threads used for S3 multipart uploads _per file_.
        :param threads: Number of files to process concurrently.
        :param enable_metrics: Whether to collect and log throughput/latency metrics.
        :param background_tee: Whether to use background threads for stream branching.
        :param update_db: Whether to update the central database state during processing.
        """
        self.config = configuration
        self.inbox = inbox
        self._source_s3_options = inbox.s3
        self._update_db = update_db
        self._clean_inbox = clean_inbox
        self._log_dir = status_file_path.parent

        log.debug("Loading crypt4gh keys...")

        s3_pool_size = max(10, threads * (1 + max_concurrent_uploads) + 1)
        log.debug(f"Configuring S3 client pool size: {s3_pool_size}")

        self._setup = RunSetupCoordinator(
            config=self.config,
            consented_pub_key=get_public_key(configuration.archives.consented.public_key_path),
            non_consented_pub_key=get_public_key(configuration.archives.non_consented.public_key_path),
            s3_client_cache=S3ClientCache(pool_size=s3_pool_size),
        )
        self._pipeline_executor = FilePipelineExecutor(
            source_s3=init_s3_client(s3_options=self._source_s3_options, max_pool_connections=s3_pool_size),
            source_bucket=self._source_s3_options.bucket,
            private_key=Crypt4GH.retrieve_private_key(inbox.private_key_path, passphrase=inbox.private_key_passphrase),
            progress_logger=FileProgressLogger[ProcessingState](status_file_path),
            threads=threads,
            max_concurrent_uploads=max_concurrent_uploads,
            enable_metrics=enable_metrics,
            background_tee=background_tee,
            qc_local_storage=self.config.detailed_qc.local_storage,
        )

    def _upload_final_metadata(self, submission_metadata: SubmissionMetadata, run_state: SubmissionRunState) -> None:
        dest_key = f"{run_state.submission_id}/metadata/metadata.json"
        redacted_metadata = submission_metadata.content.to_redacted_dict()

        run_state.target_s3.put_object(
            Bucket=run_state.target_bucket,
            Key=dest_key,
            Body=json.dumps(redacted_metadata).encode("utf-8"),
        )

    def _maybe_cleanup_inbox(self, run_state: SubmissionRunState) -> None:
        if not self._clean_inbox:
            return

        with DbContext(
            self.config.model_dump(by_alias=True),
            run_state.submission_id,
            start_state=SubmissionStateEnum.CLEANING,
            end_state=SubmissionStateEnum.CLEANED,
            enabled=self._update_db,
        ):
            _clean_submission_from_bucket(
                self._source_s3_options.bucket,
                CleanConfig.model_validate({"s3": self._source_s3_options.model_dump(by_alias=True)}),
                run_state.submission_id,
            )

    def _upload_redacted_logs(self, submission_metadata: SubmissionMetadata, run_state: SubmissionRunState) -> None:
        if not self._log_dir.exists():
            return

        redaction_patterns = submission_metadata.content.create_redaction_patterns()
        for file_path in sorted(self._log_dir.rglob("*")):
            if not file_path.is_file():
                continue

            key_suffix = file_path.relative_to(self._log_dir).as_posix()
            dest_key = f"{run_state.submission_id}/logs/{key_suffix}"

            if redaction_patterns:
                with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as tmp_file:
                    tmp_path = Path(tmp_file.name)
                    redact_file(file_path, tmp_path, redaction_patterns)
                    body = tmp_path.read_bytes()
            else:
                body = file_path.read_bytes()

            run_state.target_s3.put_object(Bucket=run_state.target_bucket, Key=dest_key, Body=body)

    def run(self, submission_metadata: SubmissionMetadata) -> None:
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
        submission_run = self._setup.new_submission_run(submission_metadata)
        self._pipeline_executor.process_submission_files(submission_run)

        if submission_run.context.has_errors:
            raise RuntimeError("Submission failed consistency checks or validation.")

        self._upload_final_metadata(submission_metadata, submission_run)
        self._upload_redacted_logs(submission_metadata, submission_run)
        log.info(f"Submission {submission_run.submission_id} processed successfully.")
        self._maybe_cleanup_inbox(submission_run)
