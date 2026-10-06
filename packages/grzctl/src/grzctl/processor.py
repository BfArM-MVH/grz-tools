import io
import json
import logging
import stat
import subprocess
import tempfile
from concurrent.futures import Future, ThreadPoolExecutor
from contextlib import ExitStack, suppress
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

from boto3.s3.transfer import TransferConfig
from botocore.exceptions import ClientError
from grz_common.constants import TQDM_DEFAULTS
from grz_common.exceptions import (
    DetailedQCError,
    DownloadError,
    MissingObjectError,
    MissingSubmissionFileError,
    UploadError,
)
from grz_common.pipeline.components import (
    DevNullSink,
    ObserverWithMetrics,
    PipelineError,
    ReadStream,
    Tee,
    TqdmObserver,
    WriteStream,
)
from grz_common.pipeline.components.crypt4gh import Crypt4GHDecryptor, Crypt4GHEncryptor
from grz_common.pipeline.components.perf import StreamMetricsRegistry
from grz_common.pipeline.components.s3 import S3Downloader, S3MultipartUploader, calculate_s3_part_size
from grz_common.pipeline.components.validation import (
    BamValidator,
    ChecksumValidator,
    FastqValidator,
)
from grz_common.pipeline.context import ReadPairConsistencyValidator, SubmissionContext
from grz_common.progress import FileProgressLogger, ProcessingState
from grz_common.transfer import head_object, init_s3_client, s3_errors
from grz_common.utils.redaction import redact_file
from grz_common.workers.submission import SubmissionMetadata
from grz_db.errors import DuplicateInitialSubmissionError
from grz_db.models.submission import SubmissionDb
from grz_pydantic_models.submission.metadata import File, FileType
from grz_pydantic_models.submission.thresholds import Thresholds
from tqdm.auto import tqdm

from .dbcontext import FilesFailedError
from .models.config import GrzctlConfig, InboxTarget


def _close_undriven(stage: WriteStream) -> None:
    """Close a stage that the pipeline may never have driven, ignoring what it reports.

    Closing joins the validator's thread, which is the point of calling it. A stage the pipeline
    drove was closed there already and reports nothing new. A stage it never reached saw no data,
    so what it has to report is about that rather than about the file.
    """
    with suppress(PipelineError):
        stage.close()


log = logging.getLogger(__name__)

if TYPE_CHECKING:
    from types_boto3_s3.client import S3Client
else:
    S3Client = Any


def _inbox_file_key(submission_id: str, file_meta: File) -> str:
    """S3 key of an encrypted file in the inbox."""
    return f"{submission_id}/files/{file_meta.encrypted_file_path()}"


def _archive_file_key(submission_id: str, file_meta: File) -> str:
    """S3 key of an encrypted file in the interrogation bucket and the archives."""
    return f"{submission_id}/files/{file_meta.encrypted_file_path()}"


def _archive_metadata_key(submission_id: str) -> str:
    """S3 key of the redacted metadata in the interrogation bucket and the archives."""
    return f"{submission_id}/metadata/metadata.json"


def _qc_file_path(local_storage: str, submission_id: str, file_meta: File) -> Path:
    """Local path of a decrypted file for detailed QC."""
    return Path(local_storage) / submission_id / "files" / file_meta.file_path


QC_DIRECTORY_MODE = 0o770
"""Mode of the directories below ``detailed_qc.local_storage``, which hold decrypted files and the unredacted metadata."""


def _ensure_directory_mode(local_storage: Path, directory: Path, mode: int) -> None:
    """Create ``directory`` below ``local_storage``, and give it and its parents below ``local_storage`` the mode ``mode``.

    ``mkdir`` applies the umask and leaves an existing directory as it is, so a directory with another mode gets ``chmod``.
    """
    directory.mkdir(parents=True, exist_ok=True)
    path = local_storage
    for part in directory.relative_to(local_storage).parts:
        path /= part
        if stat.S_IMODE(path.stat().st_mode) != mode:
            path.chmod(mode)


def _copy_object(
    s3_client: S3Client, source_bucket: str, target_bucket: str, key: str, preferred_part_size: int
) -> None:
    """Copy an object between buckets in at most ``MULTIPART_MAX_PARTS`` parts.

    Without a ``TransferConfig``, boto3 copies in 8 MiB parts and allows up to 10000 of them.
    An object above 7.81 GiB would then exceed the 1000 parts that Ceph/Bluestore/Quobyte allow.
    """
    size = head_object(s3_client, source_bucket, key)["ContentLength"]
    config = TransferConfig(multipart_chunksize=calculate_s3_part_size(size, preferred_part_size))
    s3_client.copy(CopySource={"Bucket": source_bucket, "Key": key}, Bucket=target_bucket, Key=key, Config=config)


def _put_object(s3_client: S3Client, bucket: str, key: str, body: bytes) -> None:
    """Upload a small object through ``S3MultipartUploader``, which checks its MD5 as it does for the staged files."""
    ReadStream(io.BytesIO(body)) >> S3MultipartUploader(s3_client, bucket, key)


@dataclass
class SubmissionRunState:
    """Holds all the state needed to run the streaming pipeline for one submission.

    Encapsulates the resolved S3 clients/buckets for the interrogation (staging)
    and final archive, the target re-encryption key, and the per-submission
    pipeline ``context`` that components share as they stream files through.
    """

    submission_metadata: SubmissionMetadata
    interrogation_s3: S3Client
    interrogation_bucket: str
    interrogation_part_size: int
    final_s3: S3Client
    final_bucket: str
    final_part_size: int
    target_public_key: bytes
    context: SubmissionContext = field(default_factory=SubmissionContext)
    consistency_validator: ReadPairConsistencyValidator = field(init=False)

    def __post_init__(self) -> None:
        partner_map = ReadPairConsistencyValidator.get_partner_map(self.submission_metadata)
        self.consistency_validator = ReadPairConsistencyValidator(self.context, partner_map)

    @property
    def submission_id(self) -> str:
        return self.submission_metadata.content.submission_id


class FilePipelineExecutor:
    """Streams a submission's files from the inbox into their outputs.

    Each file has two outputs, each with its own progress log:

    - its validated, re-encrypted copy, staged in the interrogation bucket (``staging_log``)
    - its decrypted copy on local storage for detailed QC (``local_log``)

    An output counts as done while its log records it and its copy is still there. A call to
    :meth:`process_files` streams each file at most once, into the outputs it is missing.
    """

    def __init__(  # noqa: PLR0913, PLR0917
        self,
        source_s3: S3Client,
        source_bucket: str,
        private_key: bytes,
        sender_private_key: bytes,
        threads: int,
        max_concurrent_uploads: int,
        qc_local_storage: str,
        staging_log: FileProgressLogger[ProcessingState],
        local_log: FileProgressLogger[ProcessingState],
    ):
        self._source_s3 = source_s3
        self._source_bucket = source_bucket
        self._private_key = private_key
        self._sender_private_key = sender_private_key
        self._threads = threads
        self._max_concurrent_uploads = max_concurrent_uploads
        self._qc_local_storage = qc_local_storage
        self._staging_log = staging_log
        self._local_log = local_log

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

    def process_files(self, run_state: SubmissionRunState, stage: bool, write_local: bool) -> None:
        """Stream each file of the submission into the requested outputs it is missing.

        :param run_state: The submission's run state.
        :param stage: Validate each file and stage its re-encrypted copy in the interrogation bucket.
        :param write_local: Write each file's decrypted copy to local storage.
        """
        files_map = run_state.submission_metadata.files
        total_bytes = sum(f.file_size_in_bytes for f in files_map.values())
        thresholds = self.get_thresholds(run_state.submission_metadata)

        log.info(f"Processing {len(files_map)} files ({total_bytes / (1024**3):.2f} GB)...")
        with (
            tqdm(total=total_bytes, desc="Total     ", position=0, **TQDM_DEFAULTS) as pbar_global,  # type: ignore[call-overload]
            ThreadPoolExecutor(max_workers=self._threads) as pool,
        ):
            try:
                futures: list[Future] = [
                    pool.submit(
                        self._process_file,
                        run_state=run_state,
                        file_meta=file_meta,
                        threshold=thresholds.get(file_meta.file_path),
                        pbar_global=pbar_global,
                        stage=stage,
                        write_local=write_local,
                    )
                    for file_meta in files_map.values()
                ]
                for future in futures:
                    future.result()
            except BaseException:
                # Ctrl-C and SIGTERM reach only this thread: the running files finish, the queued ones never start
                pool.shutdown(wait=True, cancel_futures=True)
                raise

    def _process_file(  # noqa: PLR0913, PLR0917
        self,
        run_state: SubmissionRunState,
        file_meta: File,
        threshold: Thresholds | None,
        pbar_global: Any,
        stage: bool,
        write_local: bool,
    ) -> None:
        """Stream one file of the submission into the outputs it is still missing.

        Whatever goes wrong here fails this file alone. The error goes to the run's context, which
        fails the run once every file is done and lets the other threads stop early. Only a file that
        finished correctly is marked completed. The read-pair check compares a pair once both files
        are completed, so it leaves the partner of a failed file unchecked.

        :param run_state: The submission's run state.
        :param file_meta: The file's entry in the submission metadata.
        :param threshold: The QC thresholds that apply to the file, if it has any.
        :param pbar_global: The progress bar over all files of the submission.
        :param stage: Validate the file and stage its re-encrypted copy in the interrogation bucket.
        :param write_local: Write the file's decrypted copy to local storage.
        """
        inbox_key = _inbox_file_key(run_state.submission_id, file_meta)
        file_path_str = str(file_meta.file_path)

        # the size and the modification time of the inbox object tell a progress log record which
        # version of the file it describes
        try:
            head = head_object(
                self._source_s3, self._source_bucket, inbox_key, missing_error=MissingSubmissionFileError
            )
            s3_size = head["ContentLength"]
            s3_mtime = head["LastModified"].timestamp()
        except Exception as e:
            log.exception(
                "Failed to access source file",
                extra={
                    "submission_id": run_state.submission_id,
                    "file_path": file_meta.file_path,
                    "src_key": inbox_key,
                },
            )
            run_state.context.add_error(file_path_str, e)
            # the file could not be read, so neither its size nor its modification time is known
            failure: ProcessingState = {"processing_successful": False, "errors": [str(e)]}
            self._record(file_meta, failure, size=-1, mtime=-1.0, staging=stage, local=write_local)
            return

        # the except clause below records the failure for these outputs, also when a lookup fails
        needs_staging = stage
        needs_local_copy = write_local

        try:
            # an earlier run may have written some of this file's outputs already
            staging_record = self._staging_record(run_state, file_meta, s3_size, s3_mtime) if stage else None
            if staging_record is not None:
                # the read-pair check of this file's partner compares these stats
                run_state.context.record_stats(file_path_str, staging_record.get("stats", {}))
            needs_staging = stage and staging_record is None
            needs_local_copy = write_local and not self._has_local_copy(run_state, file_meta, s3_size, s3_mtime)

            if not needs_staging and not needs_local_copy:
                # every output this run asks for is there, so only the progress bar moves
                log.info(f"Skipping {file_meta.file_path}, already processed.")
                with TqdmObserver.lock:
                    pbar_global.update(file_meta.file_size_in_bytes)
                run_state.context.mark_completed(file_path_str)
                return

            if run_state.context.has_errors:
                # another file failed and the run fails with it, so this one is not transferred
                return

            # download, decrypt, validate, and write the outputs that are missing
            self._stream_file(
                run_state,
                file_meta,
                inbox_key,
                threshold,
                pbar_global,
                stage=needs_staging,
                write_local=needs_local_copy,
            )

            if needs_staging:
                # a read pair is compared once both files are marked completed, so whichever of the
                # two finishes second is the one that reports a mismatch
                run_state.consistency_validator.check(file_meta.file_path)

            success: ProcessingState = {
                "processing_successful": True,
                "errors": [],
                "stats": run_state.context.get_stats(file_path_str),
            }
            self._record(
                file_meta, success, size=s3_size, mtime=s3_mtime, staging=needs_staging, local=needs_local_copy
            )
            run_state.context.mark_completed(file_path_str)

        except BaseException as e:
            # grz-check reports a Rust panic as pyo3's PanicException, which is no Exception.
            # A worker thread never receives KeyboardInterrupt, so this catches no interrupt.
            log.exception(
                "Failed processing file",
                extra={
                    "submission_id": run_state.submission_id,
                    "file_path": file_meta.file_path,
                    "src_key": inbox_key,
                },
            )
            run_state.context.add_error(file_path_str, e)
            # the outputs this run was going to write are recorded as failed, so a rerun redoes them
            failure = {"processing_successful": False, "errors": [str(e)]}
            self._record(
                file_meta, failure, size=s3_size, mtime=s3_mtime, staging=needs_staging, local=needs_local_copy
            )

    def _record(  # noqa: PLR0913
        self, file_meta: File, state: ProcessingState, *, size: int, mtime: float, staging: bool, local: bool
    ) -> None:
        """Record ``state`` of a file in the progress logs of the selected outputs."""
        for progress_log, selected in ((self._staging_log, staging), (self._local_log, local)):
            if selected:
                progress_log.set_state(str(file_meta.file_path), file_meta, state, size=size, mtime=mtime)

    def _staging_record(
        self, run_state: SubmissionRunState, file_meta: File, s3_size: int, s3_mtime: float
    ) -> ProcessingState | None:
        """Return the staging record of a file whose staged copy is still in the interrogation bucket."""
        state = self._staging_log.get_state(str(file_meta.file_path), file_meta, size=s3_size, mtime=s3_mtime)
        if not state or not state.get("processing_successful"):
            return None
        try:
            head_object(
                run_state.interrogation_s3,
                run_state.interrogation_bucket,
                _archive_file_key(run_state.submission_id, file_meta),
            )
        except MissingObjectError:
            log.info(f"The staged copy of {file_meta.file_path} is gone, staging it again.")
            return None
        return state

    def _has_local_copy(self, run_state: SubmissionRunState, file_meta: File, s3_size: int, s3_mtime: float) -> bool:
        """Whether a file's decrypted copy is recorded and still on local storage."""
        state = self._local_log.get_state(str(file_meta.file_path), file_meta, size=s3_size, mtime=s3_mtime)
        if not state or not state.get("processing_successful"):
            return False
        return _qc_file_path(self._qc_local_storage, run_state.submission_id, file_meta).is_file()

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

    def _stream_file(  # noqa: PLR0913, PLR0917
        self,
        run_state: SubmissionRunState,
        file_meta: File,
        inbox_key: str,
        threshold: Thresholds | None,
        pbar_global: Any,
        stage: bool,
        write_local: bool,
    ) -> None:
        """Download and decrypt one file, check its checksum, and write it into the requested outputs.

        :param stage: Also validate the file's format and stage its re-encrypted copy in the interrogation bucket.
        :param write_local: Also write the decrypted copy to local storage.
        """
        metrics = StreamMetricsRegistry()

        with (
            tqdm(  # type: ignore[call-overload]
                total=file_meta.file_size_in_bytes,
                desc="Processing",
                postfix={"file": str(file_meta.file_path).rsplit("/", maxsplit=1)[-1]},
                leave=False,
                **TQDM_DEFAULTS,
            ) as pbar_local,
            ExitStack() as stack,
        ):
            # download and decrypt
            source = S3Downloader(
                self._source_s3, self._source_bucket, inbox_key, missing_error=MissingSubmissionFileError
            )

            # A validator runs a thread from the moment it is built, so it is built once the
            # download has opened and registered for closing straight away. Leaving the run before
            # the pipeline drives it would otherwise hold that thread until the process exits.
            checksum_validator = ChecksumValidator(expected_checksum=file_meta.file_checksum)
            format_validator = self.build_format_validator(file_meta=file_meta, threshold=threshold) if stage else None
            validation_chain = checksum_validator | metrics.measure("3a_Checksum")
            if format_validator:
                validation_chain |= format_validator | metrics.measure("3b_Format")
            stack.callback(_close_undriven, validation_chain)

            pipeline = (
                source
                | metrics.measure("1_Source")
                | Crypt4GHDecryptor(private_key=self._private_key)
                | metrics.measure("2_Decrypt")
            )

            # tee to local storage for detailed QC
            if write_local:
                path = _qc_file_path(self._qc_local_storage, run_state.submission_id, file_meta)
                _ensure_directory_mode(Path(self._qc_local_storage), path.parent, QC_DIRECTORY_MODE)
                writer = stack.enter_context(open(path, "wb"))
                pipeline |= Tee(metrics.measure("2b_Write")(writer))

            pipeline |= Tee(validation_chain)

            if not stage:
                pipeline |= Tee(TqdmObserver([pbar_global, pbar_local]))
                pipeline >> DevNullSink()
            else:
                # re-encrypt
                pipeline = (
                    pipeline
                    | Crypt4GHEncryptor(
                        recipient_pubkey=run_state.target_public_key,
                        sender_privkey=self._sender_private_key,
                    )
                    | metrics.measure("4_Encrypt")
                )
                pipeline |= Tee(TqdmObserver([pbar_global, pbar_local]))

                # upload to the interrogation bucket, under the file's archive key.
                # Size the parts by the inbox object: the re-encrypted object has the same payload and a
                # one-packet header, the smallest header an inbox object can have, so it is never larger.
                uploader = S3MultipartUploader(
                    run_state.interrogation_s3,
                    run_state.interrogation_bucket,
                    _archive_file_key(run_state.submission_id, file_meta),
                    part_size=calculate_s3_part_size(source.length, run_state.interrogation_part_size),
                    max_threads=self._max_concurrent_uploads,
                    content_type="application/octet-stream",
                )
                pipeline >> uploader

        stats = checksum_validator.metrics
        if format_validator:
            stats.update(format_validator.metrics)

        log.info(f"Performance for {file_meta.file_path}: {metrics.report()}")

        run_state.context.record_stats(str(file_meta.file_path), stats)


class SubmissionProcessor:
    """
    Orchestrates the streaming submission pipeline:
    Inbox -> Decrypt -> Validate -> Re-Encrypt -> Archive
    """

    def __init__(
        self,
        configuration: GrzctlConfig,
        inbox: InboxTarget,
        log_dir: Path,
        max_concurrent_uploads: int = 1,
        threads: int = 1,
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
        :param log_dir: Directory for the progress logs, which are archived with the submission.
        :param max_concurrent_uploads: Number of threads used for S3 multipart uploads _per file_.
        :param threads: Number of files to process concurrently.
        """
        self.config = configuration
        self.inbox = inbox
        self._source_s3_options = inbox.s3
        self._log_dir = log_dir

        log.debug("Loading crypt4gh keys...")
        # the pipeline stages work with the raw 32 bytes of each key
        self._consented_pub_key = configuration.archives.consented.load_public_key().public_bytes_raw()
        self._non_consented_pub_key = configuration.archives.non_consented.load_public_key().public_bytes_raw()

        self._s3_pool_size = max(10, threads * (1 + max_concurrent_uploads) + 1)
        log.debug(f"Configuring S3 client pool size: {self._s3_pool_size}")

        # The inbox key decrypts the files, and it also signs their re-encrypted copies,
        # as the step-by-step ``encrypt`` command does.
        inbox_private_key = inbox.load_private_key().private_bytes_raw()

        self._pipeline_executor = FilePipelineExecutor(
            source_s3=init_s3_client(s3_options=self._source_s3_options, max_pool_connections=self._s3_pool_size),
            source_bucket=self._source_s3_options.bucket,
            private_key=inbox_private_key,
            sender_private_key=inbox_private_key,
            threads=threads,
            max_concurrent_uploads=max_concurrent_uploads,
            qc_local_storage=self.config.detailed_qc.local_storage,
            staging_log=FileProgressLogger[ProcessingState](log_dir / "progress_staging.cjson"),
            local_log=FileProgressLogger[ProcessingState](log_dir / "progress_local.cjson"),
        )

    def _new_run_state(self, submission_metadata: SubmissionMetadata) -> SubmissionRunState:
        """Resolve the archive and the re-encryption key from the consent status."""
        submission_date = submission_metadata.content.submission.submission_date
        is_research_consented = submission_metadata.content.consents_to_research(submission_date)
        target_archive = self.config.archives.consented if is_research_consented else self.config.archives.non_consented
        interrogation_archive = self.config.archives.interrogation

        run_state = SubmissionRunState(
            submission_metadata=submission_metadata,
            interrogation_s3=init_s3_client(interrogation_archive.s3, max_pool_connections=self._s3_pool_size),
            interrogation_bucket=interrogation_archive.s3.bucket,
            interrogation_part_size=interrogation_archive.s3.multipart_chunksize,
            final_s3=init_s3_client(target_archive.s3, max_pool_connections=self._s3_pool_size),
            final_bucket=target_archive.s3.bucket,
            final_part_size=target_archive.s3.multipart_chunksize,
            target_public_key=self._consented_pub_key if is_research_consented else self._non_consented_pub_key,
        )

        log.info(f"Consent Status: {'Consented' if is_research_consented else 'Non-Consented'}")
        log.info(f"Target Archive: {run_state.final_bucket} (via Interrogation: {run_state.interrogation_bucket})")

        return run_state

    def _qc_selection_enabled(self) -> bool:
        return self.config.detailed_qc.target_percentage > 0.0

    def _predict_qc(self, db: SubmissionDb, submission_id: str) -> bool:
        """Guess the detailed QC decision before basic QC has passed."""
        if not self._qc_selection_enabled():
            return False
        detailed_qc = self.config.detailed_qc
        try:
            return db.should_qc(submission_id, detailed_qc.target_percentage, detailed_qc.salt, predict=True)
        except Exception as e:
            log.warning(f"Could not predict detailed QC for {submission_id}, not prefetching: {e}")
            return False

    def _determine_qc_flag(self, db: SubmissionDb, submission_id: str) -> bool:
        detailed_qc = self.config.detailed_qc
        should_qc = self._qc_selection_enabled() and db.should_qc(
            submission_id, detailed_qc.target_percentage, detailed_qc.salt
        )
        if should_qc:
            log.info(f"Submission {submission_id} selected for detailed QC.")

        return should_qc

    def _discard_qc_prefetch(self, run_state: SubmissionRunState) -> None:
        """Delete the decrypted files that the main pass wrote to local storage."""
        log.info(f"Deleting the files prefetched for detailed QC of {run_state.submission_id}...")
        for file_meta in run_state.submission_metadata.files.values():
            path = _qc_file_path(self.config.detailed_qc.local_storage, run_state.submission_id, file_meta)
            try:
                path.unlink(missing_ok=True)
            except OSError as e:
                log.warning(f"Could not delete prefetched QC file {path}: {e}")

    def _upload_final_metadata(self, submission_metadata: SubmissionMetadata, run_state: SubmissionRunState) -> None:
        redacted_metadata = submission_metadata.content.to_redacted_dict()

        key = _archive_metadata_key(run_state.submission_id)
        body = json.dumps(redacted_metadata).encode("utf-8")
        _put_object(run_state.interrogation_s3, run_state.interrogation_bucket, key, body)

    def _log_files(self, submission_id: str) -> dict[str, Path]:
        """Map the archive key of each local log file to its path."""
        if not self._log_dir.exists():
            return {}

        return {
            f"{submission_id}/logs/{file_path.relative_to(self._log_dir).as_posix()}": file_path
            for file_path in sorted(self._log_dir.rglob("*"))
            if file_path.is_file()
        }

    def _upload_redacted_logs(self, submission_metadata: SubmissionMetadata, run_state: SubmissionRunState) -> None:
        redaction_patterns = submission_metadata.content.create_redaction_patterns()
        for dest_key, file_path in self._log_files(run_state.submission_id).items():
            if redaction_patterns:
                with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as tmp_file:
                    tmp_path = Path(tmp_file.name)
                    redact_file(file_path, tmp_path, redaction_patterns)
                    body = tmp_path.read_bytes()
            else:
                body = file_path.read_bytes()

            _put_object(run_state.interrogation_s3, run_state.interrogation_bucket, dest_key, body)

    def _get_expected_keys(self, run_state: SubmissionRunState) -> list[str]:
        """The keys of a submission in the order they are archived, with the metadata last.

        An archived metadata object marks a finished archival, as it does for ``grzctl archive``.
        """
        # a set, because lab data may share a file such as target_regions.bed
        file_keys = {
            _archive_file_key(run_state.submission_id, f) for f in run_state.submission_metadata.files.values()
        }
        log_keys = list(self._log_files(run_state.submission_id))
        return [*sorted(file_keys), *log_keys, _archive_metadata_key(run_state.submission_id)]

    @staticmethod
    def _is_archived(run_state: SubmissionRunState) -> bool:
        """Whether the target archive holds the submission's metadata, and log a warning if it does.

        The archive copy ends with the metadata, so an archived metadata object means an earlier run finished.
        If S3 denies access to it, the submission counts as not archived, as it does for ``grzctl archive``.
        """
        key = _archive_metadata_key(run_state.submission_id)
        try:
            head_object(run_state.final_s3, run_state.final_bucket, key)
        except MissingObjectError:
            return False
        except DownloadError as e:
            cause = e.__cause__
            if not (isinstance(cause, ClientError) and cause.response["Error"]["Code"] == "AccessDenied"):
                raise
            log.warning(
                f"Cannot check whether submission '{run_state.submission_id}' is already archived, "
                f"because S3 denies access to '{key}'. Archiving it anyway."
            )
            return False
        log.warning(
            f"Submission '{run_state.submission_id}' is already archived in '{run_state.final_bucket}'. "
            "Nothing is streamed or copied again."
        )
        return True

    def _commit_to_archive(self, run_state: SubmissionRunState) -> None:
        """Copy all staged files from the interrogation bucket to the final archive.

        After a successful copy the source files are deleted from the interrogation
        bucket.  If the copy fails midway, the exception propagates and the caller
        is responsible for cleaning up the interrogation bucket (via
        ``_handle_interrogation_failure``).

        .. warning:: A partial copy failure leaves some objects in the final archive
           that cannot be automatically rolled back.  The operator must remove them
           manually.
        """
        expected_keys = self._get_expected_keys(run_state)
        log.info(f"Copying {len(expected_keys)} files from interrogation bucket to final archive...")
        for key in tqdm(expected_keys, desc="Copying to final archive", leave=False, **TQDM_DEFAULTS):  # type: ignore[call-overload]
            log.debug(f"Copying {key}...")
            with s3_errors(f"Copy of {key} to s3://{run_state.final_bucket}", UploadError):
                _copy_object(
                    run_state.final_s3,
                    run_state.interrogation_bucket,
                    run_state.final_bucket,
                    key,
                    run_state.final_part_size,
                )

        log.info("Copy complete. Removing files from interrogation bucket...")
        self._delete_staged_copies(run_state, expected_keys)
        log.info("Finished removing temporary files from interrogation bucket.")

    @staticmethod
    def _delete_staged_copies(run_state: SubmissionRunState, keys: list[str]) -> None:
        """Delete staged copies from the interrogation bucket, and log each delete that fails instead of raising it.

        A staged copy left behind only costs storage, so a failed delete never decides the outcome of a run.
        """
        for key in tqdm(keys, desc="Cleaning staging area", leave=False, **TQDM_DEFAULTS):  # type: ignore[call-overload]
            try:
                run_state.interrogation_s3.delete_object(Bucket=run_state.interrogation_bucket, Key=key)
            except Exception as e:
                log.warning(f"Failed to delete {key} from the interrogation bucket: {e}")

    def _handle_interrogation_failure(self, run_state: SubmissionRunState) -> None:
        if self.config.archives.interrogation.keep_failed:
            log.info("Interrogation config keep_failed is True. Leaving failed files in interrogation bucket.")
            return

        log.warning("Cleaning up interrogation bucket due to failure...")
        self._delete_staged_copies(run_state, self._get_expected_keys(run_state))

    @staticmethod
    def _raise_on_file_errors(run_state: SubmissionRunState, message: str) -> None:
        """Log the file errors of a pass and raise them as one error.

        :param run_state: The submission's run state.
        :param message: What failed, such as ``"Processing failed"``.
        :raises FilesFailedError: If the pass recorded a file error. Its ``__cause__`` is the decisive one.
        """
        metadata_order = {str(f.file_path): i for i, f in enumerate(run_state.submission_metadata.files.values())}
        file_errors = sorted(run_state.context.errors, key=lambda file_error: metadata_order[file_error.file])
        if not file_errors:
            return
        log.error(f"{message}: {'; '.join(f'{file_error.file}: {file_error.error}' for file_error in file_errors)}")
        error = FilesFailedError(message, file_errors)
        raise error from error.decisive.error

    def run(self, submission_metadata: SubmissionMetadata) -> None:
        """
        Execute the processing pipeline for a single submission.

        High-level view:
        1. Determine the target archive from the consent status at the submission date.
           If it already holds the submission's metadata, stop.
        2. Guess whether the submission will be selected for detailed QC.
        3. Spawn threads to process files (Download -> Decrypt -> Validate -> Encrypt -> Archive).
           If the guess is yes, also write the decrypted data to local QC storage.
        4. Mark basic QC as passed in the database.
        5. Determine detailed QC eligibility and update the database accordingly.
        6. If detailed QC is needed, fetch the files that step 3 did not write to local storage.
           If it is not needed, delete the files that step 3 wrote.
        7. Stage redacted metadata.
        8. Copy files from interrogation bucket to target archive bucket, then clean up.
        9. Optionally clean the inbox.

        See ``packages/grzctl/docs/process.md`` for the flow, including parallel runs.

        The interrogation bucket serves as a staging area: files are first uploaded there
        during processing, then copied to the final archive on success. If processing fails,
        files are either cleaned up or left in the interrogation bucket depending on the
        ``keep_failed`` configuration. A lifecycle rule on the interrogation bucket is
        recommended to automatically clean up orphaned files from incomplete transfers.

        :param submission_metadata: The parsed metadata object containing donor and file information.
        :raises FilesFailedError: If processing a file fails, in the main pass or in the detailed QC pass.
            Its ``__cause__`` is the decisive file error.
        """
        submission_run = self._new_run_state(submission_metadata)
        if self._is_archived(submission_run):
            return
        db = SubmissionDb(self.config.db.database_url, self.config.db.signing_author)
        prefetch = self._predict_qc(db, submission_run.submission_id)
        if prefetch:
            log.info(
                f"Submission {submission_run.submission_id} is likely to be selected for detailed QC, "
                "writing its decrypted files to QC storage during processing."
            )

        selected_for_qc = False
        try:
            self._pipeline_executor.process_files(submission_run, stage=True, write_local=prefetch)

            self._raise_on_file_errors(submission_run, "Processing failed")

            # validation passed, so mark basic QC as passed in the database.
            try:
                db.modify_submission(submission_run.submission_id, "basic_qc_passed", True)
            except DuplicateInitialSubmissionError as e:
                # another initial submission of this case passed basic QC while this one was processed
                log.warning(
                    f"Submission '{submission_run.submission_id}' data validated, but {e} "
                    "Failing basic QC for this submission."
                )
                db.modify_submission(submission_run.submission_id, "basic_qc_passed", False)
                raise

            # determine whether to perform detailed QC (now that basic QC is marked as passed).
            selected_for_qc = self._determine_qc_flag(db, submission_run.submission_id)
            if not selected_for_qc and prefetch:
                self._discard_qc_prefetch(submission_run)

            if selected_for_qc:
                log.info(f"Running detailed QC pass for {submission_run.submission_id}...")
                detailed_qc = self.config.detailed_qc

                # writes only the local copies that the main pass did not write
                self._pipeline_executor.process_files(submission_run, stage=False, write_local=True)
                self._raise_on_file_errors(submission_run, "Writing the files for detailed QC failed")

                # write metadata to local storage for the QC workflow
                submission_basepath = Path(detailed_qc.local_storage) / submission_run.submission_id
                metadata_dir = submission_basepath / "metadata"
                _ensure_directory_mode(Path(detailed_qc.local_storage), metadata_dir, QC_DIRECTORY_MODE)
                metadata_file = metadata_dir / "metadata.json"
                metadata_file.write_text(json.dumps(submission_metadata.content.get_raw_dict(), indent=2))
                log.info(f"Wrote submission metadata to {metadata_file}")

                if detailed_qc.auto_run:
                    output_basepath = submission_basepath / "qc"
                    shell_command = detailed_qc.shell_command.format(
                        submission_basepath=str(submission_basepath),
                        output_basepath=str(output_basepath),
                        submission_id=submission_run.submission_id,
                    )
                    log.info(f"Running detailed QC workflow: {shell_command}")
                    try:
                        subprocess.run(shell_command, shell=True, check=True)  # noqa: S602
                    except subprocess.CalledProcessError as e:
                        raise DetailedQCError(f"Detailed QC workflow failed: {e}") from e
                    log.info(f"Detailed QC workflow completed for {submission_run.submission_id}.")

            self._upload_final_metadata(submission_metadata, submission_run)
            self._upload_redacted_logs(submission_metadata, submission_run)

            # Copy files from interrogation bucket to final archive.
            self._commit_to_archive(submission_run)

            log.info(f"Submission {submission_run.submission_id} processed successfully.")
        except Exception:
            # Validation failed, or a step after the upload to the interrogation bucket did
            # (copy to final archive, metadata upload, QC workflow). Either way, clean up the
            # staged files so they don't linger.
            self._handle_interrogation_failure(submission_run)
            # Decrypted files stay only for a submission that is selected for detailed QC.
            if not selected_for_qc and prefetch:
                self._discard_qc_prefetch(submission_run)
            raise
