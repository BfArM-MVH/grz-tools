import logging
from concurrent.futures import Future, ThreadPoolExecutor
from contextlib import ExitStack
from datetime import date
from io import BytesIO
from pathlib import Path
from typing import TYPE_CHECKING, Any

from crypt4gh.keys import get_public_key
from grz_common.constants import TQDM_DEFAULTS
from grz_common.workers.submission import SubmissionMetadata
from grz_db.models.submission import SubmissionDb
from grz_pydantic_models.submission.metadata import File, FileType
from grz_pydantic_models.submission.thresholds import Thresholds
from grzctl.models.config import InboxTarget, ProcessConfig
from pydantic import AnyHttpUrl
from tqdm.auto import tqdm

from ..models.s3 import S3Options
from ..progress import FileProgressLogger, ProcessingState
from ..transfer import init_s3_client
from ..utils.crypt import Crypt4GH
from .components import ObserverWithMetrics, Tee, TqdmObserver
from .components.crypt4gh import Crypt4GHDecryptor, Crypt4GHEncryptor
from .components.perf import MeasuringWriteStream, MeasuringReadStream, MetricsRegistry
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
        max_concurrent_uploads: int = 1,
        threads: int = 1,
        enable_metrics: bool = True,
        background_tee: bool = False,
    ):
        self.config = configuration
        self.source_s3_options = inbox.s3

        log.debug("Loading crypt4gh keys...")
        self.private_key = Crypt4GH.retrieve_private_key(
            inbox.private_key_path, passphrase=inbox.private_key_passphrase
        )
        self.consented_pub_key = get_public_key(configuration.archives.consented.public_key_path)
        self.non_consented_pub_key = get_public_key(configuration.archives.non_consented.public_key_path)

        self.redact_patterns = redact_patterns or []
        self.enable_metrics = enable_metrics
        self.background_tee = background_tee

        self.context = SubmissionContext()
        self.consistency = ConsistencyValidator(self.context)
        self.progress_logger = FileProgressLogger[ProcessingState](status_file_path)

        self.threads = threads
        self.max_concurrent_uploads = max_concurrent_uploads

        self.pool_size = max(10, threads * (1 + max_concurrent_uploads) + 1)
        log.debug(f"Configuring S3 client pool size: {self.pool_size}")

        self.source_s3 = init_s3_client(s3_options=self.source_s3_options, max_pool_connections=self.pool_size)
        self._target_s3_map: dict[tuple[AnyHttpUrl, str], S3Client] = {}

        self.submission_id: str | None = None
        self.should_qc: bool = False
        self.target_s3: Any = None
        self.target_bucket: str | None = None
        self.target_public_key: bytes | None = None
        self.partner_map: dict[str, str] = {}

    def _get_target_s3(self, options: S3Options):
        """Get or create an S3 client for the target options."""
        key = (options.endpoint_url, options.bucket)
        if key not in self._target_s3_map:
            self._target_s3_map[key] = init_s3_client(s3_options=options, max_pool_connections=self.pool_size)
        return self._target_s3_map[key]

    def _measure(self, stream: Any, name: str, registry: MetricsRegistry | None, is_observer: bool = False) -> Any:
        """Conditionally wrap a stream with measuring instrumentation."""
        if not self.enable_metrics or not registry:
            return stream
        if is_observer:
            wrapper = MeasuringWriteStream(name, registry)
            wrapper.set_sink(stream)
            return wrapper

        return MeasuringReadStream(stream, name, registry)

    def run(self, submission_metadata: SubmissionMetadata):
        self.submission_id = submission_metadata.content.submission_id
        is_consented = submission_metadata.content.consents_to_research(date.today())
        db = SubmissionDb(self.config.db.database_url, self.config.db.author)  # type: ignore[arg-type]
        target_percentage = self.config.detailed_qc.target_percentage
        self.should_qc = (
            db.should_qc(self.submission_id, target_percentage, self.config.detailed_qc.salt)
            if target_percentage > 0.0
            else False
        )

        target_archive = self.config.archives.consented if is_consented else self.config.archives.non_consented
        self.target_bucket = target_archive.s3.bucket
        self.target_s3 = self._get_target_s3(target_archive.s3)
        self.target_public_key = self.consented_pub_key if is_consented else self.non_consented_pub_key

        log.info(f"Consent Status: {'Consented' if is_consented else 'Non-Consented'}")
        log.info(f"Target Archive: {self.target_bucket}")

        files_map = submission_metadata.files
        total_bytes = sum(f.file_size_in_bytes for f in files_map.values())

        self.partner_map.clear()
        for _donor, _lab_datum, pairs, _thresholds in submission_metadata.iter_paired_end_fastqs():
            for fq1, fq2 in pairs:
                key1 = fq1.file_path
                key2 = fq2.file_path
                self.partner_map[key1] = key2
                self.partner_map[key2] = key1

        thresholds: dict[str, Thresholds] = {}
        for _, _, files, t in submission_metadata.iter_single_end_fastqs():
            for f in files:
                thresholds[f.file_path] = t
        for _, _, pairs, t in submission_metadata.iter_paired_end_fastqs():
            for fq1, fq2 in pairs:
                thresholds[fq1.file_path] = t
                thresholds[fq2.file_path] = t

        log.info(f"Processing {len(files_map)} files ({total_bytes / (1024**3):.2f} GB)...")

        with (
            tqdm(total=total_bytes, desc="Total     ", position=0, **TQDM_DEFAULTS) as pbar_global,  # type: ignore[call-overload]
            ThreadPoolExecutor(max_workers=self.threads) as pool,
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
        log.info(f"Submission {self.submission_id} processed successfully.")
        # TODO: cleanup inbox if successful and requested
        self.should_qc = False
        self.submission_id = None

    def _process_one_file(self, file_meta: File, threshold: Thresholds | None, pbar_global: Any):
        src_key = f"{self.submission_id}/files/{file_meta.file_path}.c4gh"
        file_path_str = str(file_meta.file_path)
        file_name = file_path_str.rsplit("/", maxsplit=1)[-1]

        try:
            head = self.source_s3.head_object(Bucket=self.source_s3_options.bucket, Key=src_key)
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

            partner = self.partner_map.get(file_meta.file_path)
            if partner and not self.consistency.check_pair(file_meta.file_path, partner):
                raise RuntimeError(f"Consistency Check Failed: {file_meta.file_path}")

            self.progress_logger.set_state(
                file_path_str,
                file_meta,
                {"processing_successful": True, "errors": []},
                size=s3_size,
                mtime=s3_mtime,
            )

        except Exception as e:
            log.error(f"Failed {file_meta.file_path}: {e}")
            self.context.add_error(str(e))
            self.progress_logger.set_state(
                file_path_str,
                file_meta,
                {"processing_successful": False, "errors": [str(e)]},
                size=s3_size,
                mtime=s3_mtime,
            )

    def _run_pipeline(self, file_meta: File, threshold: Thresholds | None, pbar_global: Any, file_name: str):
        if not self.target_public_key:
            raise RuntimeError("Target public key not set.")
        if not self.target_bucket:
            raise RuntimeError("Target bucket not set.")
        if not self.submission_id:
            raise RuntimeError("Submission id not set.")

        src_key = f"{self.submission_id}/files/{file_meta.file_path}.c4gh"
        dest_key = f"{self.submission_id}/files/{file_meta.encrypted_file_path()}"

        part_size = calculate_s3_part_size(file_meta.file_size_in_bytes, None)
        metrics = MetricsRegistry() if self.enable_metrics else None

        checksum_validator = ChecksumValidator(expected_checksum=file_meta.file_checksum)
        validation_chain = self._measure(checksum_validator, "3a_Checksum", metrics, is_observer=True)

        format_validator: ObserverWithMetrics | None = None
        if file_meta.file_type == FileType.fastq:
            format_validator = FastqValidator(
                mean_read_length_threshold=threshold.mean_read_length if threshold else None
            )
        elif file_meta.file_type == FileType.bam:
            format_validator = BamValidator()

        if format_validator:
            wrapped_format = self._measure(format_validator, "3b_Format", metrics, is_observer=True)
            validation_chain = validation_chain | wrapped_format

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
            pipeline = S3Downloader(self.source_s3, self.source_s3_options.bucket, src_key)
            pipeline = self._measure(pipeline, "1_Source", metrics)

            pipeline = pipeline | Crypt4GHDecryptor(private_key=self.private_key)
            pipeline = self._measure(pipeline, "2_Decrypt", metrics)

            if self.should_qc:
                path = Path(self.config.detailed_qc.local_storage) / self.submission_id / "files" / file_meta.file_path
                path.parent.mkdir(parents=True, exist_ok=True)
                writer = stack.enter_context(open(path, "wb"))
                writer = self._measure(writer, "2b_Write", metrics, is_observer=True)
                pipeline = pipeline | Tee(writer, threaded=self.background_tee)

            pipeline = pipeline | Tee(validation_chain, threaded=self.background_tee)

            pipeline = pipeline | Crypt4GHEncryptor(recipient_pubkey=self.target_public_key)
            pipeline = self._measure(pipeline, "4_Encrypt", metrics)

            pipeline = pipeline | Tee(TqdmObserver([pbar_global, pbar_local]), threaded=self.background_tee)

            uploader = S3MultipartUploader(
                self.target_s3,
                self.target_bucket,
                dest_key,
                part_size=part_size,
                max_threads=self.max_concurrent_uploads,
            )

            pipeline >> uploader

        stats = checksum_validator.metrics
        if format_validator:
            stats.update(format_validator.metrics)

        if metrics:
            log.info(f"Performance for {file_meta.file_path}: {metrics.report()}")

        self.context.record_stats(file_meta.file_path, stats)

    def _upload_final_metadata(self, submission_metadata: SubmissionMetadata):
        """Redacts and uploads the final metadata."""
        dest_key = f"{self.submission_id}/metadata/metadata.json"
        redacted_metadata = submission_metadata.content.to_redacted_dict()

        self.target_s3.put_object(
            Bucket=self.target_bucket, Key=dest_key, Body=BytesIO(str(redacted_metadata).encode("utf-8"))
        )
