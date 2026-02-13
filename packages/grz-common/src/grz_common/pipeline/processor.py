import logging
import shutil
from concurrent.futures import ThreadPoolExecutor
from contextlib import ExitStack
from datetime import date
from io import BytesIO
from pathlib import Path
from typing import Any

from grz_common.constants import TQDM_DEFAULTS
from grz_common.workers.submission import SubmissionMetadata
from grzctl.models.config import ProcessConfig
from tqdm.auto import tqdm

from ..models.s3 import S3Options
from ..progress import FileProgressLogger, ProcessingState
from ..transfer import init_s3_client
from .components import (
    Crypt4GHDecryptor,
    Crypt4GHEncryptor,
    MeasuringStream,
    MetricsRegistry,
    S3Downloader,
    S3MultipartUploader,
    TqdmObserver,
    ValidatorObserver,
    calculate_s3_part_size,
)
from .context import ConsistencyValidator, SubmissionContext

log = logging.getLogger(__name__)


class SubmissionProcessor:
    """
    Orchestrates the streaming submission pipeline:
    Inbox -> Decrypt -> Validate -> Re-Encrypt -> Archive
    """

    def __init__(
        self,
        configuration: ProcessConfig,
        source_s3_options: S3Options,
        keys: dict[str, bytes],
        status_file_path: Path,
        redact_patterns: list[tuple[str, str]] | None = None,
        max_concurrent_uploads: int = 1,
        threads: int = 1,
        enable_metrics: bool = False,
    ):
        self.config = configuration
        self.source_s3_options = source_s3_options
        self.keys = keys
        self.redact_patterns = redact_patterns or []
        self.enable_metrics = enable_metrics

        self.context = SubmissionContext()
        self.consistency = ConsistencyValidator(self.context)
        self.progress_logger = FileProgressLogger[ProcessingState](status_file_path)

        self.threads = threads
        self.max_concurrent_uploads = max_concurrent_uploads

        self.pool_size = max(10, threads * (1 + max_concurrent_uploads) + 1)
        log.debug(f"Configuring S3 client pool size: {self.pool_size}")

        self.source_s3 = init_s3_client(s3_options=source_s3_options, max_pool_connections=self.pool_size)
        self._target_s3_map = {}

        self.submission_id: str | None = None
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

    def _measure(self, stack: ExitStack, stream: Any, name: str, registry: MetricsRegistry | None) -> Any:
        """Conditionally wrap a stream with measuring instrumentation."""
        if self.enable_metrics and registry:
            return stack.enter_context(MeasuringStream(stream, name, registry))
        return stream

    def run(self, submission_metadata: SubmissionMetadata):
        self.submission_id = submission_metadata.content.submission_id
        is_consented = submission_metadata.content.consents_to_research(date.today())

        target_s3_options = self.config.consented_archive_s3 if is_consented else self.config.non_consented_archive_s3
        self.target_bucket = target_s3_options.bucket
        self.target_s3 = self._get_target_s3(target_s3_options)
        self.target_public_key = self.keys["consented_public"] if is_consented else self.keys["non_consented_public"]

        log.info(f"Consent Status: {'Consented' if is_consented else 'Non-Consented'}")
        log.info(f"Target Archive: {self.target_bucket}")

        files_map = submission_metadata.files
        total_bytes = sum(f.file_size_in_bytes for f in files_map.values())

        # TODO: Implement helper in metadata model
        # for r1, r2 in submission_metadata.get_paired_fastqs():
        #     self.partner_map[r1.file_path] = r2.file_path
        #     self.partner_map[r2.file_path] = r1.file_path

        log.info(f"Processing {len(files_map)} files ({total_bytes / (1024**3):.2f} GB)...")

        with (
            tqdm(total=total_bytes, desc="Total     ", position=0, **TQDM_DEFAULTS) as pbar_global,
            ThreadPoolExecutor(max_workers=self.threads) as pool,
        ):
            futures = [
                pool.submit(self._process_one_file, file_meta=file_meta, pbar_global=pbar_global)
                for file_meta in files_map.values()
            ]
            for f in futures:
                f.result()

        if self.context.has_errors:
            raise RuntimeError("Submission failed consistency checks or validation.")

        self._upload_final_metadata(submission_metadata)
        log.info(f"Submission {self.submission_id} processed successfully.")

    def _process_one_file(self, file_meta: Any, pbar_global: Any):
        src_key = f"{self.submission_id}/files/{file_meta.file_path}.c4gh"
        file_path_str = str(file_meta.file_path)
        file_name = file_path_str.split("/")[-1]

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
            with tqdm(
                total=file_meta.file_size_in_bytes,
                desc="Processing",
                postfix={"file": file_name},
                leave=False,
                **TQDM_DEFAULTS,
            ) as pbar_local:
                self._run_pipeline(file_meta, pbar=[pbar_global, pbar_local])

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

    def _run_pipeline(self, file_meta: Any, pbar: Any):
        src_key = f"{self.submission_id}/files/{file_meta.file_path}.c4gh"
        dest_key = f"{self.submission_id}/files/{file_meta.encrypted_file_path()}"

        part_size = calculate_s3_part_size(file_meta.file_size_in_bytes, None)
        metrics = MetricsRegistry() if self.enable_metrics else None

        with ExitStack() as stack:
            source = stack.enter_context(S3Downloader(self.source_s3, self.source_s3_options.bucket, src_key))
            source = self._measure(stack, source, "1_Source", metrics)

            decrypted = stack.enter_context(Crypt4GHDecryptor(source, self.keys["private"]))
            decrypted = self._measure(stack, decrypted, "2_Decrypt", metrics)

            validator = stack.enter_context(
                ValidatorObserver(decrypted, file_type=file_meta.file_type, expected_checksum=file_meta.file_checksum)
            )
            validator_stream = self._measure(stack, validator, "3_Validate", metrics)

            encrypted = stack.enter_context(Crypt4GHEncryptor(validator_stream, self.target_public_key))
            encrypted = self._measure(stack, encrypted, "4_Encrypt", metrics)

            monitored = stack.enter_context(TqdmObserver(encrypted, pbar))
            monitored = self._measure(stack, monitored, "5_Monitor", metrics)

            uploader = stack.enter_context(
                S3MultipartUploader(
                    self.target_s3,
                    self.target_bucket,
                    dest_key,
                    part_size=part_size,
                    max_threads=self.max_concurrent_uploads,
                )
            )

            shutil.copyfileobj(monitored, uploader, length=part_size)

            validator.verify()

            if metrics:
                log.info(f"Performance for {file_meta.file_path}: {metrics.report()}")

            self.context.record_stats(file_meta.file_path, validator.metrics)

    def _upload_final_metadata(self, submission_metadata: SubmissionMetadata):
        """Redacts and uploads the final metadata."""
        dest_key = f"{self.submission_id}/metadata/metadata.json"
        # TODO: Implement proper redaction in SubmissionMetadata model
        redacted_metadata = submission_metadata.content.to_redacted_dict()

        self.target_s3.put_object(
            Bucket=self.target_bucket, Key=dest_key, Body=BytesIO(str(redacted_metadata).encode("utf-8"))
        )
