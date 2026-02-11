import logging
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
        max_concurrent_uploads: int = 4,
        threads: int = 1,
    ):
        self.config = configuration
        self.source_s3_options = source_s3_options
        self.keys = keys
        self.redact_patterns = redact_patterns or []

        self.context = SubmissionContext()
        self.consistency = ConsistencyValidator(self.context)
        self.progress_logger = FileProgressLogger[ProcessingState](status_file_path)

        self.threads = threads
        self.max_concurrent_uploads = max_concurrent_uploads
        # Calculate required connection pool size to prevent exhaustion
        # Each thread needs: 1 connection for source (download) + N connections for target (multipart upload)
        # We add a small buffer (+5) for metadata operations
        self.pool_size = max(10, threads * (1 + max_concurrent_uploads) + 5)
        log.debug(f"Configuring S3 client pool size: {self.pool_size}")

        self.source_s3 = init_s3_client(s3_options=source_s3_options, max_pool_connections=self.pool_size)
        self._target_s3_map = {}

    def _get_target_s3(self, options: S3Options):
        """Get or create an S3 client for the target options."""
        key = (options.endpoint_url, options.bucket)
        if key not in self._target_s3_map:
            self._target_s3_map[key] = init_s3_client(s3_options=options, max_pool_connections=self.pool_size)
        return self._target_s3_map[key]

    def run(self, submission_metadata: SubmissionMetadata):
        submission_id = submission_metadata.content.submission_id
        is_consented = submission_metadata.content.consents_to_research(date.today())
        target_s3_options = self.config.consented_archive_s3 if is_consented else self.config.non_consented_archive_s3
        target_s3 = self._get_target_s3(target_s3_options)

        target_public_key = self.keys["consented_public"] if is_consented else self.keys["non_consented_public"]

        log.info(f"Consent Status: {'Consented' if is_consented else 'Non-Consented'}")
        log.info(f"Target Archive: {target_s3_options.bucket}")

        files_map = submission_metadata.files
        total_bytes = sum(f.file_size_in_bytes for f in files_map.values())

        partner_map = {}
        # for r1, r2 in submission_metadata.TODO_get_paired_fastqs():
        #     partner_map[r1.file_path] = r2.file_path
        #     partner_map[r2.file_path] = r1.file_path

        log.info(f"Processing {len(files_map)} files ({total_bytes / (1024**3):.2f} GB)...")

        with (
            tqdm(total=total_bytes, desc="Total     ", position=0, **TQDM_DEFAULTS) as pbar_global,
            ThreadPoolExecutor(max_workers=self.threads) as pool,
        ):
            futures = []
            for _rel_path, file_meta in files_map.items():
                futures.append(
                    pool.submit(
                        self._process_one_file,
                        file_meta=file_meta,
                        partner_map=partner_map,
                        submission_id=submission_id,
                        source_bucket=self.source_s3_options.bucket,
                        target_s3=target_s3,
                        target_bucket=target_s3_options.bucket,
                        target_public_key=target_public_key,
                        pbar_global=pbar_global,
                    )
                )
            for f in futures:
                f.result()

        if self.context.has_errors:
            raise RuntimeError("Submission failed consistency checks or validation.")

        # upload metadata.json as completion marker
        self._upload_final_metadata(
            submission_metadata, submission_id, target_s3, target_s3_options.bucket, target_public_key
        )

        log.info(f"Submission {submission_id} processed successfully.")

    def _process_one_file(  # noqa: PLR0913
        self,
        file_meta: Any,
        partner_map: dict[str, str],
        submission_id: str,
        source_bucket: str,
        target_s3: Any,
        target_bucket: str,
        target_public_key: bytes,
        pbar_global: Any,
    ):
        src_key = f"{submission_id}/files/{file_meta.file_path}.c4gh"
        file_path_str = str(file_meta.file_path)
        file_name = file_path_str.split("/")[-1]

        try:
            head = self.source_s3.head_object(Bucket=source_bucket, Key=src_key)
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
                self._run_pipeline(
                    file_meta,
                    submission_id,
                    source_bucket,
                    target_s3,
                    target_bucket,
                    target_public_key,
                    pbar=[pbar_global, pbar_local],
                )

            partner = partner_map.get(file_meta.file_path)
            if partner and not self.consistency.check_pair(file_meta.file_path, partner):
                raise RuntimeError(f"Consistency Check Failed: {file_meta.file_path}")

            self.progress_logger.set_state(
                file_path_str, file_meta, {"processing_successful": True, "errors": []}, size=s3_size, mtime=s3_mtime
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

    def _run_pipeline(  # noqa: PLR0913
        self,
        file_meta: Any,
        sub_id: str,
        src_bucket: str,
        target_s3: Any,
        target_bucket: str,
        target_pub_key: bytes,
        pbar: Any,
    ):
        src_key = f"{sub_id}/files/{file_meta.file_path}.c4gh"
        dest_key = f"{sub_id}/files/{file_meta.encrypted_file_path()}"

        with ExitStack() as stack:
            source = stack.enter_context(S3Downloader(self.source_s3, src_bucket, src_key))
            decrypted = stack.enter_context(Crypt4GHDecryptor(source, self.keys["private"]))

            validator = stack.enter_context(
                ValidatorObserver(decrypted, file_type=file_meta.file_type, expected_checksum=file_meta.file_checksum)
            )

            encrypted = stack.enter_context(Crypt4GHEncryptor(validator, target_pub_key))
            monitored = stack.enter_context(TqdmObserver(encrypted, pbar))

            uploader = S3MultipartUploader(
                target_s3,
                target_bucket,
                dest_key,
                part_size=calculate_s3_part_size(file_meta.file_size_in_bytes, None),
                max_threads=self.max_concurrent_uploads,
            )

            uploader.upload(monitored)
            validator.verify()

            self.context.record_stats(file_meta.file_path, validator.metrics)

    def _upload_final_metadata(
        self, submission_metadata: SubmissionMetadata, sub_id: str, s3_client, bucket: str, public_key: bytes
    ):
        """Redacts and uploads the final metadata."""
        dest_key = f"{sub_id}/metadata/metadata.json"
        redacted_metadata = submission_metadata.content.to_redacted_dict()  # TODO: use redaction
        s3_client.put_object(Bucket=bucket, Key=dest_key, Body=BytesIO(str(redacted_metadata).encode("utf-8")))
