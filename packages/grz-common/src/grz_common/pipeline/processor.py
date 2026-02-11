import logging
from concurrent.futures import ThreadPoolExecutor
from contextlib import ExitStack
from datetime import date
from typing import Any

from grz_common.constants import TQDM_DEFAULTS
from grz_common.workers.submission import SubmissionMetadata
from tqdm.auto import tqdm

from ..models.s3 import S3Options
from ..transfer import init_s3_client
from .context import ConsistencyValidator, SubmissionContext
from .crypt4gh import Crypt4GHDecryptor, Crypt4GHEncryptor
from .s3 import S3Downloader, S3MultipartUploader
from .tqdm import TqdmObserver
from .validators import ValidatorObserver

log = logging.getLogger(__name__)


class SubmissionProcessor:
    """
    Orchestrates the streaming submission pipeline:
    Inbox -> Decrypt -> Validate -> Re-Encrypt -> Archive
    """

    def __init__(
        self,
        configuration: Any,
        source_s3_options: S3Options,
        keys: dict[str, bytes],
        redact_patterns: list[tuple[str, str]] | None = None,
    ):
        self.config = configuration
        self.source_s3_options = source_s3_options
        self.keys = keys
        self.redact_patterns = redact_patterns or []

        self.context = SubmissionContext()
        self.consistency = ConsistencyValidator(self.context)

        self.source_s3 = init_s3_client(source_s3_options)
        self._target_s3_map = {}

    def _get_target_s3(self, options: S3Options):
        """Get or create an S3 client for the target options."""
        key = (options.endpoint_url, options.bucket)
        if key not in self._target_s3_map:
            self._target_s3_map[key] = init_s3_client(options)
        return self._target_s3_map[key]

    def run(self, submission_metadata: SubmissionMetadata, threads: int = 1):
        """
        Run the pipeline.
        """
        submission_id = submission_metadata.content.submission_id
        is_consented = submission_metadata.content.consents_to_research(date.today())
        target_s3_options = self.config.consented_archive_s3 if is_consented else self.config.non_consented_archive_s3
        target_s3 = self._get_target_s3(target_s3_options)

        target_public_key = self.keys["consented_public"] if is_consented else self.keys["non_consented_public"]

        log.info(f"Consent Status: {'Consented' if is_consented else 'Non-Consented'}")
        log.info(f"Target Archive: {target_s3_options.bucket}")

        files_map = submission_metadata.files
        total_bytes = sum(f.file_size_in_bytes for f in files_map.values())

        partner_map = {} #self._build_partner_map(submission_metadata.content)

        log.info(f"Processing {len(files_map)} files ({total_bytes / (1024**3):.2f} GB)...")

        with (
            tqdm(total=total_bytes, desc="Processing", **TQDM_DEFAULTS) as pbar_global,
            ThreadPoolExecutor(max_workers=threads) as pool,
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
                        pbar=pbar_global,
                    )
                )
            for f in futures:
                f.result()

        # 6. Final Consistency Check
        if self.context.has_errors:
            raise RuntimeError("Submission failed consistency checks or validation.")

        # 7. Upload Metadata Completion Marker (Encrypted)
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
        pbar: Any,
    ):
        if self.context.has_errors:
            return

        try:
            self._run_pipeline(
                file_meta, submission_id, source_bucket, target_s3, target_bucket, target_public_key, pbar
            )

            # Fail-Fast Consistency
            partner = partner_map.get(file_meta.file_path)
            if partner and not self.consistency.check_pair(file_meta.file_path, partner):
                raise RuntimeError(f"Consistency Check Failed: {file_meta.file_path}")

        except Exception as e:
            log.error(f"Failed {file_meta.file_path}: {e}")
            self.context.add_error(str(e))

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
            # Pipeline Construction
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
                part_size=self.config.chunk_size,
                max_threads=self.config.max_concurrent_uploads,
            )

            uploader.upload(monitored)

            # Record Validation Metrics
            self.context.record_stats(file_meta.file_path, validator.metrics)

    def _upload_final_metadata(
        self, submission_metadata: SubmissionMetadata, sub_id: str, s3_client, bucket: str, public_key: bytes
    ):
        """Redacts and uploads the final metadata."""
        dest_key = f"{sub_id}/metadata/metadata.json"
        redacted_metadata = submission_metadata  # TODO: use redaction
        s3_client.put_object(Bucket=bucket, Key=dest_key, Body=redacted_metadata)
