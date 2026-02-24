"""Command for processing a submission."""

import json
import logging
import re
import time
from pathlib import Path
from typing import Any

import click
import grz_common.cli as grzcli
from grz_common.models.s3 import S3Options
from grz_common.pipeline.processor import SubmissionProcessor
from grz_common.workers.download import S3BotoDownloadWorker
from grz_common.workers.submission import EncryptedSubmission, SubmissionMetadata
from grz_db.errors import DuplicateSubmissionError, DuplicateTanGError
from grz_db.models.submission import SubmissionStateEnum

from ..dbcontext import DbContext
from ..models.config import InboxTarget, ProcessConfig
from ..models.pruefbericht import PruefberichtModel
from .db.cli import get_submission_db_instance
from .pruefbericht import _generate_pruefbericht_from_metadata
from .pruefbericht import _try_submit as _try_submit_pruefbericht

log = logging.getLogger(__name__)


@click.command()
@grzcli.configuration
@grzcli.submission_id
@grzcli.output_dir
@grzcli.threads
@grzcli.update_db
@click.option(
    "--submit-pruefbericht/--no-submit-pruefbericht",
    default=False,
    help="Submit Prüfbericht to BfArM after successful processing.",
)
@click.option(
    "--save-pruefbericht",
    type=click.Path(),
    default=None,
    help="Save generated Prüfbericht to the specified path.",
)
@click.option(
    "--redact-logs/--no-redact-logs",
    default=True,
    help="Redact sensitive information from logs before archiving.",
)
@click.option(
    "--concurrent-uploads",
    type=int,
    default=4,
    help="Maximum concurrent part uploads per file.",
)
@click.option(
    "--inbox-bucket",
    default=None,
    help="Inbox bucket name to use. Required when a submitter has multiple inboxes configured.",
)
@click.option(
    "--clean-inbox/--no-clean-inbox",
    default=True,
    help="Clean submission from inbox bucket after successful processing.",
)
def process(  # noqa: PLR0913
    configuration: dict[str, Any],
    submission_id: str,
    output_dir: str,
    threads: int,
    update_db: bool,
    submit_pruefbericht: bool,
    save_pruefbericht: str | None,
    redact_logs: bool,
    concurrent_uploads: int,
    inbox_bucket: str | None = None,
    clean_inbox: bool = True,
    **kwargs,
):
    """
    Process a submission through the streaming pipeline.
    """
    config = ProcessConfig.model_validate(configuration)

    inbox = _resolve_inbox_target(config, submission_id, inbox_bucket)
    _, metadata_dir, log_dir, encrypted_files_dir = _setup_directories(output_dir)

    log.info(f"Starting streaming pipeline for submission: {submission_id}")

    # first, download metadata to understand the submission structure
    log.info("Downloading metadata...")
    download_worker = S3BotoDownloadWorker(
        inbox.s3, status_file_path=log_dir / "progress_download.cjson", threads=threads
    )
    download_worker.download_metadata(submission_id, metadata_dir, metadata_file_name="metadata.json")
    local_metadata_path = metadata_dir / "metadata.json"

    submission_metadata = SubmissionMetadata(local_metadata_path)

    # TODO see https://github.com/BfArM-MVH/grz-tools/pull/517
    redact_patterns = []
    if redact_logs:
        redact_patterns.append((re.escape(submission_metadata.content.submission.tan_g), "REDACTED_TAN"))
        redact_patterns.append(
            (re.escape(submission_metadata.content.submission.local_case_id), "REDACTED_LOCAL_CASE_ID")
        )

    # register submission in db if not yet registered
    if update_db and config.db:
        db_service = get_submission_db_instance(config.db.database_url)
        try:
            if not db_service.get_submission(submission_id):
                _db_submission = db_service.add_submission(submission_id)
        except (DuplicateSubmissionError, DuplicateTanGError) as e:
            raise click.Abort() from e
        except Exception as e:
            raise click.ClickException(f"Failed to add submission: {e}") from e

    status_file_path = log_dir / "progress_processing.cjson"

    processor = SubmissionProcessor(
        configuration=config,
        inbox=inbox,
        redact_patterns=redact_patterns,
        status_file_path=status_file_path,
        threads=threads,
        max_concurrent_uploads=concurrent_uploads,
        clean_inbox=clean_inbox,
        update_db=update_db,
    )

    with DbContext(
        configuration=configuration,
        submission_id=submission_id,
        start_state=SubmissionStateEnum.PROCESSING,
        end_state=SubmissionStateEnum.PROCESSED,
        enabled=update_db,
    ):
        processor.run(submission_metadata)

    # TODO: have a pruefbericht function that doesn't need an EncryptedSubmission instance?
    encrypted_submission = EncryptedSubmission(
        metadata_dir=metadata_dir,
        encrypted_files_dir=encrypted_files_dir,
        log_dir=log_dir,
    )

    _handle_pruefbericht(
        config=config,
        configuration=configuration,
        encrypted_submission=encrypted_submission,
        submission_id=submission_id,
        log_dir=log_dir,
        submit_pruefbericht=submit_pruefbericht,
        save_pruefbericht=save_pruefbericht,
        update_db=update_db,
    )


def _resolve_inbox_target(config: ProcessConfig, submission_id: str, requested_bucket: str | None) -> InboxTarget:
    le_id = submission_id.split("_", maxsplit=1)[0]

    if le_id not in config.s3.inboxes:
        available = ", ".join(config.s3.inboxes.keys())
        raise click.ClickException(f"Submitter ID '{le_id}' not found in configuration. Available: {available}")

    submitter_inboxes = config.s3.inboxes[le_id]

    if requested_bucket:
        if requested_bucket not in submitter_inboxes:
            raise click.ClickException(f"Inbox bucket '{requested_bucket}' not found for '{le_id}'.")
        bucket_name = requested_bucket
        inbox_config = submitter_inboxes[requested_bucket]
    elif len(submitter_inboxes) == 1:
        bucket_name, inbox_config = next(iter(submitter_inboxes.items()))
    else:
        available_buckets = ", ".join(submitter_inboxes.keys())
        raise click.ClickException(
            f"Multiple inboxes found for '{le_id}' ({available_buckets}). Please specify --inbox-bucket."
        )

    s3_options = S3Options(bucket=bucket_name, **inbox_config.model_dump())
    return InboxTarget(s3=s3_options, private_key_path=inbox_config.private_key_path)


def _setup_directories(output_dir: str) -> tuple[Path, Path, Path, Path]:
    """Create and return required directories."""
    base_dir = Path(output_dir)
    metadata_dir = base_dir / "metadata"
    log_dir = base_dir / "logs"
    encrypted_files_dir = base_dir / "encrypted_files"

    for d in [base_dir, metadata_dir, log_dir, encrypted_files_dir]:
        d.mkdir(mode=0o770, parents=True, exist_ok=True)

    return base_dir, metadata_dir, log_dir, encrypted_files_dir


def _prepare_redact_patterns(encrypted_submission: EncryptedSubmission) -> list[tuple[str, str]]:
    """
    Prepare redaction patterns for tanG and localCaseId.

    :param encrypted_submission: The encrypted submission
    :returns: List of (pattern, replacement) tuples
    """
    patterns = encrypted_submission.metadata.content.create_redaction_patterns()
    if patterns:
        log.info(f"Log redaction enabled with {len(patterns)} pattern(s)")
    return patterns


def _handle_pruefbericht(  # noqa: C901, PLR0913
    config: ProcessConfig,
    configuration: dict[str, Any],
    encrypted_submission: EncryptedSubmission,
    submission_id: str,
    log_dir: Path,
    submit_pruefbericht: bool,
    save_pruefbericht: str | None,
    update_db: bool,
    max_retries: int = 10,
) -> None:
    """
    Generate and optionally submit Prüfbericht to BfArM.

    This implements steps 1.8 of the SOP: Prüfbericht generation and submission.
    """
    metadata = encrypted_submission.metadata.content

    # generate the Prüfbericht
    log.info("Generating Prüfbericht...")
    try:
        pruefbericht = _generate_pruefbericht_from_metadata(metadata, failed=False)
        log.info("Prüfbericht generated successfully")
    except Exception as e:
        log.error(f"Failed to generate Prüfbericht: {e}")
        if submit_pruefbericht:
            raise
        return

    # save Prüfbericht with redacted TAN if requested
    if save_pruefbericht:
        save_path = Path(save_pruefbericht)
        redacted_data = pruefbericht.model_dump(by_alias=True, mode="json")
        redacted_data["SubmittedCase"]["tan"] = "<REDACTED>"
        with open(save_path, "w") as f:
            json.dump(redacted_data, f, indent=2)
        log.info(f"Saved Prüfbericht (with redacted TAN) to: {save_path}")

    # also save a copy to the logs directory (with redacted TAN)
    pruefbericht_log_path = log_dir / "pruefbericht.json"
    redacted_for_log = pruefbericht.model_dump(by_alias=True, mode="json")
    redacted_for_log["SubmittedCase"]["tan"] = "<REDACTED>"
    with open(pruefbericht_log_path, "w") as f:
        json.dump(redacted_for_log, f, indent=2)
    log.info(f"Saved Prüfbericht copy to logs: {pruefbericht_log_path}")

    # submit Prüfbericht if requested
    if submit_pruefbericht:
        # pruefbericht_config is guaranteed to exist as it's mandatory in ProcessConfig
        pruefbericht_config: PruefberichtModel = config.pruefbericht

        if (auth_url := pruefbericht_config.authorization_url) is None:
            raise ValueError("pruefbericht.authorization_url is required but not configured")
        if (client_id := pruefbericht_config.client_id) is None:
            raise ValueError("pruefbericht.client_id is required but not configured")
        if (client_secret := pruefbericht_config.client_secret) is None:
            raise ValueError("pruefbericht.client_secret is required but not configured")
        if (api_base_url := pruefbericht_config.api_base_url) is None:
            raise ValueError("pruefbericht.api_base_url is required but not configured")

        log.info("Submitting Prüfbericht to BfArM...")

        initial_delay = 30.0
        backoff_factor = 2.0

        with DbContext(
            configuration=configuration,
            submission_id=submission_id,
            start_state=SubmissionStateEnum.REPORTING,
            end_state=SubmissionStateEnum.REPORTED,
            enabled=update_db,
        ):
            for attempt in range(1, max_retries + 2):
                try:
                    _expiry, _token = _try_submit_pruefbericht(
                        pruefbericht=pruefbericht,
                        api_base_url=str(api_base_url),
                        auth_url=str(auth_url),
                        client_id=client_id,
                        client_secret=client_secret,
                        token="",
                    )
                    break
                except Exception as e:
                    if attempt > max_retries:
                        log.error(f"Prüfbericht submission failed after {max_retries} retries.")
                        raise e
                    wait_time = initial_delay * (backoff_factor ** (attempt - 1))

                    log.warning(
                        f"Prüfbericht submission attempt {attempt} failed with error: {e}. Retrying in {wait_time} seconds..."
                    )
                    time.sleep(wait_time)

        log.info("Prüfbericht submitted successfully!")
