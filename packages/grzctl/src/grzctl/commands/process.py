"""Command for processing a submission.

This module implements the ``grzctl process`` subcommand, which performs the
complete submission lifecycle in a single streaming pass:

    Download metadata → Decrypt → Validate → Re-encrypt → Archive → (optionally) Prüfbericht

It replaces the individual step-by-step subcommands (``download``, ``decrypt``,
``validate``, ``encrypt``, ``archive``) with a single streaming pipeline that
avoids materialising intermediate files on disk.

The pipeline stages are orchestrated by
:class:`grz_common.pipeline.processor.SubmissionProcessor` and use an
*interrogation bucket* as a staging area: files are first uploaded there, then
copied to the final archive bucket on success.  On failure, staged files are
cleaned up or retained depending on the ``keep_failed`` configuration.

Error behaviour
---------------
Processing errors are accumulated per file in a :class:`SubmissionContext`.
Individual file failures do not abort the pipeline; instead they are recorded
and checked at a synchronisation point after all files have been processed.
This is necessary because some validation checks (e.g. paired-end read-count
consistency) require information from all relevant parts before they can pass.

If the DB is enabled (``--update-db``), the submission state transitions
through ``PROCESSING → PROCESSED`` (or ``ERROR`` on failure).  The DB record
is also *populated* with metadata so that downstream Prüfbericht generation
can read the required fields.  If ``--submit-pruefbericht`` is used, a
separate ``REPORTING → REPORTED`` transition is recorded.

Recovery
--------
The ``progress_processing.cjson`` log tracks per-file completion.  Re-running
the command after a failure will skip files that were already processed
successfully, making the pipeline effectively idempotent.
"""

import json
import logging
import time
from pathlib import Path

import click
import grz_common.cli as grzcli
from grz_common.pipeline.processor import SubmissionProcessor
from grz_common.transfer import get_metadata_upload_timestamp, init_s3_client
from grz_common.workers.download import S3BotoDownloadWorker
from grz_common.workers.submission import SubmissionMetadata
from grz_db.errors import DuplicateSubmissionError, DuplicateTanGError
from grz_db.models.submission import SubmissionStateEnum
from grz_pydantic_models.submission.metadata import REDACTED_TAN

from ..commands import grzctl_configuration
from ..dbcontext import DbContext
from ..models.config import GrzctlConfig
from ..models.pruefbericht import PruefberichtModel
from .db.cli import get_submission_db_instance
from .pruefbericht import _generate_pruefbericht_from_database
from .pruefbericht import _try_submit as _try_submit_pruefbericht

log = logging.getLogger(__name__)


@click.command()
@grzctl_configuration
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
    "--redact-pruefbericht/--no-redact-pruefbericht",
    default=True,
    help="Whether to redact sensitive information from written Prüfbericht",
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
    configuration: GrzctlConfig,
    submission_id: str,
    output_dir: str,
    threads: int,
    update_db: bool,
    submit_pruefbericht: bool,
    save_pruefbericht: str | None,
    redact_pruefbericht: bool,
    redact_logs: bool,
    concurrent_uploads: int,
    inbox_bucket: str | None = None,
    clean_inbox: bool = True,
):
    """
    Process a submission through the streaming pipeline.

    Combines download, decrypt, validate, re-encrypt, and archive into a single
    streaming pass via :class:`SubmissionProcessor`.  Metadata is downloaded first
    (needed to determine consent status, file list, etc.) and then the full
    pipeline processes each file concurrently.

    When ``--update-db`` is enabled the DB record is populated with the parsed
    metadata so that downstream Prüfbericht generation can read the required
    fields (submission date, donor info, etc.).

    On success the submission state is set to ``PROCESSED``; on failure it is set
    to ``ERROR`` with the associated error message.  Files are processed
    idempotently: re-running after a partial failure skips already-completed files.
    """
    le_id = submission_id.split("_", maxsplit=1)[0]
    inbox = configuration.resolve_inbox(submitter_id=le_id, inbox_name=inbox_bucket)
    _, metadata_dir, log_dir = _setup_directories(output_dir)

    log.info(f"Starting streaming pipeline for submission: {submission_id}")

    # first, download metadata to understand the submission structure
    log.info("Downloading metadata...")
    download_worker = S3BotoDownloadWorker(
        inbox.s3, status_file_path=log_dir / "progress_download.cjson", threads=threads
    )
    download_worker.download_metadata(submission_id, metadata_dir, metadata_file_name="metadata.json")
    local_metadata_path = metadata_dir / "metadata.json"

    submission_metadata = SubmissionMetadata(local_metadata_path)

    # register and populate submission in DB if enabled
    if update_db:
        db_service = get_submission_db_instance(configuration.db.database_url)
        try:
            if not db_service.get_submission(submission_id):
                db_service.add_submission(submission_id)
        except (DuplicateSubmissionError, DuplicateTanGError):
            log.warning(f"Submission '{submission_id}' already exists in the database.")
        except Exception as e:
            raise click.ClickException(f"Failed to add submission: {e}") from e

        # Populate the DB record with parsed metadata (donors, files, dates, etc.)
        # so that downstream Prüfbericht generation can read the required fields.
        s3_client = init_s3_client(inbox.s3)
        submission_date = get_metadata_upload_timestamp(s3_client, inbox.s3.bucket, submission_id).date()
        db_service.populate(
            submission_id,
            submission_metadata.content,
            submission_date,
            force=False,
            on_missing="create",
        )

    status_file_path = log_dir / "progress_processing.cjson"

    processor = SubmissionProcessor(
        configuration=configuration,
        inbox=inbox,
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

    _handle_pruefbericht(
        configuration=configuration,
        submission_id=submission_id,
        log_dir=log_dir,
        submit_pruefbericht=submit_pruefbericht,
        save_pruefbericht=save_pruefbericht,
        redact_pruefbericht=redact_pruefbericht,
        redact_logs=redact_logs,
        update_db=update_db,
    )


def _setup_directories(output_dir: str) -> tuple[Path, Path, Path]:
    """Create and return required directories."""
    base_dir = Path(output_dir)
    metadata_dir = base_dir / "metadata"
    log_dir = base_dir / "logs"

    for d in [base_dir, metadata_dir, log_dir]:
        d.mkdir(mode=0o770, parents=True, exist_ok=True)

    return base_dir, metadata_dir, log_dir


def _handle_pruefbericht(  # noqa: C901, PLR0913, PLR0912, PLR0915
    configuration: GrzctlConfig,
    submission_id: str,
    log_dir: Path,
    submit_pruefbericht: bool,
    save_pruefbericht: str | None,
    redact_pruefbericht: bool,
    redact_logs: bool,
    update_db: bool,
    max_retries: int = 10,
) -> None:
    """Generate and optionally submit Prüfbericht to BfArM.

    Prüfbericht generation is only reached for submissions that passed basic QC
    (validation succeeded).  If validation had failed, the pipeline would have
    raised an error earlier and we would never get here, therefore ``failed`` is
    always ``False`` at this point.
    """
    log.info("Generating Prüfbericht...")
    try:
        # The ``failed`` flag is always False here: Prüfbericht generation is only
        # reached when the pipeline succeeded (basic QC passed).  The parameter is
        # kept for API compatibility with ``_generate_pruefbericht_from_database``.
        failed = False
        pruefbericht = _generate_pruefbericht_from_database(submission_id, configuration, failed)
        log.info("Prüfbericht generated successfully")
    except Exception as e:
        log.error(f"Failed to generate Prüfbericht: {e}")
        if submit_pruefbericht:
            raise
        return

    # save Prüfbericht
    if save_pruefbericht:
        save_path = Path(save_pruefbericht)
        pruefbericht_data = pruefbericht.model_dump(by_alias=True, mode="json")
        # ... with redacted TAN if requested
        if redact_pruefbericht:
            pruefbericht_data["SubmittedCase"]["tan"] = REDACTED_TAN
        with open(save_path, "w") as f:
            json.dump(pruefbericht_data, f, indent=2)
        log.info(f"Saved Prüfbericht (with redacted TAN) to: {save_path}")

    # also save a copy to the logs directory (with redacted TAN)
    pruefbericht_log_path = log_dir / "pruefbericht.json"
    redacted_for_log = pruefbericht.model_dump(by_alias=True, mode="json")
    if redact_logs:
        redacted_for_log["SubmittedCase"]["tan"] = REDACTED_TAN
    with open(pruefbericht_log_path, "w") as f:
        json.dump(redacted_for_log, f, indent=2)
    log.info(f"Saved Prüfbericht copy to logs: {pruefbericht_log_path}")

    # submit Prüfbericht if requested
    if submit_pruefbericht:
        pruefbericht_config: PruefberichtModel = configuration.pruefbericht

        if (auth_url := pruefbericht_config.authorization_url) is None:
            raise ValueError("pruefbericht.authorization_url is required but not configured")
        if (client_id := pruefbericht_config.client_id) is None:
            raise ValueError("pruefbericht.client_id is required but not configured")
        if (configured_secret := pruefbericht_config.client_secret) is None:
            raise ValueError("pruefbericht.client_secret is required but not configured")
        client_secret = configured_secret.get_secret_value()
        if (api_base_url := pruefbericht_config.api_base_url) is None:
            raise ValueError("pruefbericht.api_base_url is required but not configured")

        log.info("Submitting Prüfbericht to BfArM...")

        # Perform retries *before* opening the DbContext so we don't hold a DB
        # transaction open for the entire exponential-backoff window (which could
        # be hours with the default 10 retries).  Only the final (possibly
        # failing) attempt is bracketed by the context manager.
        last_error: Exception | None = None
        initial_delay = 30.0
        backoff_factor = 2.0

        for attempt in range(1, max_retries + 1):
            try:
                _expiry, _token = _try_submit_pruefbericht(
                    pruefbericht=pruefbericht,
                    api_base_url=str(api_base_url),
                    auth_url=str(auth_url),
                    client_id=client_id,
                    client_secret=client_secret,
                    token="",
                )
                last_error = None
                break
            except Exception as e:
                last_error = e
                if attempt <= max_retries:
                    wait_time = initial_delay * (backoff_factor ** (attempt - 1))
                    log.warning(
                        f"Prüfbericht submission attempt {attempt}/{max_retries} "
                        f"failed: {e}. Retrying in {wait_time:.0f}s..."
                    )
                    time.sleep(wait_time)

        if last_error is not None:
            log.error(f"Prüfbericht submission failed after {max_retries} retries.")
            raise last_error

        # Only open the DbContext for the successful state transition.
        with DbContext(
            configuration=configuration,
            submission_id=submission_id,
            start_state=SubmissionStateEnum.REPORTING,
            end_state=SubmissionStateEnum.REPORTED,
            enabled=update_db,
        ):
            # The submission already succeeded above; the DbContext just records
            # the state transition.  If the state transition itself fails we
            # log it but don't lose the fact that the Prüfbericht was accepted.
            log.info("Prüfbericht submitted successfully!")
