"""Command for processing a submission."""

import json
import logging
import re
from pathlib import Path
from typing import Any

import click
import crypt4gh.keys
import grz_common.cli as grzcli
from grz_common.pipeline.processor import SubmissionProcessor
from grz_common.utils.crypt import Crypt4GH
from grz_common.workers.download import S3BotoDownloadWorker
from grz_common.workers.submission import EncryptedSubmission, SubmissionMetadata
from grz_db.errors import DuplicateSubmissionError, DuplicateTanGError
from grz_db.models.submission import SubmissionStateEnum

from ..dbcontext import DbContext
from ..models.config import ProcessConfig
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
    "--validate/--no-validate",
    "validate",
    default=True,
    help="Enable or disable validation of files during processing.",
)
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
def process(
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
    **kwargs,
):
    """
    Process a submission through the streaming pipeline.
    """
    config = ProcessConfig.model_validate(configuration)

    # select correct S3 inbox options
    s3_options = config.s3
    if config.s3.inboxes:
        le_id = submission_id.split("_")[0]
        if le_id in config.s3.inboxes:
            submitter_inboxes = config.s3.inboxes[le_id]
            if inbox_bucket:
                if inbox_bucket in submitter_inboxes:
                    s3_options = config.s3.model_copy(update={"bucket": inbox_bucket})
                else:
                    raise click.ClickException(f"Inbox bucket '{inbox_bucket}' not found for submitter '{le_id}'")
            elif len(submitter_inboxes) == 1:
                # automatically use the only available bucket
                bucket_name = next(iter(submitter_inboxes))
                s3_options = config.s3.model_copy(update={"bucket": bucket_name})
            else:
                inbox_names = ", ".join(submitter_inboxes.keys())
                raise click.ClickException(
                    f"Multiple inboxes found for submitter '{le_id}' ({inbox_names}), please specify --inbox-bucket"
                )
        elif not config.s3.bucket:
            raise click.ClickException(f"No inboxes found for submitter '{le_id}' and no default bucket configured.")

    if not s3_options.bucket:
        raise click.ClickException("S3 bucket is required. Specify it in config or via --inbox-bucket.")

    log.info(f"Starting streaming pipeline for submission: {submission_id}")

    # setup directories
    submission_dir_path = Path(output_dir)
    metadata_dir = submission_dir_path / "metadata"
    log_dir = submission_dir_path / "logs"
    encrypted_files_dir = submission_dir_path / "encrypted_files"

    for d in [submission_dir_path, metadata_dir, log_dir, encrypted_files_dir]:
        d.mkdir(mode=0o770, parents=True, exist_ok=True)

    # first, download metadata to understand the submission structure
    log.info("Downloading metadata...")
    download_worker = S3BotoDownloadWorker(
        s3_options, status_file_path=log_dir / "progress_download.cjson", threads=threads
    )
    download_worker.download_metadata(submission_id, metadata_dir, metadata_file_name="metadata.json")
    local_metadata_path = metadata_dir / "metadata.json"

    keys = {
        "private": Crypt4GH.retrieve_private_key(config.keys.grz_private_key_path),
        "consented_public": crypt4gh.keys.get_public_key(config.keys.consented_archive_public_key_path),
        "non_consented_public": crypt4gh.keys.get_public_key(config.keys.non_consented_archive_public_key_path),
    }

    submission_metadata = SubmissionMetadata(local_metadata_path)

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

    processor = SubmissionProcessor(
        configuration=config,
        source_s3_options=s3_options,
        keys=keys,
        redact_patterns=redact_patterns,
    )

    with DbContext(
        configuration=configuration,
        submission_id=submission_id,
        start_state=SubmissionStateEnum.PROCESSING,
        end_state=SubmissionStateEnum.PROCESSED,
        enabled=update_db,
    ):
        # Run the streaming pipeline
        processor.run(submission_metadata, threads=threads)

    # TODO
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


def _handle_pruefbericht(  # noqa: PLR0913
    config: ProcessConfig,
    configuration: dict[str, Any],
    encrypted_submission: EncryptedSubmission,
    submission_id: str,
    log_dir: Path,
    submit_pruefbericht: bool,
    save_pruefbericht: str | None,
    update_db: bool,
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
        save_path.write_text(json.dumps(redacted_data, indent=2))
        log.info(f"Saved Prüfbericht (with redacted TAN) to: {save_path}")

    # also save a copy to the logs directory (with redacted TAN)
    pruefbericht_log_path = log_dir / "pruefbericht.json"
    redacted_for_log = pruefbericht.model_dump(by_alias=True, mode="json")
    redacted_for_log["SubmittedCase"]["tan"] = "<REDACTED>"
    pruefbericht_log_path.write_text(json.dumps(redacted_for_log, indent=2))
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

        with DbContext(
            configuration=configuration,
            submission_id=submission_id,
            start_state=SubmissionStateEnum.REPORTING,
            end_state=SubmissionStateEnum.REPORTED,
            enabled=update_db,
        ):
            _expiry, _token = _try_submit_pruefbericht(
                pruefbericht=pruefbericht,
                api_base_url=str(api_base_url),
                auth_url=str(auth_url),
                client_id=client_id,
                client_secret=client_secret,
                token="",
            )

        log.info("Prüfbericht submitted successfully!")
