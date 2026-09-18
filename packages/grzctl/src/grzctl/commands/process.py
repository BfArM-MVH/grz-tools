"""Command for processing a submission.

This module implements the ``grzctl process`` subcommand, which performs the
complete submission lifecycle in a single streaming pass:

    Download metadata → Decrypt → Validate → Re-encrypt → Archive → (optionally) Prüfbericht

It replaces the individual step-by-step subcommands (``download``, ``decrypt``,
``validate``, ``encrypt``, ``archive``) with a single streaming pipeline that
avoids materialising intermediate files on disk.

The pipeline stages are orchestrated by
:class:`grzctl.processor.SubmissionProcessor` and use an
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

The submission state transitions through ``PROCESSING → PROCESSED`` (or
``ERROR`` on failure).  The DB record is also *populated* with metadata so
that downstream Prüfbericht generation can read the required fields.  If
``--submit-pruefbericht`` is used, a separate ``REPORTING → REPORTED``
transition is recorded, or ``REPORTING → ERROR`` if BfArM never accepts it.

Recovery
--------
Each output of a file has its own progress log: ``progress_staging.cjson`` for
the validated copy staged in the interrogation bucket, ``progress_local.cjson``
for the decrypted copy on local storage.  Re-running the command after a
failure streams each file only into the outputs that are not recorded or whose
copy is gone, making the pipeline effectively idempotent.
"""

import json
import logging
import time
from pathlib import Path

import click
import grz_common.cli as grzcli
from grz_common.transfer import get_metadata_upload_timestamp, init_s3_client
from grz_common.workers.download import download_metadata_file
from grz_common.workers.submission import SubmissionMetadata
from grz_db.errors import DuplicateInitialSubmissionError
from grz_db.models.submission import CASE_LINK_KEY, SubmissionStateEnum
from grz_pydantic_models.pruefbericht.v0 import Pruefbericht
from grz_pydantic_models.submission.metadata import REDACTED_TAN

from ..commands import grzctl_configuration
from ..dbcontext import DbContext
from ..models.config import GrzctlConfig
from ..models.pruefbericht import PruefberichtModel
from ..processor import SubmissionProcessor
from .pruefbericht import _generate_pruefbericht_from_database, _get_submission_credentials
from .pruefbericht import _try_submit as _try_submit_pruefbericht
from .validate import _check_duplicate_initial

log = logging.getLogger(__name__)


@click.command()
@grzctl_configuration
@grzcli.submission_id
@grzcli.output_dir
@grzcli.threads
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
def process(  # noqa: PLR0913, PLR0917
    configuration: GrzctlConfig,
    submission_id: str,
    output_dir: str,
    threads: int,
    submit_pruefbericht: bool,
    save_pruefbericht: str | None,
    redact_pruefbericht: bool,
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

    The DB record is populated with the parsed metadata so that downstream
    Prüfbericht generation can read the required fields (submission date,
    donor info, etc.).

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
    s3_client = init_s3_client(inbox.s3)
    download_metadata_file(s3_client, inbox.s3.bucket, submission_id, metadata_dir)
    local_metadata_path = metadata_dir / "metadata.json"

    submission_metadata = SubmissionMetadata(local_metadata_path)

    processor = SubmissionProcessor(
        configuration=configuration,
        inbox=inbox,
        log_dir=log_dir,
        threads=threads,
        max_concurrent_uploads=concurrent_uploads,
        clean_inbox=clean_inbox,
    )

    with DbContext(
        configuration=configuration,
        submission_id=submission_id,
        start_state=SubmissionStateEnum.PROCESSING,
        end_state=SubmissionStateEnum.PROCESSED,
    ) as dbcontext_inst:
        # Populate the DB record with parsed metadata (donors, files, dates, etc.)
        # so that downstream Prüfbericht generation can read the required fields.
        # A rejected write, such as a duplicate tanG, then records the ERROR state.
        submission_date = get_metadata_upload_timestamp(s3_client, inbox.s3.bucket, submission_id).date()
        # A case link set by hand, with ``db case relink``, outlives a rerun. Resolving the case from
        # the metadata key again would undo that repair, and populate refuses to without --force.
        stored = dbcontext_inst.db.get_submission(submission_id)
        stored_case_id = stored.case_id if stored else None
        if stored_case_id is not None:
            log.warning(
                f"Submission '{submission_id}' keeps its case link (case {stored_case_id}); "
                "its metadata is not resolved to a case again."
            )
        dbcontext_inst.db.populate(
            submission_id,
            submission_metadata.content,
            submission_date,
            force=False,
            on_missing="create",
            ignore_fields={CASE_LINK_KEY} if stored_case_id is not None else None,
        )
        try:
            _check_duplicate_initial(dbcontext_inst.db, submission_metadata.content)
        except DuplicateInitialSubmissionError as e:
            # fail basic QC before any file is streamed, as ``grzctl validate`` does
            log.warning(f"{e} Failing basic QC for submission '{submission_id}' without processing.")
            dbcontext_inst.db.modify_submission(submission_id, "basic_qc_passed", False)
            raise
        processor.run(submission_metadata)

    _handle_pruefbericht(
        configuration=configuration,
        submission_id=submission_id,
        log_dir=log_dir,
        submit_pruefbericht=submit_pruefbericht,
        save_pruefbericht=save_pruefbericht,
        redact_pruefbericht=redact_pruefbericht,
    )


def _setup_directories(output_dir: str) -> tuple[Path, Path, Path]:
    """Create and return required directories."""
    base_dir = Path(output_dir)
    metadata_dir = base_dir / "metadata"
    log_dir = base_dir / "logs"

    for d in [base_dir, metadata_dir, log_dir]:
        d.mkdir(mode=0o770, parents=True, exist_ok=True)

    return base_dir, metadata_dir, log_dir


def _handle_pruefbericht(  # noqa: PLR0913, PLR0917
    configuration: GrzctlConfig,
    submission_id: str,
    log_dir: Path,
    submit_pruefbericht: bool,
    save_pruefbericht: str | None,
    redact_pruefbericht: bool,
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
        # the run archived the submission, so a Prüfbericht that cannot be generated is the one
        # thing left undone, and swallowing it would report the run as complete
        log.error(f"Failed to generate Prüfbericht: {e}")
        raise

    _save_pruefbericht(pruefbericht, log_dir, save_pruefbericht, redact_pruefbericht)

    if submit_pruefbericht:
        # Entering the context writes REPORTING and commits it, so the retries below hold no
        # transaction open. A Prüfbericht that never gets through is recorded as an error, and a
        # later ``grzctl pruefbericht submit`` records the reporting states again.
        with DbContext(
            configuration=configuration,
            submission_id=submission_id,
            start_state=SubmissionStateEnum.REPORTING,
            end_state=SubmissionStateEnum.REPORTED,
        ):
            _submit_pruefbericht_with_retries(pruefbericht, configuration.pruefbericht)
            log.info("Prüfbericht submitted successfully!")


def _save_pruefbericht(
    pruefbericht: Pruefbericht, log_dir: Path, save_pruefbericht: str | None, redact_pruefbericht: bool
) -> None:
    """Write the Prüfbericht to ``save_pruefbericht`` if given, and a copy with redacted TAN to ``log_dir``."""
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
    redacted_for_log["SubmittedCase"]["tan"] = REDACTED_TAN
    with open(pruefbericht_log_path, "w") as f:
        json.dump(redacted_for_log, f, indent=2)
    log.info(f"Saved Prüfbericht copy to logs: {pruefbericht_log_path}")


def _submit_pruefbericht_with_retries(pruefbericht: Pruefbericht, pruefbericht_config: PruefberichtModel) -> None:
    """Submit the Prüfbericht, retrying with exponential backoff; re-raise the last error."""
    auth_url, client_id, client_secret, api_base_url = _get_submission_credentials(pruefbericht_config)

    log.info("Submitting Prüfbericht to BfArM...")

    max_attempts = 10
    initial_delay = 30.0
    backoff_factor = 2.0

    for attempt in range(1, max_attempts + 1):
        try:
            _try_submit_pruefbericht(
                pruefbericht=pruefbericht,
                api_base_url=api_base_url,
                auth_url=auth_url,
                client_id=client_id,
                client_secret=client_secret,
                token="",
            )
            return
        except Exception as e:
            if attempt == max_attempts:
                log.error(f"Prüfbericht submission failed after {max_attempts} attempts.")
                raise
            wait_time = initial_delay * (backoff_factor ** (attempt - 1))
            log.warning(
                f"Prüfbericht submission attempt {attempt}/{max_attempts} failed: {e}. Retrying in {wait_time:.0f}s..."
            )
            time.sleep(wait_time)
