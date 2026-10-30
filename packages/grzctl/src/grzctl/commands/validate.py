"""Command for validating a submission."""

import logging
from pathlib import Path

import click
import grz_common.cli as grzcli
from grz_common.workers.worker import Worker
from grz_db.errors import DuplicateInitialSubmissionError
from grz_db.models.submission import SubmissionStateEnum
from grz_pydantic_models.submission.metadata import GrzSubmissionMetadata

from ..commands import grzctl_configuration
from ..commands.db.cli import get_submission_db_instance
from ..dbcontext import DbContext
from ..models.config import GrzctlConfig
from .duplicate_checks import check_duplicate_initial

log = logging.getLogger(__name__)


def _warn_on_duplicate_initial(configuration: GrzctlConfig, metadata: GrzSubmissionMetadata) -> None:
    """Best-effort, read-only duplicate-initial check when DB updates are disabled.

    Warns instead of raising, since nothing is written to the database in this mode.

    A check that was attempted and could not be answered is worth saying out loud:
    silence here reads as "not a duplicate", and that is an answer nobody got.
    """
    try:
        db = get_submission_db_instance(configuration.db.database_url, author=None)
        check_duplicate_initial(db, metadata)
    except DuplicateInitialSubmissionError as e:
        log.warning(f"{e} Basic QC would be recorded as failed; continuing because DB updates are disabled.")
    except Exception as e:
        log.warning(f"Could not check whether this is a duplicate initial submission: {e}")


@click.command()
@grzctl_configuration
@grzcli.submission_dir
@grzcli.force
@grzcli.threads
@click.option(
    "--submitter-id",
    "submitter_id",
    required=True,
    type=str,
    metavar="STRING",
    help="Expected Leistungserbringer (LE) identifier for metadata validation.",
)
@click.option(
    "--mmap/--no-mmap",
    "mmap",
    default=False,
    hidden=True,
    help="Whether to use mmap.",
)
@grzcli.update_db
def validate(  # noqa: PLR0913, PLR0917
    configuration: GrzctlConfig,
    submission_dir,
    force,
    threads,
    submitter_id,
    mmap,
    update_db,
    **kwargs,
):
    """Validate the submission (standalone with DB updates)."""
    submission_dir = Path(submission_dir)

    worker_inst = Worker(
        metadata_dir=submission_dir / "metadata",
        files_dir=submission_dir / "files",
        log_dir=submission_dir / "logs",
        encrypted_files_dir=submission_dir / "encrypted_files",
        threads=threads,
    )
    submission = worker_inst.parse_submission()
    submission_id = submission.metadata.content.submission_id

    # Missing research consent is no validation error, but it decides the archive later on.
    submission_date = submission.metadata.content.submission.submission_date
    if missing_consent := submission.metadata.content.explain_no_research_consent(submission_date):
        log.info(
            f"Submission {submission_id} gives no research consent on its submission date ({submission_date}), "
            f"so it will go to the non-consented archive: {missing_consent}"
        )

    identifiers = configuration.identifiers.model_copy(update={"le": submitter_id})

    with DbContext(
        configuration=configuration,
        submission_id=submission_id,
        start_state=SubmissionStateEnum.VALIDATING,
        end_state=SubmissionStateEnum.VALIDATED,
        enabled=update_db,
    ) as dbcontext_inst:
        if update_db:
            try:
                check_duplicate_initial(dbcontext_inst.db, submission.metadata.content)
            except DuplicateInitialSubmissionError as e:
                # Another initial submission of this case already passed basic QC, so this
                # submission fails basic QC without spending any validation effort. Re-raising
                # lets DbContext record the ERROR state with the duplicate_initial reason.
                log.warning(f"{e} Failing basic QC for submission '{submission_id}' without validating.")
                _ = dbcontext_inst.db.modify_submission(submission_id, "basic_qc_passed", "false")
                raise
        else:
            _warn_on_duplicate_initial(configuration, submission.metadata.content)

        worker_inst.validate(identifiers=identifiers, force=force, no_mmap=not mmap)

        if update_db:
            try:
                _ = dbcontext_inst.db.modify_submission(submission_id, "basic_qc_passed", "true")
            except DuplicateInitialSubmissionError as e:
                # Catches a competing initial submission of this case that passed basic QC
                # while this submission was being validated; handled the same way as the
                # pre-check above.
                log.warning(
                    f"Submission '{submission_id}' data validated, but {e} Failing basic QC for this submission."
                )
                _ = dbcontext_inst.db.modify_submission(submission_id, "basic_qc_passed", "false")
                raise
