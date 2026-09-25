"""Command for encrypting a submission."""

import logging
from pathlib import Path

import click
import grz_common.cli as grzcli
from grz_common.workers.worker import Worker
from grz_db.models.submission import SubmissionStateEnum

from ..commands import grzctl_configuration
from ..dbcontext import DbContext
from ..models.config import GrzctlConfig
from .inbox_resolution import resolve_inbox

log = logging.getLogger(__name__)


@click.command()
@grzctl_configuration
@grzcli.submission_dir
@grzcli.force
@click.option(
    "--check-validation-logs/--no-check-validation-logs",
    "check_validation_logs",
    default=True,
    help="Check validation logs before encrypting.",
)
@grzcli.update_db
def encrypt(
    configuration: GrzctlConfig,
    submission_dir,
    force,
    check_validation_logs,
    update_db,
    **kwargs,
):
    """Encrypt a submission (standalone with DB updates).

    The files are encrypted for the archive that the submission's research consent selects.
    They are signed with the private key of the inbox that the submission came from.
    grzctl looks up that inbox under the submitter named in the submission's metadata.
    It takes the inbox recorded in the database, else the submitter's only inbox.
    It reads the database only with --update-db.
    `grzctl download` records the inbox, and `grzctl db backfill` records it for older submissions.
    If no inbox resolves, the files are signed with a random key.
    """
    submission_dir = Path(submission_dir)

    worker_inst = Worker(
        metadata_dir=submission_dir / "metadata",
        files_dir=submission_dir / "files",
        log_dir=submission_dir / "logs",
        encrypted_files_dir=submission_dir / "encrypted_files",
    )
    submission = worker_inst.parse_submission()
    submission_id = submission.metadata.content.submission_id

    submission_date = submission.metadata.content.submission.submission_date
    consented = submission.metadata.content.consents_to_research(submission_date)

    archive_target = configuration.archives.consented if consented else configuration.archives.non_consented

    with DbContext(
        configuration=configuration,
        submission_id=submission_id,
        start_state=SubmissionStateEnum.ENCRYPTING,
        end_state=SubmissionStateEnum.ENCRYPTED,
        enabled=update_db,
    ) as db_context:
        submitter_id = submission.metadata.content.submission.submitter_id
        inbox_name = resolve_inbox(
            configuration, submitter_id=submitter_id, submission_id=submission_id, db_service=db_context.db
        )
        if inbox_name is None:
            log.warning(
                f"No inbox resolves for submission {submission_id}, so its files are signed with a random key. "
                "To sign them with the private key of the submission's inbox, "
                "record the inbox with 'grzctl db backfill' and encrypt with --update-db."
            )
            signing_key = None
        else:
            # raise a ConfigurationError, so that the DbContext records the failure reason configuration_error
            signing_key = configuration.inbox_target(submitter_id, inbox_name).load_private_key()
        worker_inst.encrypt(
            recipient_public_key=archive_target.load_public_key(),
            submitter_private_key=signing_key,
            force=force,
            check_validation_logs=check_validation_logs,
        )
