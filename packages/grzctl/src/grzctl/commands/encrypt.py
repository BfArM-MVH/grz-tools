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

    The files are encrypted for the archive that the submission's research consent selects,
    and signed with the key in archives.signing_key or archives.signing_key_path.
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
    ):
        worker_inst.encrypt(
            recipient_public_key=archive_target.load_public_key(),
            submitter_private_key=configuration.archives.load_signing_key(),
            force=force,
            check_validation_logs=check_validation_logs,
        )
