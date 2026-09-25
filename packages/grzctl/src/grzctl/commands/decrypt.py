"""Command for decrypting a submission."""

import logging
from pathlib import Path

import click
import grz_common.cli as grzcli
from grz_common.exceptions import ConfigurationError
from grz_common.workers.worker import Worker
from grz_db.models.submission import SubmissionStateEnum

from ..commands import grzctl_configuration, inbox_option
from ..dbcontext import DbContext
from ..models.config import GrzctlConfig
from .inbox_resolution import require_inbox

log = logging.getLogger(__name__)


@click.command()
@grzctl_configuration
@grzcli.submission_dir
@inbox_option
@grzcli.force
@grzcli.update_db
def decrypt(
    configuration: GrzctlConfig,
    submission_dir,
    inbox_name,
    force,
    update_db,
    **kwargs,
):
    """
    Decrypt a submission.

    Decrypting a submission requires the _private_ key of the original recipient.
    That is the private key of the inbox that the submission came from.
    grzctl looks up that inbox under the submitter named in the submission's metadata.
    The inbox is the one that --inbox names.
    Without --inbox, grzctl takes the inbox recorded in the database, else the submitter's only inbox.
    It reads the database only with --update-db.
    `grzctl download` records the inbox, and `grzctl db backfill` records it for older submissions.
    """
    log.info("Starting decryption...")

    submission_dir = Path(submission_dir)

    worker_inst = Worker(
        metadata_dir=submission_dir / "metadata",
        files_dir=submission_dir / "files",
        log_dir=submission_dir / "logs",
        encrypted_files_dir=submission_dir / "encrypted_files",
    )
    encrypted_submission = worker_inst.parse_encrypted_submission()
    submission_id = encrypted_submission.submission_id

    with DbContext(
        configuration=configuration,
        submission_id=submission_id,
        start_state=SubmissionStateEnum.DECRYPTING,
        end_state=SubmissionStateEnum.DECRYPTED,
        enabled=update_db,
    ) as db_context:
        submitter_id = encrypted_submission.metadata.content.submission.submitter_id
        # raise a ConfigurationError, so that the DbContext records the failure reason configuration_error
        inbox_name = require_inbox(
            configuration,
            submitter_id=submitter_id,
            submission_id=submission_id,
            inbox_name=inbox_name,
            db_service=db_context.db,
            hint="Pass --inbox, or record the inbox with 'grzctl db backfill'.",
            exc_type=ConfigurationError,
        )
        private_key = configuration.inbox_target(submitter_id, inbox_name).load_private_key()
        worker_inst.decrypt(recipient_private_key=private_key, force=force)

    log.info("Decryption successful!")
