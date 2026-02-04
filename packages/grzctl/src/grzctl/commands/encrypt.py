"""Command for encrypting a submission."""

import logging
from pathlib import Path
from typing import Any

import click
import grz_cli.commands.encrypt as encrypt_module
import grz_common.cli as grzcli
from grz_common.workers.worker import Worker
from grz_db.models.submission import SubmissionStateEnum

from grzctl.dbcontext import DbContext

log = logging.getLogger(__name__)


@click.command()
@grzcli.configuration
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
    configuration: dict[str, Any],
    submission_dir,
    force,
    check_validation_logs,
    update_db,
    **kwargs,
):
    """
    Encrypt a submission (wrapper with DB updates).
    """
    submission_dir = Path(submission_dir)

    submission_id = "unknown"
    if update_db:
        try:
            worker_inst = Worker(
                metadata_dir=submission_dir / "metadata",
                files_dir=submission_dir / "files",
                log_dir=submission_dir / "logs",
                encrypted_files_dir=submission_dir / "encrypted_files",
            )
            submission = worker_inst.parse_submission()
            submission_id = submission.metadata.content.submission_id
        except Exception:
            log.warning("Could not determine submission ID; database will not be updated.")
            update_db = False

    with DbContext(
        configuration=configuration,
        submission_id=submission_id,
        start_state=SubmissionStateEnum.ENCRYPTING,
        end_state=SubmissionStateEnum.ENCRYPTED,
        enabled=update_db,
    ):
        encrypt_module.encrypt.callback(  # type: ignore[misc]
            configuration=configuration,
            submission_dir=submission_dir,
            force=force,
            check_validation_logs=check_validation_logs,
            **kwargs,
        )
