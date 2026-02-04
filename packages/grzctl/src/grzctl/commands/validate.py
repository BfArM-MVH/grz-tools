"""Command for validating a submission."""

import logging
from pathlib import Path
from typing import Any

import click
import grz_cli.commands.validate as validate_module
import grz_common.cli as grzcli
from grz_common.workers.worker import Worker
from grz_db.models.submission import SubmissionStateEnum

from grzctl.dbcontext import DbContext

log = logging.getLogger(__name__)


@click.command()
@grzcli.configuration
@grzcli.submission_dir
@grzcli.force
@grzcli.threads
@click.option(
    "--with-grz-check/--no-grz-check",
    "with_grz_check",
    default=True,
    hidden=True,
    help="Whether to use grz-check to perform validation",
)
@grzcli.update_db
def validate(
    configuration: dict[str, Any],
    submission_dir,
    force,
    threads,
    with_grz_check,
    update_db,
    **kwargs,
):
    """
    Validate the submission (wrapper with DB updates).
    """
    submission_dir = Path(submission_dir)

    submission_id = ""
    if update_db:
        try:
            worker_inst = Worker(
                metadata_dir=submission_dir / "metadata",
                files_dir=submission_dir / "files",
                log_dir=submission_dir / "logs",
                encrypted_files_dir=submission_dir / "encrypted_files",
                threads=threads,
            )
            submission = worker_inst.parse_submission()
            submission_id = submission.metadata.content.submission_id
        except Exception:
            log.warning("Could not determine submission ID; database will not be updated.")
            update_db = False

    with DbContext(
        configuration=configuration,
        submission_id=submission_id,
        start_state=SubmissionStateEnum.VALIDATING,
        end_state=SubmissionStateEnum.VALIDATED,
        enabled=update_db,
    ):
        validate_module.validate.callback(  # type: ignore[misc]
            configuration=configuration,
            submission_dir=submission_dir,
            force=force,
            threads=threads,
            with_grz_check=with_grz_check,
            **kwargs,
        )
