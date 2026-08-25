"""Command for validating a submission."""

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
def validate(  # noqa: PLR0913
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

    identifiers = configuration.identifiers.model_copy(update={"le": submitter_id})

    with DbContext(
        configuration=configuration,
        submission_id=submission_id,
        start_state=SubmissionStateEnum.VALIDATING,
        end_state=SubmissionStateEnum.VALIDATED,
        enabled=update_db,
    ) as dbcontext_inst:
        worker_inst.validate(identifiers=identifiers, force=force, no_mmap=not mmap)

        if update_db:
            _ = dbcontext_inst.db.modify_submission(submission_id, "basic_qc_passed", "true")
