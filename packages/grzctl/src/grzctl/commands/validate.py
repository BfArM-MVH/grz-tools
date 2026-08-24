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
@grzcli.metadata_dir
@grzcli.files_dir
@grzcli.logs_dir
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
    metadata_dir,
    files_dir,
    logs_dir,
    force,
    threads,
    submitter_id,
    mmap,
    update_db,
    **kwargs,
):
    """Validate the submission (standalone with DB updates)."""
    bundled_mode = submission_dir is not None
    granular_mode = any(map(lambda v: v is not None, [metadata_dir, files_dir, logs_dir]))

    if bundled_mode and granular_mode:
        raise click.UsageError("'--submission-dir' is mutually exclusive with explicit path options.")

    if not bundled_mode and not granular_mode:
        raise click.UsageError("You must specify either '--submission-dir' or the required explicit path options.")

    _metadata_dir = Path(submission_dir) / "metadata" if bundled_mode else Path(metadata_dir)
    _files_dir = Path(submission_dir) / "files" if bundled_mode else Path(files_dir)
    _logs_dir = Path(submission_dir) / "logs" if bundled_mode else Path(logs_dir)
    _encrypted_files_dir = (
        Path(submission_dir) / "encrypted_files" if bundled_mode else _metadata_dir.parent / "encrypted_files"
    )

    worker_inst = Worker(
        metadata_dir=_metadata_dir,
        files_dir=_files_dir,
        log_dir=_logs_dir,
        encrypted_files_dir=_encrypted_files_dir,
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
