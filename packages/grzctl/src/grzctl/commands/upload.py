"""Command for uploading a submission."""

from pathlib import Path

import click
import grz_common.cli as grzcli
from grz_common.workers.worker import Worker
from grz_db.models.submission import SubmissionStateEnum

from ..commands import grzctl_configuration, inbox_option
from ..dbcontext import DbContext
from ..models.config import GrzctlConfig


@click.command()
@grzcli.submission_dir
@grzcli.threads
@grzctl_configuration
@inbox_option
@grzcli.update_db
def upload(
    configuration: GrzctlConfig,
    submission_dir,
    threads: int,
    inbox_name: str,
    update_db: bool,
    **kwargs,
):
    """Upload a submission to a GRZ/GDC (wrapper with DB updates)."""
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
    submitter_id = submission.metadata.content.submission.submitter_id

    s3_options = configuration.resolve_inbox(submitter_id=submitter_id, inbox_name=inbox_name).s3

    with DbContext(
        configuration=configuration,
        submission_id=submission_id,
        start_state=SubmissionStateEnum.UPLOADING,
        end_state=SubmissionStateEnum.UPLOADED,
        enabled=update_db,
    ):
        uploaded_id = worker_inst.upload(s3_options)
        click.echo(uploaded_id)
