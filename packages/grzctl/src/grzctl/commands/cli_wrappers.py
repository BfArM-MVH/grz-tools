"""grzctl-specific wrappers for grz-cli commands that use the unified config."""

import logging
from pathlib import Path
from typing import Any

import click
import grz_cli.commands.encrypt as encrypt_module
import grz_cli.commands.validate as validate_module
import grz_common.cli as grzcli
from grz_common.models.s3 import S3Options
from grz_common.workers.worker import Worker
from grz_db.models.submission import SubmissionStateEnum

from ..dbcontext import DbContext
from ..models.config import GrzctlConfig

log = logging.getLogger(__name__)


@click.command()
@grzcli.configuration
@grzcli.submission_dir
@grzcli.force
@grzcli.threads
@click.option(
    "--mmap/--no-mmap",
    "mmap",
    default=False,
    hidden=True,
    help="Whether to use mmap.",
)
@grzcli.update_db
def validate(  # noqa: PLR0913
    configuration: dict[str, Any],
    submission_dir,
    force,
    threads,
    mmap,
    update_db,
    **kwargs,
):
    """Validate the submission (wrapper with DB updates)."""
    submission_dir = Path(submission_dir)

    submission_id = ""
    if update_db:
        worker_inst = Worker(
            metadata_dir=submission_dir / "metadata",
            files_dir=submission_dir / "files",
            log_dir=submission_dir / "logs",
            encrypted_files_dir=submission_dir / "encrypted_files",
            threads=threads,
        )
        submission = worker_inst.parse_submission()
        submission_id = submission.metadata.content.submission_id

    with DbContext(
        configuration=configuration,
        submission_id=submission_id,
        start_state=SubmissionStateEnum.VALIDATING,
        end_state=SubmissionStateEnum.VALIDATED,
        enabled=update_db,
    ) as dbcontext_inst:
        validate_module.validate.callback(  # type: ignore[misc]
            configuration=configuration,
            submission_dir=submission_dir,
            force=force,
            threads=threads,
            mmap=mmap,
            **kwargs,
        )

        if update_db:
            _ = dbcontext_inst.db.modify_submission(submission_id, "basic_qc_passed", "true")


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
    """Encrypt a submission (wrapper with DB updates)."""
    submission_dir = Path(submission_dir)

    submission_id = "unknown"
    if update_db:
        worker_inst = Worker(
            metadata_dir=submission_dir / "metadata",
            files_dir=submission_dir / "files",
            log_dir=submission_dir / "logs",
            encrypted_files_dir=submission_dir / "encrypted_files",
        )
        submission = worker_inst.parse_submission()
        submission_id = submission.metadata.content.submission_id

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


@click.command()
@grzcli.submission_dir
@grzcli.threads
@grzcli.configuration
@click.option(
    "--inbox-bucket",
    default=None,
    help="Inbox bucket name to upload to. Required when multiple inboxes are configured.",
)
@grzcli.update_db
def upload(
    configuration: dict[str, Any],
    submission_dir,
    threads,
    inbox_bucket,
    update_db,
    **kwargs,
):
    """Upload a submission to a GRZ/GDC (wrapper with DB updates)."""
    config = GrzctlConfig.model_validate(configuration)
    inbox_s3 = config.resolve_inbox_by_bucket(inbox_bucket)

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

    s3_options = S3Options(
        endpoint_url=inbox_s3.endpoint_url,
        bucket=inbox_s3.bucket,
        access_key=inbox_s3.access_key,
        secret=inbox_s3.secret,
    )

    with DbContext(
        configuration=configuration,
        submission_id=submission_id,
        start_state=SubmissionStateEnum.UPLOADING,
        end_state=SubmissionStateEnum.UPLOADED,
        enabled=update_db,
    ):
        uploaded_id = worker_inst.upload(s3_options)
        click.echo(uploaded_id)
