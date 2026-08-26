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
@grzcli.metadata_dir
@grzcli.files_dir
@grzcli.output_encrypted_files_dir
@grzcli.logs_dir
@grzcli.force
@click.option(
    "--check-validation-logs/--no-check-validation-logs",
    "check_validation_logs",
    default=True,
    help="Check validation logs before encrypting.",
)
@grzcli.update_db
def encrypt(  # noqa: PLR0913
    configuration: GrzctlConfig,
    submission_dir,
    metadata_dir,
    files_dir,
    output_encrypted_files_dir,
    logs_dir,
    force,
    check_validation_logs,
    update_db,
    **kwargs,
):
    """Encrypt a submission (standalone with DB updates)."""
    bundled_mode = submission_dir is not None
    granular_mode = any(v is not None for v in [metadata_dir, files_dir, output_encrypted_files_dir, logs_dir])

    if bundled_mode and granular_mode:
        raise click.UsageError("'--submission-dir' is mutually exclusive with explicit path options.")

    if not bundled_mode and not granular_mode:
        raise click.UsageError("You must specify either '--submission-dir' or the required explicit path options.")

    if bundled_mode:
        base = Path(submission_dir)
        _metadata_dir = base / "metadata"
        _files_dir = base / "files"
        _logs_dir = base / "logs"
        _encrypted_files_dir = base / "encrypted_files"
    else:
        required = {
            "--metadata-dir": metadata_dir,
            "--files-dir": files_dir,
            "--output-encrypted-files-dir": output_encrypted_files_dir,
            "--logs-dir": logs_dir,
        }
        missing = [name for name, path in required.items() if path is None]
        if missing:
            raise click.UsageError(f"Granular mode requires: {', '.join(missing)}")
        _metadata_dir = Path(metadata_dir)
        _files_dir = Path(files_dir)
        _logs_dir = Path(logs_dir)
        _encrypted_files_dir = Path(output_encrypted_files_dir)

    worker_inst = Worker(
        metadata_dir=_metadata_dir,
        files_dir=_files_dir,
        log_dir=_logs_dir,
        encrypted_files_dir=_encrypted_files_dir,
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
            recipient_public_key_path=archive_target.public_key_path,
            submitter_private_key_path=configuration.keys.grz_private_key_path,
            force=force,
            check_validation_logs=check_validation_logs,
        )
