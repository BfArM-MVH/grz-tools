"""Command for decrypting a submission."""

import logging
import sys
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
@grzcli.encrypted_files_dir
@grzcli.logs_dir
@grzcli.force
@grzcli.update_db
def decrypt(  # noqa: PLR0913
    configuration: GrzctlConfig,
    submission_dir,
    metadata_dir,
    files_dir,
    encrypted_files_dir,
    logs_dir,
    force,
    update_db,
    **kwargs,
):
    """
    Decrypt a submission.

    Decrypting a submission requires the _private_ key of the original recipient.
    """
    bundled_mode = submission_dir is not None
    granular_mode = any(map(lambda v: v is not None, [metadata_dir, files_dir, encrypted_files_dir, logs_dir]))

    if bundled_mode and granular_mode:
        raise click.UsageError("'--submission-dir' is mutually exclusive with explicit path options.")

    if bundled_mode:
        base = Path(submission_dir)
        _metadata_dir = base / "metadata"
        _files_dir = base / "files"
        _encrypted_files_dir = base / "encrypted_files"
        _logs_dir = base / "logs"
    elif granular_mode:
        required = {
            "--metadata-dir": metadata_dir,
            "--files-dir": files_dir,
            "--encrypted-files-dir": encrypted_files_dir,
            "--logs-dir": logs_dir,
        }
        missing = [name for name, path in required.items() if path is None]
        if missing:
            raise click.UsageError(f"Granular mode requires: {', '.join(missing)}")
        _metadata_dir, _files_dir, _encrypted_files_dir, _logs_dir = (
            Path(metadata_dir),
            Path(files_dir),
            Path(encrypted_files_dir),
            Path(logs_dir),
        )
    else:
        raise click.UsageError("You must specify either '--submission-dir' or the required explicit path options.")

    keys = configuration.keys

    grz_privkey_path = keys.grz_private_key_path
    if not grz_privkey_path:
        log.error("GRZ private key path is required for decryption.")
        sys.exit(1)

    log.info("Starting decryption...")

    worker_inst = Worker(
        metadata_dir=_metadata_dir,
        files_dir=_files_dir,
        log_dir=_logs_dir,
        encrypted_files_dir=_encrypted_files_dir,
    )
    submission_id = worker_inst.parse_encrypted_submission().submission_id
    with DbContext(
        configuration=configuration,
        submission_id=submission_id,
        start_state=SubmissionStateEnum.DECRYPTING,
        end_state=SubmissionStateEnum.DECRYPTED,
        enabled=update_db,
    ):
        worker_inst.decrypt(grz_privkey_path, force=force)

    log.info("Decryption successful!")
