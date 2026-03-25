"""Command for downloading a submission."""

import logging
from pathlib import Path
from typing import Any

import click
import grz_common.cli as grzcli
from grz_common.workers.worker import Worker
from grz_db.models.submission import SubmissionStateEnum

from ..dbcontext import DbContext
from ..models.config import DownloadConfig

log = logging.getLogger(__name__)


@click.command()
@grzcli.configuration
@click.option("--submission-id", required=True, type=str, metavar="STRING", help="S3 submission ID")
@click.option(
    "--output-dir",
    "output_dir_base",
    metavar="PATH",
    type=grzcli.DIR_RW_C,
    required=False,
    help="Base output directory for all components.",
)
@click.option("--metadata-dir", type=grzcli.DIR_RW_C, required=False)
@click.option("--encrypted-files-dir", type=grzcli.DIR_RW_C, required=False)
@grzcli.logs_dir
@grzcli.threads
@grzcli.force
@grzcli.update_db
def download(  # noqa: PLR0913
    configuration: dict[str, Any],
    submission_id,
    output_dir_base,
    metadata_dir,
    encrypted_files_dir,
    logs_dir,
    threads,
    force,
    update_db,
    **kwargs,
):
    """
    Download a submission from a GRZ.
    """
    bundled_mode = output_dir_base is not None
    granular_mode = any([metadata_dir, encrypted_files_dir, logs_dir])

    if bundled_mode and granular_mode:
        raise click.UsageError("'--output-dir' is mutually exclusive with explicit path options.")

    if bundled_mode:
        base = Path(output_dir_base)
        if not base.is_dir():
            log.debug("Creating submission directory %s", base)
            base.mkdir(mode=0o770, parents=False, exist_ok=False)
        _metadata_dir = base / "metadata"
        _encrypted_files_dir = base / "encrypted_files"
        _logs_dir = base / "logs"
    elif granular_mode:
        required = {
            "--metadata-dir": metadata_dir,
            "--encrypted-files-dir": encrypted_files_dir,
            "--logs-dir": logs_dir,
        }
        missing = [name for name, path in required.items() if path is None]
        if missing:
            raise click.UsageError(f"Flexible mode requires: {', '.join(missing)}")
        _metadata_dir, _encrypted_files_dir, _logs_dir = Path(metadata_dir), Path(encrypted_files_dir), Path(logs_dir)
    else:
        raise click.UsageError("You must specify either '--output-dir' or all explicit path options.")

    config = DownloadConfig.model_validate(configuration)

    log.info("Starting download...")

    worker_inst = Worker(
        metadata_dir=_metadata_dir,
        files_dir=Path("/dev/null"),
        log_dir=_logs_dir,
        encrypted_files_dir=_encrypted_files_dir,
        threads=threads,
    )

    with DbContext(
        configuration=configuration,
        submission_id=submission_id,
        start_state=SubmissionStateEnum.DOWNLOADING,
        end_state=SubmissionStateEnum.DOWNLOADED,
        enabled=update_db,
    ):
        worker_inst.download(config.s3, submission_id, force=force)

    log.info("Download finished!")
