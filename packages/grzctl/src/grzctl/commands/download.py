"""Command for downloading a submission."""

import logging
from pathlib import Path
from typing import Any

import click
import grz_common.cli as grzcli
from grz_common.workers.worker import Worker

from ..models.config import DownloadConfig

log = logging.getLogger(__name__)


@click.command()
@grzcli.configuration
@grzcli.submission_id
@grzcli.output_dir
@grzcli.threads
@grzcli.force
def download(  # noqa: PLR0913
    configuration: dict[str, Any],
    config_file: tuple[Path],
    submission_id,
    output_dir,
    threads,
    force,
):
    """
    Download a submission from a GRZ.

    Downloaded metadata is stored within the `metadata` sub-folder of the submission output directory.
    Downloaded files are stored within the `encrypted_files` sub-folder of the submission output directory.
    """
    config = DownloadConfig.model_validate(configuration)

    log.info("Starting download...")

    submission_dir_path = Path(output_dir)
    if not submission_dir_path.is_dir():
        log.debug("Creating submission directory %s", submission_dir_path)
        submission_dir_path.mkdir(mode=0o770, parents=False, exist_ok=False)

    worker_inst = Worker(
        metadata_dir=submission_dir_path / "metadata",
        files_dir=submission_dir_path / "files",
        log_dir=submission_dir_path / "logs",
        encrypted_files_dir=submission_dir_path / "encrypted_files",
        threads=threads,
    )
    worker_inst.download(config.s3, submission_id, force=force)

    log.info("Download finished!")
