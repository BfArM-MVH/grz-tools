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
@grzcli.submission_id
@grzcli.output_dir
@grzcli.threads
@grzcli.force
@grzcli.update_db
@click.option(
    '--populate/--no-populate', 
    default=True, 
    help="Update the submission metadata with information from metadata.json and S3. If combined with --force, will overwrite information in db without asking."
)
def download(  # noqa: PLR0913
    configuration: dict[str, Any],
    submission_id,
    output_dir,
    threads,
    force,
    update_db,
    populate,
    **kwargs,
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

    with DbContext(
        configuration=configuration,
        submission_id=submission_id,
        start_state=SubmissionStateEnum.DOWNLOADING,
        end_state=SubmissionStateEnum.DOWNLOADED,
        enabled=update_db,
    ) as db_context:
        worker_inst.download(config.s3, submission_id, force=force)
        if populate: worker_inst.populate(config.s3, db_context.db, submission_id, force_populate=force)

    log.info("Download finished!")
