"""Command for submitting (validating, encrypting, and uploading) a submission."""

import logging

import click

from .encrypt import encrypt
from .upload import upload
from .validate import validate
from grz_common.utils.version_check import check_version_and_exit_if_needed

log = logging.getLogger(__name__)

from grz_common.cli import config_file, force, submission_dir, threads
from ..models.config import UploadConfig


@click.command("submit")
@submission_dir
@config_file
@threads
@force
@click.pass_context
def submit(ctx, submission_dir, config_file, threads, force):
    """
    Validate, encrypt, and then upload.

    This is a convenience command that performs the following steps in order:
    1. Validate the submission
    2. Encrypt the submission
    3. Upload the encrypted submission
    """

    config = UploadConfig.from_path(config_file)
    check_version_and_exit_if_needed(config.s3)

    click.echo("Starting submission process...")
    ctx.invoke(validate, submission_dir=submission_dir, config_file=config_file, force=force)
    ctx.invoke(encrypt, submission_dir=submission_dir, config_file=config_file, force=force)
    ctx.invoke(upload, submission_dir=submission_dir, config_file=config_file, threads=threads)
    click.echo("Submission finished!")
