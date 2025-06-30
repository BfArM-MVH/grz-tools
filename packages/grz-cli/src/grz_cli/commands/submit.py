"""Command for submitting (validating, encrypting, and uploading) a submission."""

import logging

import click

from .encrypt import encrypt
from .upload import upload
from .validate import validate

log = logging.getLogger(__name__)

from grz_common.cli import config_file, force, submission_dir, threads


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
    click.echo("Starting submission process...")
    ctx.invoke(validate, submission_dir=submission_dir, config_file=config_file, force=force)
    ctx.invoke(encrypt, submission_dir=submission_dir, config_file=config_file, force=force)
    ctx.invoke(upload, submission_dir=submission_dir, config_file=config_file, threads=threads)
    click.echo("Submission finished!")
