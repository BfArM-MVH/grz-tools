"""Command for submitting (validating, encrypting, and uploading) a submission."""

import logging

import click

from .encrypt import encrypt
from .upload import upload
from .validate import validate

log = logging.getLogger(__name__)

from .common import config_file, submission_dir, threads


@click.command("submit")
@submission_dir
@config_file
@threads
@click.pass_context
def submit(ctx, submission_dir, config_file, threads):
    """
    Validate, encrypt, and then upload.

    This is a convenience command that performs the following steps in order:
    1. Validate the submission
    2. Encrypt the submission
    3. Upload the encrypted submission
    """
    click.echo("Starting submission process...")
    ctx.invoke(validate, submission_dir=submission_dir)
    ctx.invoke(encrypt, submission_dir=submission_dir, config_file=config_file)
    ctx.invoke(upload, submission_dir=submission_dir, config_file=config_file, threads=threads)
    click.echo("Submission finished!")
