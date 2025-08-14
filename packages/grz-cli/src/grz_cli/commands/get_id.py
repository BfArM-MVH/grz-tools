"""Command for submitting (validating, encrypting, and uploading) a submission."""

import io
import logging

import click
from grz_pydantic_models.submission.metadata import GrzSubmissionMetadata

log = logging.getLogger(__name__)


@click.command("get-id")
@click.argument("metadata", type=click.File("r"))
def get_id(metadata: io.TextIOWrapper):
    """
    Compute the submission ID for a given metadata JSON file.
    """
    click.echo(GrzSubmissionMetadata.model_validate_json(metadata.read()).submission_id)
