"""Command for publishing the bundled version.json policy to an inbox bucket."""

import logging

import click
from grz_common.models.version import VersionFile
from grz_common.transfer import init_s3_resource

from ..commands import grzctl_configuration, inbox_option, submitter_id_option
from ..models.config import GrzctlConfig

log = logging.getLogger(__name__)


@click.command()
@grzctl_configuration
@submitter_id_option
@inbox_option
def version_push(configuration: GrzctlConfig, submitter_id: str, inbox_name: str, **kwargs):
    """Publish the bundled grz-cli version-compatibility policy to an inbox bucket."""
    s3_options = configuration.resolve_inbox(submitter_id=submitter_id, inbox_name=inbox_name).s3

    content = VersionFile.read_bundled_text()

    s3_resource = init_s3_resource(s3_options)
    bucket = s3_resource.Bucket(s3_options.bucket)
    bucket.put_object(Key="version.json", Body=content.encode("utf-8"))

    message = f"Published version.json to s3://{s3_options.bucket}/version.json"
    log.info(message)
    click.echo(message)
