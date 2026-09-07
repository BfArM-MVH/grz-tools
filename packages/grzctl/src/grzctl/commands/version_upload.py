"""Command for publishing the version.json policy file to an inbox bucket."""

import logging
from typing import Any

import click
import grz_common.cli as grzcli
import requests
from grz_common.models.s3 import S3ConfigModel
from grz_common.models.version import VersionFile
from grz_common.transfer import init_s3_resource

log = logging.getLogger(__name__)

VERSION_FILE_URL = "https://raw.githubusercontent.com/BfArM-MVH/grz-tools/main/version_file/version.json"


@click.command()
@grzcli.configuration
def version_upload(configuration: dict[str, Any], **kwargs):
    """Download the canonical version.json policy from GitHub and upload it to an inbox bucket."""
    config = S3ConfigModel.model_validate(configuration)

    log.info(f"Fetching version policy from {VERSION_FILE_URL}")
    response = requests.get(VERSION_FILE_URL, timeout=30)
    response.raise_for_status()
    content = response.text

    # Validate before publishing to ensure that the version.json adheres to the schema.
    VersionFile.model_validate_json(content)

    s3_resource = init_s3_resource(config.s3)
    bucket = s3_resource.Bucket(config.s3.bucket)
    bucket.put_object(Key="version.json", Body=content.encode("utf-8"))

    message = f"Uploaded version.json to s3://{config.s3.bucket}/version.json"
    log.info(message)
    click.echo(message)
