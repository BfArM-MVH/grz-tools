"""Command for publishing the bundled version.json policy to every configured inbox bucket."""

import logging
import sys

import botocore
import click
from grz_common.models.version import VERSION_FILE_KEY, VersionFile
from grz_common.transfer import init_s3_resource

from ..commands import grzctl_configuration
from ..models.config import GrzctlConfig

log = logging.getLogger(__name__)


@click.command()
@grzctl_configuration
def version_push(configuration: GrzctlConfig, **kwargs):
    """Publish the bundled grz-cli version-compatibility policy to every configured inbox bucket."""
    content = VersionFile.read_bundled_text()

    failures: list[tuple[str, str]] = []
    for le_id, entry in configuration.leistungserbringer.items():
        for inbox_name in entry.inbox_buckets:
            s3_options = configuration.resolve_inbox(submitter_id=le_id, inbox_name=inbox_name).s3
            target = f"s3://{s3_options.bucket}/{VERSION_FILE_KEY}"
            try:
                s3_resource = init_s3_resource(s3_options)
                s3_resource.Bucket(s3_options.bucket).put_object(Key=VERSION_FILE_KEY, Body=content.encode("utf-8"))
            except botocore.exceptions.ClientError as e:
                log.error(f"Failed to publish to {target} (LE {le_id}, inbox {inbox_name}): {e}")
                failures.append((le_id, inbox_name))
                continue

            click.echo(f"Published version.json to {target} (LE {le_id}, inbox {inbox_name})")

    if failures:
        failed_str = ", ".join(f"{le_id}/{inbox_name}" for le_id, inbox_name in failures)
        log.error(f"Failed to publish to {len(failures)} inbox(es): {failed_str}")
        sys.exit(1)
