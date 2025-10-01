import sys
from importlib.metadata import version

import click
from grz_common.models.s3 import S3Options
from grz_common.transfer import get_version_info


def check_version_and_exit_if_needed(s3_options: S3Options, version_file_path: str = "version.json") -> None:
    """Check grz-cli version against requirements and exit if update is required"""
    # fetch the version information
    version_info = get_version_info(s3_options, version_file_path)
    current_version = version("grz-cli")

    # handle missing version info gracefully
    if version_info is None:
        click.echo("Warning: Could not retrieve version information from S3.")
        click.echo("Skipping version check.")
        return

    # normal version checks
    if version_info.minimal_version < current_version < version_info.latest_version:
        click.echo("Consider upgrading grz-cli to the latest version")
    elif current_version <= version_info.minimal_version:
        click.echo("Aborting! grz-cli version is too old and must be upgraded!")
        sys.exit(1)
    else:
        click.echo("grz-cli is up to date.")

