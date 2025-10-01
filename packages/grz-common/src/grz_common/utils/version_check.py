import sys  
from importlib.metadata import version  
  
import click  
from grz_common.models.s3 import S3Options  
  

def check_version_and_exit_if_needed(s3_options: S3Options, version_file_path: str = "version.json") -> None:
    """Check grz-cli version against requirements and exit if update is required"""
    from grz_common.transfer import get_version_info

    # fetch the version information
    minimal_version, latest_version = get_version_info(s3_options, version_file_path)
    current_version = version("grz-cli")

    # handle missing version info gracefully
    if minimal_version is None or latest_version is None:
        click.echo("Warning: Could not retrieve version information from S3.")
        click.echo("Skipping version check.")
        return

    # normal version checks
    if minimal_version < current_version < latest_version:
        click.echo("Consider upgrading grz-cli to the latest version")
    elif current_version <= minimal_version:
        click.echo("Aborting! grz-cli version is too old and must be upgraded!")
        sys.exit(1)