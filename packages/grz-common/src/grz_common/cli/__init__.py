"""
Common click options for the CLI commands.
"""

from os import sched_getaffinity
from pathlib import Path

import click
import platformdirs

DEFAULT_CONFIG_PATH = Path(platformdirs.user_config_dir("grz-cli")) / "config.yaml"

# Aliases for path types for click options
# Naming convention: {DIR,FILE}_{Read,Write}_{Exists,Create}
DIR_R_E = click.Path(
    exists=True,
    file_okay=False,
    dir_okay=True,
    readable=True,
    writable=False,
    resolve_path=True,
)
DIR_RW_C = click.Path(
    exists=False,
    file_okay=False,
    dir_okay=True,
    readable=True,
    writable=True,
    resolve_path=True,
)
FILE_R_E = click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True)

submission_dir = click.option(
    "--submission-dir",
    metavar="PATH",
    type=DIR_R_E,
    required=False,
    help="Base directory for all submission components. Mutually exclusive with explicit path options.",
)

metadata_dir = click.option(
    "--metadata-dir",
    metavar="PATH",
    type=DIR_R_E,
    required=False,
    help="Path to the directory containing 'metadata.json'.",
)

files_dir = click.option(
    "--files-dir",
    metavar="PATH",
    type=DIR_R_E,
    required=False,
    help="Path to the directory containing the unencrypted data files.",
)

encrypted_files_dir = click.option(
    "--encrypted-files-dir",
    metavar="PATH",
    type=DIR_R_E,
    required=False,
    help="Path to the directory containing the 'encrypted_files/'.",
)

logs_dir = click.option(
    "--logs-dir",
    metavar="PATH",
    type=DIR_R_E,
    required=False,
    help="Path to the directory containing the log files.",
)

output_files_dir = click.option(
    "--output-files-dir",
    metavar="PATH",
    type=DIR_RW_C,
    required=False,
    help="Output directory where the 'files/' subdirectory will be created.",
)

output_encrypted_files_dir = click.option(
    "--output-encrypted-files-dir",
    metavar="PATH",
    type=DIR_RW_C,
    required=False,
    help="Output directory where the 'encrypted_files/' subdirectory will be created.",
)

output_logs_dir = click.option(
    "--output-logs-dir",
    metavar="PATH",
    type=DIR_RW_C,
    required=False,
    help="Output directory where the 'logs/' subdirectory will be created.",
)

config_file = click.option(
    "--config-file",
    metavar="STRING",
    type=FILE_R_E,
    required=False,
    default=DEFAULT_CONFIG_PATH,
    help="Path to config file",
)

threads = click.option(
    "--threads",
    default=min(len(sched_getaffinity(0)), 4),
    type=int,
    show_default=True,
    help="Number of threads to use for parallel operations",
)

submission_id = click.option(
    "--submission-id",
    required=True,
    type=str,
    metavar="STRING",
    help="S3 submission ID",
)

output_dir = click.option(
    "--output-dir",
    metavar="PATH",
    type=DIR_RW_C,
    required=True,
    default=None,
    help="Path to the target submission output directory",
)

output_json = click.option("--json", "output_json", is_flag=True, help="Output JSON for machine-readability.")

show_details = click.option("--details", "show_details", is_flag=True, help="Show more detailed output.")

force = click.option("--force/--no-force", help="Overwrite files and ignore cached results (dangerous!)")