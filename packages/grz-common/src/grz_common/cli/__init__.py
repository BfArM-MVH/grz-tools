"""
Common click options for the CLI commands.
"""

from os import sched_getaffinity
from pathlib import Path
from typing import Any

import click
import platformdirs

from ..utils.config import read_and_merge_config_files

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
    required=True,
    help="Path to the submission directory containing 'metadata/', 'files/', 'encrypted_files/' and 'logs/' directories",
)

config_file = click.option(
    "--config-file",
    metavar="STRING",
    type=FILE_R_E,
    required=False,
    multiple=True,
    default=[],
    help="Path to config file",
)


def get_default_config_path() -> Path:
    """
    Get the default config file path.

    :return: Path to the default config file
    """
    return Path(platformdirs.user_config_dir("grz-cli")) / "config.yaml"


def config_files_from_ctx(ctx: click.Context) -> list[Path]:
    """
    Helper to get config files from click context.

    :param ctx: click context
    :return: List of config file paths
    """
    retval: list[Path] = []
    if "config_file" in ctx.params:
        # add config files from current context if any
        retval = [
            *[Path(p) for p in ctx.params["config_file"]],
            *retval,
        ]

    while ctx.parent is not None:
        # traverse up the context tree and add config files from parent contexts
        ctx = ctx.parent
        if "config_file" in ctx.params:
            retval = [
                *[Path(p) for p in ctx.params["config_file"]],
                *retval,
            ]

    default_config_path = get_default_config_path()
    if default_config_path.is_file():
        # prepend default config path if it exists
        retval.insert(0, default_config_path)

    return retval


def read_config_from_ctx(ctx: click.Context) -> dict[str, Any]:
    """
    Helper to read and merge config files from click context.
    :param ctx: click context
    :return: Merged config dictionary
    """
    config_files = config_files_from_ctx(ctx)
    config = read_and_merge_config_files(config_files)
    return config


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
