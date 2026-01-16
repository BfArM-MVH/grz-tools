"""Command for dumping the configuration."""

import json
import logging
from pathlib import Path

import click
from grz_common.cli import config_file, config_files_from_ctx
from grz_common.utils.config import read_and_merge_config_files

log = logging.getLogger(__name__)


@click.command()
@config_file
@click.pass_context
def dump_config(ctx: click.Context, config_file: list[Path]):
    """
    Dump the merged configuration as read from config files.
    """
    config_files = config_files_from_ctx(ctx)
    log.info(f"Configuration files to load: {json.dumps([str(p.absolute()) for p in config_files], indent=2)}")

    config = read_and_merge_config_files(config_files)
    log.info(f"Merged configuration: {json.dumps(config, indent=2)}")

    log.info(
        "Note this only dumps the merged configuration as read from the files. "
        "It does not validate or process the configuration and ignores any environment variables."
    )
