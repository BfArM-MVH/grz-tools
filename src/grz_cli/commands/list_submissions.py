"""Command for listing submissions."""

import json
import logging
import sys

import click
import rich.console
import rich.table
import rich.text
from pydantic_core import to_jsonable_python

from ..utils.config import read_config
from ..workers.download import query_submissions
from .common import config_file, output_json

log = logging.getLogger(__name__)


@click.command()
@config_file
@output_json
def list_submissions(config_file, output_json):
    """
    List submissions within an inbox from oldest to newest.
    """
    config = read_config(config_file)
    submissions = query_submissions(config)

    if output_json:
        json.dump(to_jsonable_python(submissions), sys.stdout)
    else:
        console = rich.console.Console()
        table = rich.table.Table()
        table.add_column("ID", no_wrap=True)
        table.add_column("Status", no_wrap=True)
        table.add_column("Oldest Upload", overflow="fold")
        table.add_column("Newest Upload", overflow="fold")
        for submission in submissions:
            status_text = rich.text.Text(
                "Complete" if submission.complete else "Incomplete",
                style="green" if submission.complete else "yellow",
            )
            table.add_row(
                submission.submission_id,
                status_text,
                submission.oldest_upload.astimezone().strftime("%Y-%m-%d %H:%M:%S"),
                submission.newest_upload.astimezone().strftime("%Y-%m-%d %H:%M:%S"),
            )
        console.print(table)
