"""Command for listing submissions."""

import json
import logging
import sys

import click
import rich.console
import rich.table
import rich.text
from grz_common.cli import config_file, output_json
from grz_common.workers.download import InboxSubmissionState, query_submissions
from pydantic_core import to_jsonable_python

from ..models.config import ListConfig

log = logging.getLogger(__name__)


@click.command()
@config_file
@output_json
@click.option("--show-cleaned/--hide-cleaned", help="Show cleaned submissions.")
def list_submissions(config_file, output_json, show_cleaned):
    """
    List submissions within an inbox from oldest to newest.
    """
    config = ListConfig.from_path(config_file)
    submissions = query_submissions(config.s3, show_cleaned)

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
            match submission.state:
                case InboxSubmissionState.INCOMPLETE:
                    status_text = rich.text.Text("Incomplete", style="yellow")
                case InboxSubmissionState.COMPLETE:
                    status_text = rich.text.Text("Complete", style="green")
                case InboxSubmissionState.CLEANING:
                    status_text = rich.text.Text("Cleaning", style="yellow")
                case InboxSubmissionState.CLEANED:
                    status_text = rich.text.Text("Cleaned", style="sky_blue1")
                case InboxSubmissionState.ERROR:
                    status_text = rich.text.Text("Error", style="red")
                case _:
                    status_text = rich.text.Text("Unknown", style="red")
            table.add_row(
                submission.submission_id,
                status_text,
                submission.oldest_upload.astimezone().strftime("%Y-%m-%d %H:%M:%S"),
                submission.newest_upload.astimezone().strftime("%Y-%m-%d %H:%M:%S"),
            )
        console.print(table)
