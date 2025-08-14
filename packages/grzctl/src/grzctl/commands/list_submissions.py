"""Command for listing submissions."""

import json
import logging
import sys
from pathlib import Path

import click
import rich.console
import rich.table
import rich.text
from grz_common.cli import config_file, output_json
from grz_common.workers.download import InboxSubmissionState, InboxSubmissionSummary, query_submissions
from grz_db.models.submission import SubmissionDb, SubmissionStateEnum
from pydantic_core import to_jsonable_python

from ..models.config import ListConfig
from .db import get_submission_db_instance

log = logging.getLogger(__name__)


def _get_latest_state_txt(submission_db: SubmissionDb, submission_id: str) -> rich.text.Text:
    """
    Gets the latest database state of a submission ID as a Rich Text object,
    while also handling missing submissions and submissions without a state yet.
    """
    submission_from_db = submission_db.get_submission(submission_id)
    if submission_from_db:
        latest_state = submission_from_db.get_latest_state()
        if latest_state is None:
            latest_state_txt = rich.text.Text("none", style="italic yellow")
        else:
            latest_state_str = latest_state.state.value
            latest_state_txt = (
                rich.text.Text(latest_state_str, style="red")
                if latest_state_str == SubmissionStateEnum.ERROR
                else rich.text.Text(latest_state_str)
            )
    else:
        latest_state_txt = rich.text.Text("missing", style="italic yellow")

    return latest_state_txt


def _prepare_table(summaries: list[InboxSubmissionSummary], config: ListConfig) -> rich.table.Table:
    """
    Constructs a nice Rich Table to display inbox status and database state of submissions in the inbox.
    """
    table = rich.table.Table()
    table.add_column("ID", no_wrap=True)
    table.add_column("Inbox Status", no_wrap=True, justify="center")
    if config.db is not None:
        table.add_column("Database State", no_wrap=True, justify="center", style="green")
        submission_db = get_submission_db_instance(db_url=config.db.database_url)
    table.add_column("Oldest Upload", overflow="fold")
    table.add_column("Newest Upload", overflow="fold")
    for summary in summaries:
        match summary.state:
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
        row: list[rich.console.RenderableType] = [
            summary.submission_id,
            status_text,
            summary.oldest_upload.astimezone().strftime("%Y-%m-%d %H:%M:%S"),
            summary.newest_upload.astimezone().strftime("%Y-%m-%d %H:%M:%S"),
        ]
        if submission_db is not None:
            latest_state_txt = _get_latest_state_txt(submission_db, summary.submission_id)
            row.insert(2, latest_state_txt)
        table.add_row(*row)
    return table


def _validate_limit(ctx: click.Context, param: click.Parameter, value: int):
    if value < 0:
        raise click.BadParameter("limit must be a positive integer")

    return value


@click.command()
@config_file
@output_json
@click.option("--show-cleaned/--hide-cleaned", help="Show cleaned submissions.")
@click.option("--limit", type=int, default=10, callback=_validate_limit)
def list_submissions(config_file: Path, output_json: bool, show_cleaned: bool, limit: int):
    """
    List submissions within an inbox from oldest to newest, up to the requested limit.
    """
    config = ListConfig.from_path(config_file)
    submissions = query_submissions(config.s3, show_cleaned)

    if output_json:
        json.dump(to_jsonable_python(submissions), sys.stdout)
    else:
        console = rich.console.Console()
        table = _prepare_table(submissions[:limit], config)
        if len(submissions) > limit:
            console.print(f"[yellow]Limiting display to {limit} out of {len(submissions)} total submissions.[/yellow]")
        console.print(table)
