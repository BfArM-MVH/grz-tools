"""Command for listing submissions."""

import datetime
import json
import logging
import sys
from typing import Any

import click
import grz_common.cli as grzcli
import rich.console
import rich.table
import rich.text
from grz_common.workers.download import InboxSubmissionState, InboxSubmissionSummary, query_submissions
from grz_db.models.submission import SubmissionDb
from pydantic_core import to_jsonable_python

from ..models.config import GrzctlConfig
from . import limit
from .db.cli import get_submission_db_instance

log = logging.getLogger(__name__)


BYTES_PER_GIGABYTE = 1_000_000_000


def _get_latest_state_str(submission_db: SubmissionDb, submission_id: str) -> str | None:
    submission_from_db = submission_db.get_submission(submission_id)
    if submission_from_db:
        latest_state_log = submission_from_db.get_latest_state()
        latest_state = None if latest_state_log is None else latest_state_log.state.value
    else:
        latest_state = "missing"

    return latest_state


def _get_latest_state_txt(latest_state: str | None) -> rich.text.Text:
    """
    Gets the latest database state of a submission ID as a Rich Text object,
    while also handling missing submissions and submissions without a state yet.
    """
    if latest_state is None:
        latest_state_txt = rich.text.Text("none", style="italic yellow")
    elif latest_state == "missing":
        latest_state_txt = rich.text.Text("missing", style="italic yellow")
    else:
        latest_state_txt = (
            rich.text.Text(latest_state, style="red") if latest_state == "Error" else rich.text.Text(latest_state)
        )

    return latest_state_txt


def _format_upload_duration(duration: datetime.timedelta) -> rich.text.Text:
    upload_hours_remainder = duration.seconds // 3600
    upload_minutes_remainder, upload_seconds_remainder = divmod(duration.seconds - (upload_hours_remainder * 3600), 60)

    upload_duration_str = f"{duration.days}D" if duration.days else ""
    if upload_hours_remainder or upload_duration_str:
        upload_duration_str += f"{upload_hours_remainder}H"
    if upload_minutes_remainder or upload_duration_str:
        upload_duration_str += f"{upload_minutes_remainder}M"
    if upload_seconds_remainder or upload_duration_str:
        upload_duration_str += f"{upload_seconds_remainder}S"

    return (
        rich.text.Text(upload_duration_str)
        if upload_duration_str
        else rich.text.Text("Instantaneous 🚀", style="sky_blue1")
    )


def _prepare_table(
    summaries: list[InboxSubmissionSummary], database_states: dict[str, str | None] | None
) -> rich.table.Table:
    """
    Constructs a nice Rich Table to display inbox status and database state of submissions in the inbox.
    """
    table = rich.table.Table()
    table.add_column("ID", no_wrap=True)
    table.add_column("Inbox Status", no_wrap=True, justify="center")
    if database_states is not None:
        table.add_column("Database State", no_wrap=True, justify="center", style="green")
    table.add_column("Upload Duration", overflow="fold", justify="center")
    table.add_column("Newest Upload", overflow="fold")
    table.add_column("Size (GB)", justify="center")
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
            _format_upload_duration(summary.newest_upload - summary.oldest_upload),
            summary.newest_upload.astimezone().strftime("%Y-%m-%d %H:%M:%S"),
            f"{summary.total_size_bytes / BYTES_PER_GIGABYTE:.1f}",
        ]
        if database_states:
            latest_state_txt = _get_latest_state_txt(database_states[summary.submission_id])
            row.insert(2, latest_state_txt)
        table.add_row(*row)
    return table


@click.command()
@grzcli.configuration
@grzcli.output_json
@click.option("--show-cleaned/--hide-cleaned", help="Show cleaned submissions.")
@click.option(
    "--inbox-bucket",
    default=None,
    help="Inbox bucket name to list. Required when multiple inboxes are configured.",
)
@limit
def list_submissions(
    configuration: dict[str, Any],
    output_json: bool,
    show_cleaned: bool,
    inbox_bucket,
    limit: int,
    **kwargs,
):
    """
    List submissions within an inbox from oldest to newest, up to the requested limit.
    """
    config = GrzctlConfig.from_configuration(configuration)
    s3_options = config.resolve_inbox_by_bucket(inbox_bucket)

    submissions = query_submissions(s3_options, show_cleaned)

    database_states: dict[str, str | None] | None = None
    try:
        submission_db = get_submission_db_instance(db_url=config.db.database_url)
        database_states = {}
        for submission in submissions:
            database_states[submission.submission_id] = _get_latest_state_str(submission_db, submission.submission_id)
    except Exception:
        database_states = None
        log.warning("Could not query database for submission states. Run 'grzctl db upgrade' to initialize.")

    if output_json:
        submissions_jsonable = to_jsonable_python(submissions[:limit])
        if database_states is not None:
            for submission in submissions_jsonable:
                submission["database_state"] = database_states[submission["submission_id"]]
        json.dump(submissions_jsonable, sys.stdout)
    else:
        console = rich.console.Console()
        table = _prepare_table(submissions[:limit], database_states)
        if len(submissions) > limit:
            console.print(f"[yellow]Limiting display to {limit} out of {len(submissions)} total submissions.[/yellow]")
        console.print(table)
