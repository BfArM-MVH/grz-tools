"""Command for generating various reports related to GRZ activities."""

import datetime
import logging

import click
from grz_common.cli import config_file
from grz_db.models.submission import SubmissionStateEnum

from ..models.config import DbConfig
from .db import get_submission_db_instance

log = logging.getLogger(__name__)


@click.group()
@config_file
@click.pass_context
def report(ctx: click.Context, config_file: str):
    """
    Generate various reports related to GRZ activities.
    """
    config = DbConfig.from_path(config_file).db
    if not config:
        raise ValueError("DB config not found")

    ctx.obj = {"db_url": config.database_url}


@report.command()
@click.option(
    "--since",
    "since",
    type=datetime.date.fromisoformat,
    help="First date on which to include submissions (default: a week before 'until').",
)
@click.option(
    "--until",
    "until",
    type=datetime.date.fromisoformat,
    help="Last date on which to include submissions (default: today).",
)
@click.option("-s", "--separator", type=str, default="\t", help="Separator between columns (default: tab).")
@click.pass_context
def processed(ctx: click.Context, since: datetime.date | None, until: datetime.date | None, separator: str):
    """
    Generate a report of processed submissions.
    Generally, this is for regular reporting to LEs.
    """
    db = ctx.obj["db_url"]
    submission_db = get_submission_db_instance(db)

    if until is None:
        # default to today
        until = datetime.date.today()

    if since is None:
        # default to a week before 'until'
        since = until - datetime.timedelta(weeks=1)

    submissions = submission_db.list_processed_between(start=since, end=until)

    status_map: dict[bool | None, str] = {
        True: "yes",
        False: "no",
        None: "",
    }

    click.echo(f"# Submissions processed between {since} and {until}")
    click.echo(separator.join(["Submission ID", "Basic QC Passed", "Detailed QC Passed", "Prüfbericht Submitted"]))
    for submission in submissions:
        last_reported_state_change = submission.get_latest_state(filter_to_type=SubmissionStateEnum.REPORTED)
        if not last_reported_state_change:
            continue
        click.echo(
            separator.join(
                [
                    submission.id,
                    status_map[submission.basic_qc_passed],
                    status_map[submission.detailed_qc_passed],
                    last_reported_state_change.timestamp.date().isoformat(),
                ]
            )
        )
