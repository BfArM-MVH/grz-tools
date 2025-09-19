"""Command for generating various reports related to GRZ activities."""

import datetime
import logging
from pathlib import Path

import click
from grz_common.cli import config_file
from grz_db.models.submission import SubmissionStateEnum

from ..models.config import DbConfig
from .db.cli import get_submission_db_instance

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


@report.command()
@click.option(
    "--quarter",
    "quarter",
    type=click.IntRange(min=1, max=4),
    help="Quarter to generate report for.",
)
@click.option(
    "--year",
    "year",
    type=click.IntRange(min=2025),
    help="Year to generate report for.",
)
@click.option(
    "--outdir",
    "output_directory",
    type=click.Path(
        exists=True, file_okay=False, dir_okay=True, readable=True, writable=True, resolve_path=True, path_type=Path
    ),
    default=Path.cwd(),
    help="Directory to output TSV files. Defaults to current directory.",
)
@click.pass_context
def quarterly(ctx: click.Context, year: int | None, quarter: int | None, output_directory: Path):
    """
    Generate the tables for the quarterly report.
    """
    db = ctx.obj["db_url"]
    submission_db = get_submission_db_instance(db)

    if bool(year) != bool(quarter):
        raise click.UsageError("Both year and quarter must be provided or omitted.")
    elif (year and quarter) is None:
        today = datetime.date.today()
        quarter = ((today.month - 1) % 3) + 1
        # default to last quarter if ended less than 15 days ago otherwise current quarter
        if today <= datetime.date(year=today.year, month=1, day=15):
            year = today.year - 1
            quarter = 4
        else:
            year = today.year

        if (today.month in {4, 7, 10}) and (today.day <= 15):
            quarter -= 1

    log.info("Generating quarterly report for Q%d %d", quarter, year)
