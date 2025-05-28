"""Command for managing a submission database"""

import json
import sys

import click
import rich.console
import rich.table
from grz_db import (
    DatabaseConfigurationError,
    DuplicateSubmissionError,
    DuplicateTanGError,
    SubmissionDb,
    SubmissionNotFoundError,
    SubmissionStateEnum,
    SubmissionStateLog,
)

from .common import output_json

console = rich.console.Console()

DATABASE_URL = "sqlite:///test.sqlite"


def get_submission_db_instance(db_url: str | None) -> SubmissionDb:
    """Creates and returns an instance of SubmissionDb."""
    db_url = db_url or DATABASE_URL
    return SubmissionDb(db_url=db_url)


@click.group(help="Database operations")
def db():
    """Database operations"""
    pass


@db.group()
def submission():
    """Submission operations"""
    pass


db_url_option = click.option("--db", default=DATABASE_URL, envvar="GRZ_DATABASE")


@db.command()
@db_url_option
@click.option("--revision", default="head", help="Alembic revision to upgrade to (default: 'head').")
@click.option(
    "--alembic-ini",
    type=click.Path(exists=True, dir_okay=False, resolve_path=True),
    default="alembic.ini",
    help="Override path to alembic.ini file.",
)
def init(revision: str, alembic_ini: str | None, db: str):
    """Initializes or upgrades the database schema using Alembic."""
    submission_db = get_submission_db_instance(db)

    console.print(f"[cyan]Using alembic configuration: {alembic_ini}[/cyan]")

    try:
        console.print(f"[cyan]Attempting to upgrade database to revision: {revision}...[/cyan]")
        _ = submission_db.initialize_schema(
            alembic_ini_path=alembic_ini,
            revision=revision,
        )

    except (DatabaseConfigurationError, RuntimeError) as e:
        console.print(f"[red]Error during schema initialization: {e}[/red]")
        if isinstance(e, RuntimeError):
            console.print(
                "[yellow]Ensure your database is running and accessible, and alembic.ini is configured correctly.[/yellow]"
            )
            console.print(
                "[yellow]You might need to create an initial migration if this is the first time: 'alembic revision -m \"initial\" --autogenerate'[/yellow]"
            )
        raise click.ClickException(str(e)) from e
    except Exception as e:
        console.print(f"[red]An unexpected error occurred during 'db init': {type(e).__name__} - {e}[/red]")
        raise click.ClickException(str(e)) from e


@db.command("list")
@db_url_option
@output_json
def list_submissions(db: str, output_json: bool = False):
    """Lists all submissions in the database with their latest state."""
    db_service = get_submission_db_instance(db)
    submissions = db_service.list_submissions()

    if not submissions:
        console.print("[yellow]No submissions found in the database.[/yellow]")
        return

    table = rich.table.Table(title="All Submissions")
    table.add_column("ID", style="dim", width=12)
    table.add_column("tanG", style="cyan")
    table.add_column("Pseudonym", style="magenta")
    table.add_column("Latest State", style="green")
    table.add_column("Last State Timestamp (UTC)", style="yellow")

    submission_dicts = []

    for submission in submissions:
        latest_state_obj: SubmissionStateLog | None = None
        if submission.states:
            latest_state_obj = max(submission.states, key=lambda s: s.timestamp)

        if output_json:
            submission_dict = {
                "id": submission.id,
                "tan_g": submission.tan_g,
                "pseudonym": submission.pseudonym,
                "latest_state": None,
            }
            if latest_state_obj:
                submission_dict["latest_state"] = {
                    "state": latest_state_obj.state.value,
                    "timestamp": latest_state_obj.timestamp.isoformat(),
                    "data": latest_state_obj.data,
                }
            submission_dicts.append(submission_dict)
        else:
            latest_state_str = latest_state_obj.state.value if latest_state_obj else "N/A"
            latest_state_str = (
                f"[red]{latest_state_str}[/red]" if latest_state_str == SubmissionStateEnum.ERROR else latest_state_str
            )
            latest_timestamp_str = latest_state_obj.timestamp.isoformat() if latest_state_obj else "N/A"

            table.add_row(
                submission.id,
                submission.tan_g if submission.tan_g is not None else "N/A",
                submission.pseudonym if submission.pseudonym is not None else "N/A",
                latest_state_str,
                latest_timestamp_str,
            )

    if output_json:
        json.dump(submission_dicts, sys.stdout)
    else:
        console.print(table)


@submission.command()
@click.argument("submission_id", type=str)
@db_url_option
@click.option("--tan-g", "tan_g", type=str, default=None, help="The tanG for the submission.")
@click.option("--pseudonym", type=str, default=None, help="The pseudonym for the submission.")
def add(db: str, submission_id: str, tan_g: str | None, pseudonym: str | None):
    """
    Add a submission to the database.
    """
    db_service = get_submission_db_instance(db)
    try:
        db_submission = db_service.add_submission(submission_id, tan_g, pseudonym)
        console.print(f"[green]Submission '{db_submission.id}' added successfully.[/green]")
        console.print(f"  tanG: {db_submission.tan_g}, Pseudonym: {db_submission.pseudonym}")
    except (DuplicateSubmissionError, DuplicateTanGError) as e:
        console.print(f"[red]Error: {e}[/red]")
        raise click.Abort() from e
    except Exception as e:
        console.print(f"[red]An unexpected error occurred: {e}[/red]")
        raise click.ClickException(f"Failed to add submission: {e}") from e


@submission.command()
@click.argument("submission_id", type=str)
@click.argument(
    "state_str", metavar="STATE", type=click.Choice([s.value for s in SubmissionStateEnum], case_sensitive=False)
)
@db_url_option
@click.option("--data", "data_json", type=str, default=None, help='Additional JSON data (e.g., \'{"k":"v"}\').')
def update(submission_id: str, state_str: str, db: str, data_json: str | None):
    """Update a submission to the given state. Optionally accepts additional JSON data to associate with the log entry."""
    db_service = get_submission_db_instance(db)
    try:
        state_enum = SubmissionStateEnum(state_str)
    except ValueError as e:
        console.print(f"[red]Error: Invalid state value '{state_str}'.[/red]")
        raise click.Abort() from e

    parsed_data = None
    if data_json:
        try:
            parsed_data = json.loads(data_json)
        except json.JSONDecodeError as e:
            console.print(f"[red]Error: Invalid JSON string for --data: {data_json}[/red]")
            raise click.Abort() from e
    try:
        new_state_log = db_service.update_submission_state(submission_id, state_enum, parsed_data)
        console.print(
            f"[green]Submission '{submission_id}' updated to state '{new_state_log.state.value}'. Log ID: {new_state_log.id}[/green]"
        )
        if new_state_log.data:
            console.print(f"  Data: {new_state_log.data}")

        if state_enum == SubmissionStateEnum.REPORTED:
            updated_submission = db_service.get_submission(submission_id)
            if updated_submission:
                console.print(f"  Submission tanG is now: {updated_submission.tan_g}")

    except SubmissionNotFoundError as e:
        console.print(f"[red]Error: {e}[/red]")
        console.print(f"You might need to add it first: grz-cli db add-submission {submission_id}")
        raise click.Abort() from e
    except Exception as e:
        console.print(f"[red]An unexpected error occurred: {e}[/red]")
        raise click.ClickException(f"Failed to update submission state: {e}") from e


@submission.command("show")
@click.argument("submission_id", type=str)
@db_url_option
def show(submission_id: str, db: str):
    """
    Show details of a submission.
    """
    db_service = get_submission_db_instance(db)
    submission = db_service.get_submission(submission_id)
    if not submission:
        console.print(f"[red]Error: Submission with ID '{submission_id}' not found.[/red]")
        raise click.Abort()

    console.print(f"\n[bold blue]Submission Details for ID: {submission.id}[/bold blue]")
    console.print(f"  tanG: {submission.tan_g if submission.tan_g is not None else 'N/A'}")
    console.print(f"  Pseudonym: {submission.pseudonym if submission.pseudonym is not None else 'N/A'}")

    if submission.states:
        table = rich.table.Table(title=f"State History for Submission {submission.id}")
        table.add_column("Log ID", style="dim", width=12)
        table.add_column("Timestamp (UTC)", style="yellow")
        table.add_column("State", style="green")
        table.add_column("Data", style="cyan", overflow="ellipsis")

        sorted_states = sorted(submission.states, key=lambda s: s.timestamp)
        for state_log in sorted_states:
            data_str = json.dumps(state_log.data) if state_log.data else ""
            state = state_log.state.value
            state_str = f"[red]{state}[/red]" if state == SubmissionStateEnum.ERROR else state
            table.add_row(str(state_log.id), state_log.timestamp.isoformat(), state_str, data_str)
        console.print(table)
    else:
        console.print("[yellow]No state history found for this submission.[/yellow]")
