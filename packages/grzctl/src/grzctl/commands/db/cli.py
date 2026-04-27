"""Command for managing a submission database"""

import csv
import json
import logging
import sys
import traceback
from collections import namedtuple
from datetime import UTC, date, datetime, timedelta
from enum import StrEnum
from pathlib import Path
from typing import Any

import click
import grz_common.cli as grzcli
import rich.console
import rich.padding
import rich.panel
import rich.table
import rich.text
import textual.logging
from cryptography.hazmat.primitives.serialization import load_ssh_public_key
from grz_common.cli import output_json
from grz_common.logging import LOGGING_DATEFMT, LOGGING_FORMAT
from grz_common.workers.download import query_submissions
from grz_db.errors import (
    DatabaseConfigurationError,
    SubmissionError,
    SubmissionNotFoundError,
)
from grz_db.models.author import Author
from grz_db.models.submission import (
    ChangeRequestEnum,
    ChangeRequestLog,
    DetailedQCResult,
    FieldDiff,
    Submission,
    SubmissionDb,
    SubmissionDiffCollection,
    SubmissionStateEnum,
    SubmissionStateFilterModeEnum,
    SubmissionStateLog,
)
from grz_pydantic_models.common import StrictBaseModel
from grz_pydantic_models.submission.metadata import (
    REDACTED_LOCAL_CASE_ID,
    REDACTED_TAN,
    GenomicStudySubtype,
    GrzSubmissionMetadata,
    LibraryType,
    SequenceSubtype,
    SequenceType,
)
from pydantic import Field

from ...models.config import DbConfig, ListConfig
from .. import limit
from . import SignatureStatus, _verify_signature
from .sync import sync_submissions
from .tui import DatabaseBrowser

console = rich.console.Console()
console_err = rich.console.Console(stderr=True)
log = logging.getLogger(__name__)
_TEXT_MISSING = rich.text.Text("missing", style="italic yellow")


def get_submission_db_instance(db_url: str, author: Author | None = None) -> SubmissionDb:
    """Creates and returns an instance of SubmissionDb."""
    return SubmissionDb(db_url=db_url, author=author)


@click.group(help="Database operations")
@grzcli.configuration
@click.pass_context
def db(
    ctx: click.Context,
    configuration: dict[str, Any],
    **kwargs,
):
    """Database operations"""
    # set up context object
    ctx.ensure_object(dict)

    config = DbConfig.model_validate(configuration)
    db_config = config.db
    if not db_config:
        raise DatabaseConfigurationError("DB config not found")
    author_name = db_config.author.name

    if path := db_config.author.private_key_path:
        with open(path, "rb") as f:
            private_key_bytes = f.read()
    elif key := db_config.author.private_key:
        private_key_bytes = key.encode("utf-8")
    else:
        raise DatabaseConfigurationError("Either private_key or private_key_path must be provided.")

    log.debug("Reading known public keys...")
    KnownKeyEntry = namedtuple("KnownKeyEntry", ["key_format", "public_key_base64", "comment"])
    with open(db_config.known_public_keys) as f:
        public_key_list = list(map(lambda v: KnownKeyEntry(*v), map(lambda s: s.strip().split(), f.readlines())))
        public_keys = {
            comment: load_ssh_public_key(f"{fmt}\t{key}\t{comment}".encode()) for fmt, key, comment in public_key_list
        }
        for comment in public_keys:
            log.debug(f"Found public key labeled '{comment}'")

    author = Author(
        name=author_name,
        private_key_bytes=private_key_bytes,
        private_key_passphrase=db_config.author.private_key_passphrase,
    )
    ctx.obj.update(
        {
            "author": author,
            "public_keys": public_keys,
            "db_url": db_config.database_url,
        }
    )


@db.group()
@click.pass_context
def submission(ctx: click.Context):
    """Submission operations"""
    pass


@db.command()
@click.pass_context
def init(ctx: click.Context):
    """Initializes the database schema using Alembic."""
    db = ctx.obj["db_url"]
    submission_db = get_submission_db_instance(db, author=ctx.obj["author"])
    console_err.print(f"[cyan]Initializing database {db}[/cyan]")
    submission_db.initialize_schema()


@db.command()
@click.option("--revision", default="head", help="Alembic revision to upgrade to (default: 'head').")
@click.pass_context
def upgrade(
    ctx: click.Context,
    revision: str,
):
    """
    Upgrades the database schema using Alembic.
    """
    db = ctx.obj["db_url"]
    submission_db = get_submission_db_instance(db, author=ctx.obj["author"])

    try:
        revision_desc = "latest revision" if revision == "head" else f"revision '{revision}'"
        console_err.print(f"[cyan]Attempting to upgrade database to {revision_desc}...[/cyan]")
        _ = submission_db.upgrade_schema(revision=revision)
        console_err.print(f"[green]Successfully upgraded database to {revision_desc}![/green]")

    except (DatabaseConfigurationError, RuntimeError) as e:
        console_err.print(f"[red]Error during schema initialization: {e}[/red]")
        if isinstance(e, RuntimeError):
            console_err.print("[yellow]Ensure your database is running and accessible.[/yellow]")
            console_err.print(
                "[yellow]You might need to create an initial migration if this is the first time: 'alembic revision -m \"initial\" --autogenerate'[/yellow]"
            )
        raise click.ClickException(str(e)) from e
    except Exception as e:
        console_err.print(f"[red]An unexpected error occurred during 'db upgrade': {type(e).__name__} - {e}[/red]")
        raise click.ClickException(str(e)) from e


@db.command("list")
@grzcli.output_json
@limit
@click.option(
    "--state",
    "state_filters",
    type=click.Choice(SubmissionStateEnum.list(), case_sensitive=False),
    multiple=True,
    help="Filter by submission state. Can be passed multiple times.",
)
@click.option(
    "--filter-mode",
    type=click.Choice(SubmissionStateFilterModeEnum.list(), case_sensitive=False),
    default=SubmissionStateFilterModeEnum.LATEST.value,
    show_default=True,
    help="How --state is evaluated: 'latest' or 'any' state in history.",
)
@click.pass_context
def list_submissions(
    ctx: click.Context, output_json: bool, limit: int, state_filters: tuple[str, ...], filter_mode: str
):
    """Lists all submissions in the database with their latest state."""
    db = ctx.obj["db_url"]
    db_service = get_submission_db_instance(db)
    parsed_state_filters = tuple(SubmissionStateEnum(state) for state in state_filters) if state_filters else None
    parsed_filter_mode = SubmissionStateFilterModeEnum(filter_mode)

    try:
        submissions = db_service.list_submissions(
            limit=limit,
            state_filters=parsed_state_filters,
            state_filter_mode=parsed_filter_mode,
        )
    except Exception as e:
        raise click.ClickException(str(e)) from e

    if not submissions:
        console_err.print("[yellow]No submissions found in the database.[/yellow]")
        return

    table_title = "All Submissions" if not state_filters else f"Submissions ({', '.join(state_filters)})"
    table = rich.table.Table(title=table_title)
    table.add_column("ID", style="dim", min_width=29, width=29)
    table.add_column("tanG", style="cyan")
    table.add_column("Pseudonym", style="magenta")
    table.add_column("Latest State", style="green")
    table.add_column("Last State Timestamp (UTC)", style="yellow")
    table.add_column("Data Steward")
    table.add_column("Signature Status")

    submission_dicts = []

    for submission in submissions:
        latest_state_obj: SubmissionStateLog | None = None
        if submission.states:
            latest_state_obj = max(submission.states, key=lambda s: s.timestamp)

        latest_state_str = "N/A"
        latest_timestamp_str = "N/A"
        author_name_str = "N/A"
        signature_status = SignatureStatus.UNKNOWN
        verifying_key_comment = None

        if latest_state_obj:
            latest_state_str = latest_state_obj.state.value
            latest_state_str = (
                f"[red]{latest_state_str}[/red]" if latest_state_str == SubmissionStateEnum.ERROR else latest_state_str
            )
            latest_timestamp_str = latest_state_obj.timestamp.isoformat()
            author_name_str = latest_state_obj.author_name

            signature_status, verifying_key_comment = _verify_signature(
                ctx.obj["public_keys"], author_name_str, latest_state_obj
            )

        if output_json:
            submission_dict = _build_submission_dict_from(latest_state_obj, submission, signature_status)
            submission_dicts.append(submission_dict)
        else:
            table.add_row(
                submission.id,
                submission.tan_g[:8] + "…" if submission.tan_g is not None else _TEXT_MISSING,
                submission.pseudonym if submission.pseudonym is not None else _TEXT_MISSING,
                latest_state_str,
                latest_timestamp_str,
                author_name_str,
                signature_status.rich_display(verifying_key_comment),
            )

    if output_json:
        json.dump(submission_dicts, sys.stdout)
    else:
        console.print(table)


@db.command("list-change-requests")
@grzcli.output_json
@click.pass_context
def list_change_requests(ctx: click.Context, output_json: bool = False):
    """Lists all submissions in the database that have a change request."""
    db = ctx.obj["db_url"]
    db_service = get_submission_db_instance(db)
    submissions = db_service.list_change_requests()

    if not submissions:
        console_err.print("[yellow]No submissions found in the database.[/yellow]")
        return

    table = rich.table.Table(title="Submissions with change requests")
    table.add_column("ID", style="dim", width=12)
    table.add_column("tanG", style="cyan")
    table.add_column("Pseudonym", style="magenta")
    table.add_column("Change", style="green")
    table.add_column("Last State Timestamp (UTC)", style="yellow")
    table.add_column("Data Steward")
    table.add_column("Signature Status")

    submission_dicts = []

    for submission in submissions:
        for latest_change_request_obj in submission.changes:
            latest_change_str = "N/A"
            latest_timestamp_str = "N/A"
            author_name_str = "N/A"
            signature_status = SignatureStatus.UNKNOWN

            if latest_change_request_obj:
                latest_change_str = latest_change_request_obj.change.value
                latest_timestamp_str = latest_change_request_obj.timestamp.isoformat()
                author_name_str = latest_change_request_obj.author_name

                signature_status, verifying_key_comment = _verify_signature(
                    ctx.obj["public_keys"], author_name_str, latest_change_request_obj
                )

            if output_json:
                submission_dict = _build_submission_dict_from(latest_change_request_obj, submission, signature_status)
                submission_dicts.append(submission_dict)
            else:
                table.add_row(
                    submission.id,
                    submission.tan_g[:8] + "…" if submission.tan_g is not None else _TEXT_MISSING,
                    submission.pseudonym if submission.pseudonym is not None else _TEXT_MISSING,
                    latest_change_str,
                    latest_timestamp_str,
                    author_name_str,
                    signature_status.rich_display(verifying_key_comment),
                )

    if output_json:
        json.dump(submission_dicts, sys.stdout)
    else:
        console.print(table)


@db.command("tui")
@click.pass_context
def tui(ctx: click.Context):
    """Starts the interactive terminal user interface to the database."""
    db_url = ctx.obj["db_url"]
    public_keys = ctx.obj["public_keys"]
    database = get_submission_db_instance(db_url)

    # Prevent log messages from writing to stderr and messing up TUI. Since the
    # TUI is pretty much its own CLI context, it's fine to override the global
    # logging behavior here. TextualHandler() will make sure to still write log
    # messages visible to devtools.
    root_logger = logging.getLogger()
    for handler in root_logger.handlers:
        root_logger.removeHandler(handler)
    textual_handler = textual.logging.TextualHandler()
    # handlers define the format, so make sure Textual knows our project format
    textual_handler.setFormatter(logging.Formatter(fmt=LOGGING_FORMAT, datefmt=LOGGING_DATEFMT))
    root_logger.addHandler(textual_handler)

    app = DatabaseBrowser(database=database, public_keys=public_keys)
    app.run()


@db.command("should-qc")
@click.argument("submission_id")
@click.option(
    "--target-percentage",
    "target_percentage",
    type=click.FloatRange(0.0, 100.0),
    metavar="FLOAT",
    help="Minimum proportion of submissions that should be QCed (default = 2.0).",
    default=2.0,
)
@click.option(
    "--salt",
    "salt",
    help="Secret random string used as part of seed for random generator.",
    envvar="GRZCTL_SHOULD_QC_SALT",
)
@click.pass_context
def should_qc(ctx: click.Context, submission_id: str, target_percentage: float, salt: str | None):
    """Check whether a submission should be QCed."""
    database_url = ctx.obj["db_url"]
    database = get_submission_db_instance(database_url)

    try:
        result = database.should_qc(submission_id=submission_id, target_percentage=target_percentage, salt=salt)
        click.echo(str(result).lower())
    except SubmissionError as e:
        click.echo(f"Error: {e}", err=True)
        raise SystemExit(1) from e


def _build_submission_dict_from(
    log_obj: SubmissionStateLog | ChangeRequestLog | None,
    submission: Submission,
    signature_status: SignatureStatus,
) -> dict[str, Any]:
    """Serialize a submission and its latest log entry to a JSON-compatible dict.

    :param log_obj: The most recent :class:`~grz_db.models.submission.SubmissionStateLog` or
        :class:`~grz_db.models.submission.ChangeRequestLog`, or ``None`` if no log exists yet.
    :param submission: The submission ORM/Pydantic model instance.
    :param signature_status: Verification result for the log entry's author signature.
    :returns: A dictionary suitable for JSON serialisation that contains the submission identifiers
        and either a ``latest_state`` or ``latest_change_request`` key depending on *log_obj*.
    :raises TypeError: If *log_obj* is neither ``None`` nor one of the two expected log types.
    """
    submission_dict: dict[str, Any] = {
        "id": submission.id,
        "tan_g": submission.tan_g,
        "pseudonym": submission.pseudonym,
        "latest_state": None,
    }
    if log_obj:
        if isinstance(log_obj, SubmissionStateLog):
            submission_dict["latest_change_request"] = {}
            submission_dict["latest_state"] = {
                "timestamp": log_obj.timestamp.isoformat(),
                "data": log_obj.data,
                "data_steward": log_obj.author_name,
                "data_steward_signature": signature_status,
                "state": log_obj.state.value,
            }
        elif isinstance(log_obj, ChangeRequestLog):
            submission_dict["latest_state"] = {}
            submission_dict["latest_change_request"] = {
                "timestamp": log_obj.timestamp.isoformat(),
                "data": log_obj.data,
                "data_steward": log_obj.author_name,
                "data_steward_signature": signature_status,
                "change": log_obj.change.value,
            }
        else:
            raise TypeError(f"unknown type {type(log_obj)}")
    return submission_dict


@submission.command()
@click.argument("submission_id", type=str)
@click.pass_context
def add(ctx: click.Context, submission_id: str):
    """
    Add a submission to the database.
    """
    db = ctx.obj["db_url"]
    db_service = get_submission_db_instance(db)
    try:
        db_submission = db_service.add_submission(submission_id)
        console_err.print(f"[green]Submission '{db_submission.id}' added successfully.[/green]")
    except SubmissionError as e:
        console_err.print(f"[red]Error: {e}[/red]")
        raise click.Abort() from e
    except Exception as e:
        console_err.print(f"[red]An unexpected error occurred: {e}[/red]")
        raise click.ClickException(f"Failed to add submission: {e}") from e


@submission.command()
@click.argument("submission_id", type=str)
@click.argument("state_str", metavar="STATE", type=click.Choice(SubmissionStateEnum.list(), case_sensitive=False))
@click.option("--data", "data_json", type=str, default=None, help='Additional JSON data (e.g., \'{"k":"v"}\').')
@click.option("--ignore-error-state/--confirm-error-state")
@click.pass_context
def update(ctx: click.Context, submission_id: str, state_str: str, data_json: str | None, ignore_error_state: bool):  # noqa: C901
    """Update a submission to the given state. Optionally accepts additional JSON data to associate with the log entry."""
    db = ctx.obj["db_url"]
    db_service = get_submission_db_instance(db, author=ctx.obj["author"])
    try:
        state_enum = SubmissionStateEnum(state_str)
    except ValueError as e:
        console_err.print(f"[red]Error: Invalid state value '{state_str}'.[/red]")
        raise click.Abort() from e

    parsed_data = None
    if data_json:
        try:
            parsed_data = json.loads(data_json)
        except json.JSONDecodeError as e:
            console_err.print(f"[red]Error: Invalid JSON string for --data: {data_json}[/red]")
            raise click.Abort() from e
    try:
        submission = db_service.get_submission(submission_id)
        if not submission:
            raise SubmissionNotFoundError(submission_id)
        latest_state = submission.get_latest_state()
        latest_state_is_error = latest_state is not None and latest_state.state == SubmissionStateEnum.ERROR
        if (
            latest_state_is_error
            and not ignore_error_state
            and not click.confirm(
                f"Submission is currently in an 'Error' state. Are you sure you want to set it to '{state_enum}'?",
                default=False,
                show_default=True,
            )
        ):
            console_err.print(f"[yellow]Not modifying state of errored submission '{submission_id}'.[/yellow]")
            ctx.exit()

        new_state_log = db_service.update_submission_state(submission_id, state_enum, parsed_data)
        console_err.print(
            f"[green]Submission '{submission_id}' updated to state '{new_state_log.state.value}'. Log ID: {new_state_log.id}[/green]"
        )
        if new_state_log.data:
            console_err.print(f"  Data: {new_state_log.data}")

    except SubmissionNotFoundError as e:
        console_err.print(f"[red]Error: {e}[/red]")
        console_err.print(f"You might need to add it first: grz-cli db submission add {submission_id}")
        raise click.Abort() from e
    except click.exceptions.Exit as e:
        if e.exit_code != 0:
            raise e
    except Exception as e:
        console_err.print(f"[red]An unexpected error occurred: {e}[/red]")
        traceback.print_exc()
        raise click.ClickException(f"Failed to update submission state: {e}") from e


@submission.command(
    epilog="Currently available KEYs are: "
    + ", ".join(sorted(Submission.model_fields.keys() - Submission.immutable_fields))
)
@click.argument("submission_id", type=str)
@click.argument("key", metavar="KEY", type=click.Choice(Submission.model_fields.keys()))
@click.argument("value", metavar="VALUE", type=str)
@click.pass_context
def modify(ctx: click.Context, submission_id: str, key: str, value: str):
    """
    Modify a submission's database properties.
    """
    db = ctx.obj["db_url"]
    db_service = get_submission_db_instance(db, author=ctx.obj["author"])

    try:
        submission = db_service.get_submission(submission_id)
        if not submission:
            raise SubmissionNotFoundError(submission_id)
        _ = db_service.modify_submission(submission_id, key, value)
        console_err.print(f"[green]Updated {key} of submission '{submission_id}'[/green]")
    except SubmissionNotFoundError as e:
        console_err.print(f"[red]Error: {e}[/red]")
        console_err.print(f"You might need to add it first: grz-cli db submission add {submission_id}")
        raise click.Abort() from e
    except Exception as e:
        console_err.print(f"[red]An unexpected error occurred: {e}[/red]")
        traceback.print_exc()
        raise click.ClickException(f"Failed to update submission state: {e}") from e


def _prepare_submission_console_table(submission_diff: "SubmissionDiffCollection") -> rich.console.RenderableType:
    """Build a Rich renderable that shows pending submission-level metadata changes.

    :param submission_diff: :class:`SubmissionDiff` instance produced by :func:`diff_metadata`.
    :returns: A :class:`rich.table.Table` when there are pending changes, or a plain text message otherwise.
    """
    pending = [d for d in submission_diff.pending if d.key != "submission_metadata"]
    if pending:
        diff_table_tbl = rich.table.Table(title="Submission Metadata")
        diff_table_tbl.add_column("Key")
        diff_table_tbl.add_column("Before")
        diff_table_tbl.add_column("After")
        for field_diff in sorted(pending, key=lambda d: d.key):
            diff_table_tbl.add_row(
                field_diff.key,
                str(field_diff.diff.before) if field_diff.diff.before is not None else _TEXT_MISSING,
                str(field_diff.diff.after),
            )
        diff_table: rich.console.RenderableType = diff_table_tbl
    else:
        diff_table = rich.padding.Padding(rich.text.Text("No changes to submission-level metadata."), pad=(0, 0, 0, 0))
    return diff_table


def _prepare_donor_console_table(
    donor_data: list[FieldDiff], donor_id: str, status: str
) -> rich.console.RenderableType:
    """Build a Rich renderable that shows pending changes for a single donor.

    :param donor_data: List of :class:`FieldDiff` instances for the donor's fields.
    :param donor_id: Pseudonym of the donor (used in the table title).
    :param status: Human-readable database status string (e.g. ``"new"`` or ``"update"``).
    :returns: A :class:`rich.table.Table` listing only the fields whose value changed.
    """
    table_title = f"[green]Donor '{donor_id}' database status: {status}[/green]"
    diff_table = rich.table.Table(title=table_title, min_width=len(table_title), title_justify="left")
    diff_table.add_column("Key")
    diff_table.add_column("Before")
    diff_table.add_column("After")
    for field_diff in sorted(donor_data, key=lambda d: d.key):
        if field_diff.diff.before != field_diff.diff.after:
            diff_table.add_row(
                field_diff.key,
                _TEXT_MISSING if field_diff.diff.before is None else rich.pretty.Pretty(field_diff.diff.before),
                rich.pretty.Pretty(field_diff.diff.after),
            )
    return diff_table


@submission.command()
@click.argument("submission_id", type=str)
@click.argument("metadata_path", metavar="path/to/metadata.json", type=str)
@click.option(
    "--submission_date",
    type=click.DateTime(formats=["%Y-%m-%d"]),
    default=None,
    help="Submission date of the submission; overwrites submissionDate in metadata.json",
)
@click.option(
    "--confirm/--no-confirm",
    default=True,
    help="Whether to confirm changes before committing to database. (Default: confirm)",
)
@click.option(
    "--ignore-field",
    help="Do not populate the given key from the metadata to the database. Can be specified multiple times to ignore multiple keys.",
    multiple=True,
)
@click.pass_context
def populate(  # noqa: C901, PLR0913
    ctx: click.Context,
    submission_id: str,
    metadata_path: str,
    submission_date: datetime | None,
    confirm: bool,
    ignore_field: tuple[str, ...],
):
    """Populate a submission in the database based on the given metadata.json file."""
    log.debug("Ignored fields for populate: %s", ignore_field)

    if submission_date is not None:
        log.info("Submission date from provided option is used")
        if submission_date.date() >= date.today() + timedelta(days=1):
            raise RuntimeError(
                f"Submission date ({submission_date.date()}) is set to a future date (today: {date.today()}) which is not allowed"
            )
    else:
        log.warning("Submission date from metadata.json is used")

    db = ctx.obj["db_url"]
    db_service = get_submission_db_instance(db, author=ctx.obj["author"])

    try:
        submission = db_service.get_submission(submission_id)
        if not submission:
            raise SubmissionNotFoundError(submission_id)
    except SubmissionNotFoundError as e:
        console_err.print(f"[red]Error: {e}[/red]")
        console_err.print(f"You might need to add it first: grz-cli db submission add {submission_id}")
        raise click.Abort() from e
    except Exception as e:
        console_err.print(f"[red]An unexpected error occurred: {e}[/red]")
        traceback.print_exc()
        raise click.ClickException(f"Failed to update submission state: {e}") from e

    with open(metadata_path) as fd:
        metadata = GrzSubmissionMetadata.model_validate_json(fd.read())

    if metadata.submission.tan_g == REDACTED_TAN and "tan_g" not in ignore_field:
        raise ValueError(
            f"Submission {submission_id} has redacted tan_g in metadata.json: {metadata_path}. "
            "Refusing to populate a seemingly-redacted TAN. "
            "Add 'tan_g' to ignore fields and/or use 'grzctl db submission modify' directly."
        )
    if (
        not metadata.submission.local_case_id or metadata.submission.local_case_id == REDACTED_LOCAL_CASE_ID
    ) and "pseudonym" not in ignore_field:
        raise ValueError(
            f"Submission {submission_id} has missing or redacted local_case_id in metadata.json: {metadata_path}. "
            "Add 'pseudonym' to ignore fields and/or use 'grzctl db submission modify' directly."
        )

    submission_diff, donors_diff = db_service.diff(
        submission_id,
        metadata,
        submission_date,
        ignore_fields=set(ignore_field),
    )

    # build donor diff and attach Rich tables for console preview in one pass
    diff_tables: list[rich.console.RenderableType] = []
    for donor_diff in donors_diff.added + donors_diff.updated:
        diff_tables.append(
            _prepare_donor_console_table(donor_diff.changes, donor_diff.pseudonym or "", donor_diff.state)
        )
    for donor_diff in donors_diff.deleted:
        diff_tables.append(rich.text.Text(f"Donor {donor_diff.pseudonym} deleted", style="red"))

    if not submission_diff.has_pending and not donors_diff.has_pending:
        console_err.print("[green]Database is already up to date with the provided metadata![/green]")
        ctx.exit()
        return  # ctx.exit() raises SystemExit; this line satisfies static analysers

    console.print(
        rich.panel.Panel.fit(
            rich.console.Group(_prepare_submission_console_table(submission_diff), *diff_tables, fit=True),
            title="Pending Changes",
        )
    )

    if not confirm or click.confirm(
        "Are you sure you want to commit these changes to the database?",
        default=False,
        show_default=True,
    ):
        db_service.commit_changes(submission_id, submission_diff, donors_diff)
        console_err.print("[green]Database populated successfully.[/green]")


class QCStatus(StrEnum):
    PASS = "PASS"  # noqa: S105
    FAIL = "FAIL"
    TOO_LOW = "TOO LOW"
    THRESHOLD_NOT_MET = "THRESHOLD NOT MET"


class QCReportRow(StrictBaseModel):
    """Pydantic model representing a single row from a detailed QC pipeline report CSV."""

    sample_id: str
    donor_pseudonym: str
    lab_data_name: str
    library_type: LibraryType
    sequence_subtype: SequenceSubtype
    genomic_study_subtype: GenomicStudySubtype
    quality_control_status: QCStatus
    mean_depth_of_coverage: float
    mean_depth_of_coverage_provided: float
    mean_depth_of_coverage_required: float
    mean_depth_of_coverage_deviation: float
    mean_depth_of_coverage_qc_status: QCStatus = Field(alias="meanDepthOfCoverageQCStatus")
    percent_bases_above_quality_threshold: float
    quality_threshold: float
    percent_bases_above_quality_threshold_provided: float
    percent_bases_above_quality_threshold_required: float
    percent_bases_above_quality_threshold_deviation: float
    percent_bases_above_quality_threshold_qc_status: QCStatus = Field(alias="percentBasesAboveQualityThresholdQCStatus")
    targeted_regions_above_min_coverage: float
    min_coverage: float
    targeted_regions_above_min_coverage_provided: float
    targeted_regions_above_min_coverage_required: float
    targeted_regions_above_min_coverage_deviation: float
    targeted_regions_above_min_coverage_qc_status: QCStatus = Field(alias="targetedRegionsAboveMinCoverageQCStatus")


@submission.command()
@click.argument("submission_id", type=str)
@click.argument("report_csv_path", metavar="path/to/report.csv", type=grzcli.FILE_R_E)
@click.option(
    "--confirm/--no-confirm",
    default=True,
    help="Whether to confirm changes before committing to database. (Default: confirm)",
)
@click.pass_context
def populate_qc(ctx: click.Context, submission_id: str, report_csv_path: str, confirm: bool):
    """Populate the submission database from a detailed QC pipeline report."""
    db = ctx.obj["db_url"]
    db_service = get_submission_db_instance(db, author=ctx.obj["author"])

    with open(report_csv_path, encoding="utf-8", newline="") as report_csv_file:
        reader = csv.reader(report_csv_file)
        header = next(reader)
        reports = []
        for row in reader:
            reports.append(QCReportRow(**dict(zip(header, row, strict=True))))

    report_mtime = datetime.fromtimestamp(Path(report_csv_path).stat().st_mtime, tz=UTC)
    results = []
    for report in reports:
        results.append(
            DetailedQCResult(
                submission_id=submission_id,
                lab_datum_id=report.sample_id,
                pseudonym=report.donor_pseudonym,
                timestamp=report_mtime,
                sequence_type=SequenceType.dna,  # pipeline only supports DNA and doesn't pass type to report.csv
                sequence_subtype=report.sequence_subtype,
                library_type=report.library_type,
                percent_bases_above_quality_threshold_minimum_quality=report.quality_threshold,
                percent_bases_above_quality_threshold_percent=report.percent_bases_above_quality_threshold,
                percent_bases_above_quality_threshold_passed_qc=report.percent_bases_above_quality_threshold_qc_status
                == QCStatus.PASS,
                percent_bases_above_quality_threshold_percent_deviation=report.percent_bases_above_quality_threshold_deviation,
                mean_depth_of_coverage=report.mean_depth_of_coverage,
                mean_depth_of_coverage_passed_qc=report.mean_depth_of_coverage_qc_status == QCStatus.PASS,
                mean_depth_of_coverage_percent_deviation=report.mean_depth_of_coverage_deviation,
                targeted_regions_min_coverage=report.min_coverage,
                targeted_regions_above_min_coverage=report.targeted_regions_above_min_coverage,
                targeted_regions_above_min_coverage_passed_qc=report.targeted_regions_above_min_coverage_qc_status
                == QCStatus.PASS,
                targeted_regions_above_min_coverage_percent_deviation=report.targeted_regions_above_min_coverage_deviation,
            )
        )
    table = rich.table.Table(
        "Submission ID",
        "Lab Datum ID",
        "Pseudonym",
        "Timestamp",
        "Sequence Type",
        "Sequence Subtype",
        "Library Type",
        "PBaQT",
        "MDoC",
        "TRaMC",
        title="New Detailed QC Results",
    )
    for result in results:
        table.add_row(
            result.submission_id,
            result.lab_datum_id,
            result.pseudonym,
            f"{result.timestamp:%c}",
            result.sequence_type,
            result.sequence_subtype,
            result.library_type,
            rich.pretty.Pretty(result.percent_bases_above_quality_threshold_percent),
            rich.pretty.Pretty(result.mean_depth_of_coverage),
            rich.pretty.Pretty(result.targeted_regions_above_min_coverage),
        )
    console.print(table)

    if not confirm or click.confirm(
        "Are you sure you want to commit these changes to the database?", default=False, show_default=True
    ):
        for result in results:
            db_service.add_detailed_qc_result(result)


@submission.command()
@click.argument("submission_id", type=str)
@click.argument("change_str", metavar="CHANGE", type=click.Choice(ChangeRequestEnum.list(), case_sensitive=False))
@click.option("--data", "data_json", type=str, default=None, help='Additional JSON data (e.g., \'{"k":"v"}\').')
@click.pass_context
def change_request(ctx: click.Context, submission_id: str, change_str: str, data_json: str | None):
    """Register a completed change request for the given submission. Optionally accepts additional JSON data to associate with the log entry."""
    db = ctx.obj["db_url"]
    db_service = get_submission_db_instance(db, author=ctx.obj["author"])
    try:
        change_request_enum = ChangeRequestEnum(change_str)
    except ValueError as e:
        console_err.print(f"[red]Error: Invalid change request value '{change_str}'.[/red]")
        raise click.Abort() from e

    parsed_data = None
    if data_json:
        try:
            parsed_data = json.loads(data_json)
        except json.JSONDecodeError as e:
            console_err.print(f"[red]Error: Invalid JSON string for --data: {data_json}[/red]")
            raise click.Abort() from e
    try:
        new_change_request_log = db_service.add_change_request(submission_id, change_request_enum, parsed_data)
        console_err.print(
            f"[green]Submission '{submission_id}' has undergone a change request of '{new_change_request_log.change.value}'. Log ID: {new_change_request_log.id}[/green]"
        )
        if new_change_request_log.data:
            console_err.print(f"  Data: {new_change_request_log.data}")

    except SubmissionNotFoundError as e:
        console_err.print(f"[red]Error: {e}[/red]")
        console_err.print(f"You might need to add it first: grz-cli db submission add {submission_id}")
        raise click.Abort() from e
    except Exception as e:
        console_err.print(f"[red]An unexpected error occurred: {e}[/red]")
        traceback.print_exc()
        raise click.ClickException(f"Failed to update submission state: {e}") from e


@submission.command("show")
@click.argument("submission_id", type=str)
@output_json
@click.pass_context
def show(ctx: click.Context, submission_id: str, output_json: bool):
    """
    Show details of a submission.
    """
    db = ctx.obj["db_url"]
    db_service = get_submission_db_instance(db)
    submission = db_service.get_submission(submission_id)
    if not submission:
        console_err.print(f"[red]Error: Submission with ID '{submission_id}' not found.[/red]")
        raise click.Abort()

    if output_json:
        submission_dict = submission.model_dump(mode="json")
        submission_dict["states"] = []

        for state_log in sorted(submission.states, key=lambda s: s.timestamp):
            signature_status, verifying_key_comment = _verify_signature(
                ctx.obj["public_keys"], state_log.author_name, state_log
            )
            state_dict = state_log.model_dump(mode="json", include={"id", "timestamp", "state", "data"})
            state_dict["data_steward"] = state_log.author_name
            state_dict["data_steward_signature"] = signature_status
            state_dict["signature_key_comment"] = verifying_key_comment
            submission_dict["states"].append(state_dict)

        json.dump(submission_dict, sys.stdout)
        sys.stdout.write("\n")
        return

    attribute_table = rich.table.Table(box=None)
    attribute_table.add_column("Attribute", justify="right")
    attribute_table.add_column("Value")
    for label, attr_name in (
        ("tanG", "tan_g"),
        ("Pseudonym", "pseudonym"),
        ("Submission Date", "submission_date"),
        ("Submission Size", "submission_size"),
        ("Submission Type", "submission_type"),
        ("Submitter ID", "submitter_id"),
        ("Data Node ID", "data_node_id"),
        ("Disease Type", "disease_type"),
        ("Genomic Study Type", "genomic_study_type"),
        ("Genomic Study Subtype", "genomic_study_subtype"),
        ("Basic QC Passed", "basic_qc_passed"),
        ("Consented", "consented"),
        ("Selected For QC", "selected_for_qc"),
        ("Detailed QC Passed", "detailed_qc_passed"),
    ):
        attr = getattr(submission, attr_name)
        attribute_table.add_row(
            rich.text.Text(f"{label}", style="cyan"), rich.text.Text(str(attr)) if attr is not None else _TEXT_MISSING
        )

    renderables: list[rich.console.RenderableType] = [rich.padding.Padding(attribute_table, (1, 0))]
    if submission.states:
        state_table = rich.table.Table(title="State History")
        state_table.add_column("Log ID", style="dim", width=12)
        state_table.add_column("Timestamp (UTC)", style="yellow")
        state_table.add_column("State", style="green")
        state_table.add_column("Data", style="cyan", overflow="ellipsis")
        state_table.add_column("Data Steward", style="magenta")
        state_table.add_column("Signature Status")

        sorted_states = sorted(submission.states, key=lambda s: s.timestamp)
        for state_log in sorted_states:
            data_str = json.dumps(state_log.data) if state_log.data else ""
            state = state_log.state.value
            state_str = f"[red]{state}[/red]" if state == SubmissionStateEnum.ERROR else state
            data_steward_str = state_log.author_name
            signature_status, verifying_key_comment = _verify_signature(
                ctx.obj["public_keys"], data_steward_str, state_log
            )
            signature_status_str = signature_status.rich_display(verifying_key_comment)

            state_table.add_row(
                str(state_log.id),
                state_log.timestamp.isoformat(),
                state_str,
                data_str,
                data_steward_str,
                signature_status_str,
            )
        renderables.append(state_table)
    else:
        renderables.append(rich.text.Text("No state history found for this submission.", style="yellow"))

    panel = rich.panel.Panel.fit(
        rich.console.Group(*renderables),
        title=f"Submission {submission.id}",
    )
    console.print(panel)


@db.command("sync-from-inbox")
@grzcli.configuration
@click.pass_context
def sync_from_inbox(
    ctx: click.Context,
    configuration: dict[str, Any],
    **kwargs,
):
    """
    Synchronize the database with submissions found in the inbox.
    """
    try:
        list_config = ListConfig.model_validate(configuration)
    except Exception:
        console_err.print(f"[red]Error loading S3 configuration: {traceback.format_exc()}[/red]")
        sys.exit(1)

    db_url = ctx.obj["db_url"]
    author = ctx.obj["author"]
    db_service = get_submission_db_instance(db_url, author=author)

    try:
        console_err.print(f"[cyan]Scanning inbox '{list_config.s3.bucket}'...[/cyan]")
        s3_submissions = query_submissions(list_config.s3, show_cleaned=False)

        console_err.print(f"[cyan]Synchronizing {len(s3_submissions)} submissions with database...[/cyan]")
        sync_submissions(db_service, s3_submissions, author)

        console_err.print("[green]Synchronization complete.[/green]")

    except Exception:
        console_err.print(f"[red]Error during synchronization: {traceback.format_exc()}[/red]")
        traceback.print_exc()
        sys.exit(1)
