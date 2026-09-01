"""Command for managing a submission database"""

import csv
import json
import logging
import sys
import traceback
from collections import Counter, namedtuple
from collections.abc import Iterable
from datetime import UTC, date, datetime, timedelta
from enum import StrEnum
from pathlib import Path
from typing import Any, NamedTuple, NoReturn

import botocore.exceptions
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
from grz_common.transfer import init_s3_client
from grz_common.workers.download import query_submissions
from grz_db.errors import (
    CaseHasLinkedSubmissionsError,
    CaseNotFoundError,
    DatabaseConfigurationError,
    DuplicateCaseError,
    DuplicatePsnError,
    SubmissionError,
    SubmissionNotFoundError,
    SubmissionTypeInvalidForCaseError,
)
from grz_db.models.author import Author
from grz_db.models.submission import (
    Case,
    ChangeRequestEnum,
    ChangeRequestLog,
    DetailedQCResult,
    DiffState,
    FailureReasonEnum,
    FieldDiff,
    Submission,
    SubmissionBase,
    SubmissionChangeSet,
    SubmissionDb,
    SubmissionStateEnum,
    SubmissionStateFilterModeEnum,
    SubmissionStateLog,
)
from grz_pydantic_models.common import StrictBaseModel
from grz_pydantic_models.submission.metadata import (
    GenomicStudySubtype,
    GrzSubmissionMetadata,
    LibraryType,
    SequenceSubtype,
    SequenceType,
)
from pydantic import Field, ValidationError
from tqdm.auto import tqdm

from ... import get_versions
from ...commands import grzctl_configuration, inbox_options
from ...models.config import GrzctlConfig
from .. import limit
from ..change_request import resolve_and_validate_change_request
from . import SignatureStatus, _verify_signature
from .sync import sync_submissions
from .tui import DatabaseBrowser

console = rich.console.Console()
console_err = rich.console.Console(stderr=True)
log = logging.getLogger(__name__)
_TEXT_MISSING = rich.text.Text("missing", style="italic yellow")


def _attribute_table(rows: Iterable[tuple[str, Any]]) -> rich.table.Table:
    """Build the two-column Attribute/Value table the ``show`` commands render.

    A ``None`` renders as *missing* rather than as the string "None", so a column that was
    never populated cannot be mistaken for one holding that text.

    :param rows: ``(label, value)`` pairs, in display order.
    """
    table = rich.table.Table(box=None)
    table.add_column("Attribute", justify="right")
    table.add_column("Value")
    for label, value in rows:
        table.add_row(
            rich.text.Text(label, style="cyan"),
            rich.text.Text(str(value)) if value is not None else _TEXT_MISSING,
        )
    return table


def _dump_json(payload: Any) -> None:
    """Write *payload* to stdout as JSON, indented, and nothing else.

    Not through the rich console: it decides on ANSI colour when it is constructed, which for
    this module happens at import, and it honours ``FORCE_COLOR`` even when stdout is a pipe.
    A consumer then parses escape codes. ``TTY_COMPATIBLE=0`` suppresses that, but only the
    test suite sets it. Indented because a case is read by eye as often as by a program.
    """
    json.dump(payload, sys.stdout, indent=2)
    sys.stdout.write("\n")


def _abort(e: Exception) -> NoReturn:
    """Print *e* as a red error on stderr and abort the command."""
    console_err.print(f"[red]Error: {e}[/red]")
    raise click.Abort() from e


def _abort_missing_submission(e: SubmissionNotFoundError, submission_id: str) -> NoReturn:
    """Like :func:`_abort`, plus a hint on how to add the missing submission."""
    console_err.print(f"[red]Error: {e}[/red]")
    console_err.print(f"You might need to add it first: grzctl db submission add {submission_id}")
    raise click.Abort() from e


_submission_id_argument = click.argument("submission_id", type=str)


def get_submission_db_instance(db_url: str, author: Author | None = None) -> SubmissionDb:
    """Creates and returns an instance of SubmissionDb."""
    return SubmissionDb(db_url=db_url, author=author)


@click.group(help="Database operations")
@grzctl_configuration
@click.pass_context
def db(
    ctx: click.Context,
    configuration: GrzctlConfig,
    **kwargs,
):
    """Database operations"""
    # set up context object
    ctx.ensure_object(dict)

    db_config = configuration.db
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


@db.group()
@click.pass_context
def case(ctx: click.Context):
    """Case operations"""
    pass


_CASE_MUTABLE_KEYS = sorted(Case.mutable_fields())

_case_id_argument = click.argument("case_id", type=int)


@case.command("list")
@output_json
@click.pass_context
def case_list(ctx: click.Context, output_json: bool):
    """List all cases with their linked-submission counts."""
    db = ctx.obj["db_url"]
    db_service = get_submission_db_instance(db, author=ctx.obj["author"])
    cases = db_service.list_cases()

    if output_json:
        payload = [{**case_obj.model_dump(mode="json"), "submission_count": count} for case_obj, count in cases]
        _dump_json(payload)
        return

    table = rich.table.Table(title="Cases")
    table.add_column("ID", justify="right")
    table.add_column("PSN")
    table.add_column("Submitter ID")
    table.add_column("Local Case ID")
    table.add_column("Submissions", justify="right")
    for case_obj, count in cases:
        table.add_row(
            str(case_obj.id),
            case_obj.psn if case_obj.psn is not None else "-",
            case_obj.submitter_id if case_obj.submitter_id is not None else "-",
            case_obj.local_case_id if case_obj.local_case_id is not None else "-",
            str(count),
        )
    console.print(table)


@case.command("list-unlinked")
@output_json
@click.pass_context
def case_list_unlinked(ctx: click.Context, output_json: bool):
    """List submissions that carry no case, and why each carries none.

    Test submissions are excluded, since they are never case-tracked. A submission whose type
    is not known yet, or that is simply not linked yet, resolves itself on the next populate or
    backfill. The other two reasons need an operator: see 'db case list-ambiguous-keys'.
    """
    db = ctx.obj["db_url"]
    db_service = get_submission_db_instance(db, author=ctx.obj["author"])
    unlinked = db_service.list_unlinked_submissions()

    if output_json:
        _dump_json(
            [
                {
                    "submission_id": row.submission_id,
                    "submitter_id": row.submitter_id,
                    "local_case_id": row.local_case_id,
                    "submission_type": row.submission_type,
                    "reason": row.reason,
                }
                for row in unlinked
            ]
        )
        return

    if not unlinked:
        console_err.print("[green]Every case-trackable submission is linked to a case.[/green]")
        return

    table = rich.table.Table(title="Unlinked submissions")
    table.add_column("Submission ID")
    table.add_column("Type")
    table.add_column("Submitter ID")
    table.add_column("Local Case ID")
    table.add_column("Reason")
    for row in unlinked:
        table.add_row(
            row.submission_id,
            str(row.submission_type) if row.submission_type is not None else _TEXT_MISSING,
            row.submitter_id if row.submitter_id is not None else _TEXT_MISSING,
            row.local_case_id if row.local_case_id is not None else _TEXT_MISSING,
            str(row.reason).replace("_", " "),
        )
    console.print(table)

    needs_operator = sum(1 for row in unlinked if row.reason.needs_operator)
    console_err.print(f"[cyan]{len(unlinked)} unlinked, {needs_operator} of which will not resolve on a re-run.[/cyan]")


@case.command("list-ambiguous-keys")
@output_json
@click.pass_context
def case_list_ambiguous_keys(ctx: click.Context, output_json: bool):
    """List submitter-local case keys that name more than one patient.

    A patient has one initial submission that passes basic QC, so a key carrying several was
    reused across patients. Such a key opens no case and its submissions keep case_id NULL.
    This is the list the cases migration logs once during 'db upgrade'.
    """
    db = ctx.obj["db_url"]
    db_service = get_submission_db_instance(db, author=ctx.obj["author"])
    keys = db_service.list_ambiguous_case_keys()

    if output_json:
        _dump_json(
            [
                {
                    "submitter_id": key.submitter_id,
                    "local_case_id": key.local_case_id,
                    "qc_passed_initials": key.qc_passed_initials,
                    "submissions": key.submissions,
                }
                for key in keys
            ]
        )
        return

    if not keys:
        console_err.print("[green]Every submitter-local case key names a single patient.[/green]")
        return

    table = rich.table.Table(title="Keys naming more than one patient")
    table.add_column("Submitter ID")
    table.add_column("Local Case ID")
    table.add_column("QC-passed initials", justify="right")
    table.add_column("Submissions", justify="right")
    for key in keys:
        table.add_row(
            key.submitter_id,
            key.local_case_id,
            str(key.qc_passed_initials),
            str(key.submissions),
        )
    console.print(table)
    console_err.print(
        "[cyan]Each key stands for as many patients as it has QC-passed initial submissions. "
        "Create a case per patient with 'db case create', then link with 'db case relink'.[/cyan]"
    )


@case.command("show")
@_case_id_argument
@output_json
@click.pass_context
def case_show(ctx: click.Context, case_id: int, output_json: bool):
    """Show a case and its linked submissions."""
    db = ctx.obj["db_url"]
    db_service = get_submission_db_instance(db, author=ctx.obj["author"])
    case_obj = db_service.get_case(case_id)
    if case_obj is None:
        _abort(CaseNotFoundError(case_id))
    submissions = db_service.list_submissions_for_case(case_id)

    if output_json:
        payload = {
            **case_obj.model_dump(mode="json"),
            "submissions": [
                {
                    "id": s.id,
                    "submission_type": s.submission_type,
                    "basic_qc_passed": s.basic_qc_passed,
                    "latest_state": (latest.state if (latest := s.get_latest_state()) else None),
                }
                for s in submissions
            ],
        }
        _dump_json(payload)
        return

    attribute_table = _attribute_table(
        (
            ("ID", case_obj.id),
            ("PSN", case_obj.psn),
            ("Submitter ID", case_obj.submitter_id),
            ("Local Case ID", case_obj.local_case_id),
        )
    )

    submissions_table = rich.table.Table(title="Submissions")
    submissions_table.add_column("Submission ID")
    submissions_table.add_column("Type")
    # Shows which initial submission holds the case's QC-passed slot.
    submissions_table.add_column("Basic QC")
    submissions_table.add_column("Latest State")
    for s in submissions:
        latest = s.get_latest_state()
        submissions_table.add_row(
            s.id,
            str(s.submission_type) if s.submission_type is not None else "-",
            "pending" if s.basic_qc_passed is None else ("passed" if s.basic_qc_passed else "failed"),
            str(latest.state) if latest is not None else "-",
        )
    console.print(rich.panel.Panel.fit(rich.console.Group(attribute_table, submissions_table), title=f"Case {case_id}"))


@case.command("modify", epilog="Currently available KEYs are: " + ", ".join(_CASE_MUTABLE_KEYS))
@_case_id_argument
@click.argument("key", metavar="KEY", type=click.Choice(_CASE_MUTABLE_KEYS))
@click.argument("value", metavar="VALUE", type=str)
@click.pass_context
def case_modify(ctx: click.Context, case_id: int, key: str, value: str):
    """Modify a mutable field of a case."""
    db = ctx.obj["db_url"]
    db_service = get_submission_db_instance(db, author=ctx.obj["author"])
    try:
        db_service.modify_case(case_id, key, value)
        console_err.print(f"[green]Updated {key} of case {case_id}[/green]")
    except (CaseNotFoundError, DuplicateCaseError, DuplicatePsnError, ValueError) as e:
        _abort(e)


@case.command("create")
@click.argument("submitter_id", type=str)
@click.argument("local_case_id", type=str)
@click.option("--psn", default=None, help="RKI pseudonym (must be unique).")
@click.pass_context
def case_create(ctx: click.Context, submitter_id: str, local_case_id: str, psn: str | None):
    """Manually create a case (cases are normally created during populate)."""
    db = ctx.obj["db_url"]
    db_service = get_submission_db_instance(db, author=ctx.obj["author"])
    try:
        case_obj = db_service.create_case(submitter_id=submitter_id, local_case_id=local_case_id, psn=psn)
        console_err.print(
            f"[green]Created case {case_obj.id} for submitter '{submitter_id}', local case '{local_case_id}'[/green]"
        )
    except (DuplicateCaseError, DuplicatePsnError, ValueError) as e:
        _abort(e)


@case.command("relink")
@_submission_id_argument
@_case_id_argument
@click.pass_context
def case_relink(ctx: click.Context, submission_id: str, case_id: int):
    """Relink a submission to a different case for repair.

    A case may still have at most one initial submission that passed basic QC. The case
    being vacated is left as-is and may become empty.
    """
    db = ctx.obj["db_url"]
    db_service = get_submission_db_instance(db, author=ctx.obj["author"])
    try:
        db_service.set_submission_case(submission_id, case_id)
        console_err.print(f"[green]Linked submission '{submission_id}' to case {case_id}[/green]")
    # DuplicateInitialSubmissionError subclasses SubmissionTypeInvalidForCaseError, so the case's
    # one-initial rule is caught here too.
    except (SubmissionNotFoundError, CaseNotFoundError, SubmissionTypeInvalidForCaseError) as e:
        _abort(e)


@case.command("unlink")
@_submission_id_argument
@click.pass_context
def case_unlink(ctx: click.Context, submission_id: str):
    """Unlink a submission from its case, leaving it linked to none.

    The counterpart to `relink`, which can only move a submission between cases. The case
    being vacated is left as-is and may become empty. Unlinking an already unlinked
    submission changes nothing.
    """
    db = ctx.obj["db_url"]
    db_service = get_submission_db_instance(db, author=ctx.obj["author"])
    try:
        db_service.clear_submission_case(submission_id)
        console_err.print(f"[green]Unlinked submission '{submission_id}' from its case[/green]")
    except SubmissionNotFoundError as e:
        _abort(e)


@case.command("delete")
@_case_id_argument
@click.pass_context
def case_delete(ctx: click.Context, case_id: int):
    """Delete an empty case (refuses when submissions are still linked)."""
    db = ctx.obj["db_url"]
    db_service = get_submission_db_instance(db, author=ctx.obj["author"])
    try:
        db_service.delete_case(case_id)
        console_err.print(f"[green]Deleted case {case_id}[/green]")
    except (CaseNotFoundError, CaseHasLinkedSubmissionsError) as e:
        _abort(e)


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
        sys.stdout.write("\n")
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
        sys.stdout.write("\n")
    else:
        console.print(table)


@db.command("tui")
@click.pass_context
@click.option(
    "--quarter",
    type=click.IntRange(min=1, max=4),
    default=None,
    help="Quarter (1-4) for the 'Detailed QC by LE' overview panel (default: current quarter).",
)
@click.option(
    "--year",
    type=click.IntRange(min=2024, max=9999),
    default=None,
    help="Year for the selected --quarter in the 'Detailed QC by LE' overview panel (default: current year).",
)
def tui(ctx: click.Context, quarter: int | None, year: int | None):
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

    app = DatabaseBrowser(database=database, public_keys=public_keys, quarter=quarter, year=year)
    app.run()


@db.command("should-qc")
@_submission_id_argument
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
    signature_status: SignatureStatus = SignatureStatus.UNKNOWN,
) -> dict[str, Any]:
    """Serialize a submission and its latest log entry to a JSON-compatible dict.

    :param log_obj: The most recent :class:`~grz_db.models.submission.SubmissionStateLog` or
        :class:`~grz_db.models.submission.ChangeRequestLog`, or ``None`` if no log exists yet.
    :param submission: The submission ORM/Pydantic model instance.
    :param signature_status: Verification result for the log entry's author signature.
        Defaults to :attr:`~SignatureStatus.UNKNOWN` when no verification was performed.
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
@_submission_id_argument
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
        _abort(e)
    except Exception as e:
        console_err.print(f"[red]An unexpected error occurred: {e}[/red]")
        raise click.ClickException(f"Failed to add submission: {e}") from e


@submission.command()
@_submission_id_argument
@click.argument("state_str", metavar="STATE", type=click.Choice(SubmissionStateEnum.list(), case_sensitive=False))
@click.option("--data", "data_json", type=str, default=None, help='Additional JSON data (e.g., \'{"k":"v"}\').')
@click.option(
    "--failure-reason",
    type=click.Choice(FailureReasonEnum.list(), case_sensitive=False),
    help="Failure reason when state is ERROR.",
)
@click.option("--ignore-error-state/--confirm-error-state")
@click.pass_context
def update(  # noqa: C901, PLR0913
    ctx: click.Context,
    submission_id: str,
    state_str: str,
    data_json: str | None,
    ignore_error_state: bool,
    failure_reason: str | None,
):
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

        failure_reason_enum = None
        if failure_reason:
            failure_reason_enum = FailureReasonEnum(failure_reason)

        new_state_log = db_service.update_submission_state(
            submission_id,
            state_enum,
            data=parsed_data,
            failure_reason=failure_reason_enum,
            grzctl_versions={k: (v if v is not None else "unknown") for k, v in get_versions().items()},
        )

        console_err.print(
            f"[green]Submission '{submission_id}' updated to state '{new_state_log.state.value}'. Log ID: {new_state_log.id}[/green]"
        )
        if new_state_log.data:
            console_err.print(f"  Data: {new_state_log.data}")
    except SubmissionNotFoundError as e:
        _abort_missing_submission(e, submission_id)
    except click.exceptions.Exit as e:
        if e.exit_code != 0:
            raise e
    except Exception as e:
        console_err.print(f"[red]An unexpected error occurred: {e}[/red]")
        traceback.print_exc()
        raise click.ClickException(f"Failed to update submission state: {e}") from e


# Exactly what SubmissionDb.modify_submission accepts. The epilog and the KEY choices both build
# from this, so --help always advertises what the command accepts. SubmissionBase, not Submission:
# the table model adds case_id, which `db case relink` sets. Same set as --allow-overwrite.
_MODIFIABLE_SUBMISSION_KEYS = sorted(SubmissionBase.model_fields.keys() - SubmissionBase.immutable_fields)


@submission.command(
    epilog="Currently available KEYs are: "
    + ", ".join(_MODIFIABLE_SUBMISSION_KEYS)
    + ". A submission's case cannot be changed here: use 'db case relink' to change it, "
    + "or 'db case unlink' to remove it."
)
@_submission_id_argument
@click.argument("key", metavar="KEY", type=click.Choice(_MODIFIABLE_SUBMISSION_KEYS))
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
        _abort_missing_submission(e, submission_id)
    except SubmissionTypeInvalidForCaseError as e:
        _abort(e)
    except Exception as e:
        console_err.print(f"[red]An unexpected error occurred: {e}[/red]")
        traceback.print_exc()
        raise click.ClickException(f"Failed to update submission state: {e}") from e


_ignore_field_option = click.option(
    "--ignore-field",
    "ignore_field",
    type=click.Choice(
        [*(SubmissionBase.model_fields.keys() - SubmissionBase.immutable_fields), "case_id"],
        case_sensitive=False,
    ),
    help="Do not populate the given field from the metadata to the database. Can be specified multiple times. "
    "Passing --ignore-field case_id skips case resolution and linking.",
    multiple=True,
)


def _prepare_submission_console_table(changes: "SubmissionChangeSet") -> rich.console.RenderableType:
    """Build a Rich renderable that shows pending submission-level metadata changes.

    :param changes: :class:`SubmissionChangeSet` instance produced by :meth:`SubmissionDb.diff`.
    :returns: A :class:`rich.table.Table` when there are pending changes, or a plain text message otherwise.
    """
    rows = [
        (d.key, str(d.diff.before) if d.diff.before is not None else _TEXT_MISSING, str(d.diff.after))
        for d in changes.fields.pending
        if d.key != "submission_metadata"
    ]
    if (link := changes.case_link) is not None:
        rows.append(
            (
                "case_id",
                str(link.before) if link.before is not None else _TEXT_MISSING,
                str(link.after) if link.after is not None else "new case",
            )
        )
    elif changes.case_link_error is not None:
        # Shown rather than raised: the operator decides whether the rest is worth recording
        # while the case itself waits on them to resolve it.
        rows.append(("case_id", _TEXT_MISSING, f"unresolved: {changes.case_link_error}"))
    if rows:
        diff_table_tbl = rich.table.Table(title="Submission Metadata")
        diff_table_tbl.add_column("Key")
        diff_table_tbl.add_column("Before")
        diff_table_tbl.add_column("After")
        for key, before, after in sorted(rows, key=lambda row: row[0]):
            diff_table_tbl.add_row(key, before, after)
        diff_table: rich.console.RenderableType = diff_table_tbl
    else:
        diff_table = rich.padding.Padding(rich.text.Text("No changes to submission-level metadata."), pad=(0, 0, 0, 0))
    return diff_table


def _report_case_link(db_service: SubmissionDb, submission_id: str) -> None:
    """Print which case a submission is linked to; call after committing a change set that includes a case link."""
    linked = db_service.get_submission(submission_id)
    if linked is not None and linked.case_id is not None:
        console_err.print(f"[green]Submission linked to case {linked.case_id}.[/green]")


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
@_submission_id_argument
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
@_ignore_field_option
@click.pass_context
def populate(  # noqa: C901, PLR0912, PLR0913
    ctx: click.Context,
    submission_id: str,
    metadata_path: str,
    submission_date: datetime | None,
    confirm: bool,
    ignore_field: tuple[str, ...],
):
    """Populate a submission in the database based on the given metadata.json file.

    Also resolves and links the submission's case; an unresolved link is shown rather
    than blocking the rest of the diff.
    """
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
        _abort_missing_submission(e, submission_id)
    except Exception as e:
        console_err.print(f"[red]An unexpected error occurred: {e}[/red]")
        traceback.print_exc()
        raise click.ClickException(f"Failed to update submission state: {e}") from e

    with open(metadata_path) as fd:
        metadata = GrzSubmissionMetadata.model_validate_json(fd.read())

    try:
        SubmissionDb.assert_metadata_not_redacted(metadata, submission_id, set(ignore_field))
    except ValueError as e:
        raise ValueError(
            f"Refusing to populate a seemingly-redacted submission: {e} "
            f"(from {metadata_path}). "
            "Add 'tan_g'/'pseudonym' to --ignore-field to bypass, "
            "or use 'grzctl db submission modify' directly."
        ) from e

    submission_uploaded_date = (
        submission_date.date() if submission_date is not None else submission.submission_uploaded_date
    )
    if submission_date is None:
        log.warning(
            "No submission date provided and submission date is missing in the database. "
            "Will use submission date from metadata.json..."
        )
        submission_uploaded_date = metadata.submission.submission_date

    try:
        changes = db_service.diff(
            submission_id,
            metadata,
            submission_uploaded_date=submission_uploaded_date,
            ignore_fields=set(ignore_field),
        )
    except SubmissionError as e:
        _abort(e)

    # build donor diff and attach Rich tables for console preview in one pass
    diff_tables: list[rich.console.RenderableType] = []
    for donor_diff in changes.donors.added + changes.donors.updated:
        diff_tables.append(
            _prepare_donor_console_table(donor_diff.changes, donor_diff.pseudonym or "", donor_diff.state)
        )
    for donor_diff in changes.donors.deleted:
        diff_tables.append(rich.text.Text(f"Donor {donor_diff.pseudonym} deleted", style="red"))

    if not changes.has_pending:
        console_err.print(
            f"[yellow]Case link unresolved: {changes.case_link_error}[/yellow]"
            if changes.case_link_error is not None
            else "[green]Database is already up to date with the provided metadata![/green]"
        )
        return

    console.print(
        rich.panel.Panel.fit(
            rich.console.Group(_prepare_submission_console_table(changes), *diff_tables, fit=True),
            title="Pending Changes",
        )
    )

    if not confirm or click.confirm(
        "Are you sure you want to commit these changes to the database?",
        default=False,
        show_default=True,
    ):
        try:
            db_service.commit_changes(submission_id, changes)
        except SubmissionError as e:
            _abort(e)
        if changes.case_link is not None:
            _report_case_link(db_service, submission_id)
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
    # Written by GRZ_QC_Workflow >= v2.1.0. Optional so reports from older workflow
    # versions (without the column) still parse. Auto-aliased to "grzQcWorkflowVersion"
    # via StrictBaseModel's to_camel alias generator.
    grz_qc_workflow_version: str | None = None


@submission.command()
@_submission_id_argument
@click.argument("report_csv_path", metavar="path/to/report.csv", type=grzcli.FILE_R_E)
@click.option(
    "--qc-workflow-version",
    type=str,
    required=False,
    default=None,
    envvar="GRZCTL_QC_WORKFLOW_VERSION",
    help="QC workflow version to record when the report has no 'grzQcWorkflowVersion' column. "
    "If the report carries one, it takes precedence. Env: GRZCTL_QC_WORKFLOW_VERSION.",
)
@click.option(
    "--confirm/--no-confirm",
    default=True,
    help="Whether to confirm changes before committing to database. (Default: confirm)",
)
@click.pass_context
def populate_qc(
    ctx: click.Context, submission_id: str, report_csv_path: str, qc_workflow_version: str | None, confirm: bool
):
    """Populate the submission database from a detailed QC pipeline report."""
    db = ctx.obj["db_url"]
    db_service = get_submission_db_instance(db, author=ctx.obj["author"])

    with open(report_csv_path, encoding="utf-8", newline="") as report_csv_file:
        reader = csv.reader(report_csv_file)
        header = next(reader)
        reports = []
        for row in reader:
            reports.append(QCReportRow(**dict(zip(header, row, strict=True))))

    # The report (GRZ_QC_Workflow >= v2.1.0) carries its own version in the
    # 'grzQcWorkflowVersion' column; treat that as the source of truth and only fall back to
    # --qc-workflow-version for older reports that lack the column.
    report_versions = {report.grz_qc_workflow_version for report in reports if report.grz_qc_workflow_version}
    if len(report_versions) > 1:
        raise click.ClickException(f"Inconsistent grzQcWorkflowVersion values in report: {sorted(report_versions)}")
    report_version = report_versions.pop() if report_versions else None

    if report_version and qc_workflow_version and report_version != qc_workflow_version:
        raise click.ClickException(
            f"--qc-workflow-version ({qc_workflow_version}) disagrees with the report's "
            f"grzQcWorkflowVersion ({report_version}). Omit --qc-workflow-version to use the report "
            "value, which is authoritative."
        )

    effective_qc_workflow_version = report_version or qc_workflow_version
    if not effective_qc_workflow_version:
        raise click.ClickException(
            "No QC workflow version found: the report has no 'grzQcWorkflowVersion' column and "
            "--qc-workflow-version was not provided."
        )

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
                qc_workflow_version=effective_qc_workflow_version,
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
@_submission_id_argument
@click.argument("change_str", metavar="CHANGE", type=click.Choice(ChangeRequestEnum.list(), case_sensitive=False))
@click.option("--data", "data_json", type=str, default=None, help='Inline JSON data (e.g., \'{"k":"v"}\').')
@click.option(
    "--data-file",
    "data_file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Path to a JSON or YAML file with the change-request fields (see `grzctl change-request-template`).",
)
@click.option(
    "--raw-content",
    "raw_content_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Optional path to a binary file (e.g. a .pdf or .png) accompanying the request. "
        "Type is inferred from the file extension and verified by magic bytes."
    ),
)
@click.option(
    "--dry-run",
    "dry_run",
    is_flag=True,
    default=False,
    help="Validate inputs and check the submission exists, but do not write the change request.",
)
@click.pass_context
def change_request(  # noqa: PLR0913
    ctx: click.Context,
    submission_id: str,
    change_str: str,
    data_json: str | None,
    data_file: Path | None,
    raw_content_path: Path | None,
    dry_run: bool,
):
    """Register a completed change request for the given submission.

    The audit fields (requester name, email, requested-at, request content) are required.
    See ``grzctl change-request-template`` for a fill-in YAML template, and
    ``packages/grzctl/examples/demo_change_request.py`` for a runnable end-to-end
    walkthrough including the optional ``--raw-content`` (PDF/PNG) attachment path.
    """
    try:
        change_request_enum = ChangeRequestEnum(change_str)
    except ValueError as e:
        console_err.print(f"[red]Error: Invalid change request value '{change_str}'.[/red]")
        raise click.Abort() from e

    kwargs = resolve_and_validate_change_request(change_request_enum, data_json, data_file, raw_content_path)

    db = ctx.obj["db_url"]
    db_service = get_submission_db_instance(db, author=ctx.obj["author"])

    if dry_run:
        existing = db_service.get_submission(submission_id)
        if existing is None:
            console_err.print(
                f"[red]Dry run: submission '{submission_id}' not found. "
                f"You might need to add it first: grz-cli db submission add {submission_id}[/red]"
            )
            raise click.Abort()
        console_err.print(
            f"[yellow]Dry run: would register change request '{change_request_enum.value}' "
            f"for submission '{submission_id}'. No changes were written.[/yellow]"
        )
        console_err.print("[yellow]Validated fields:[/yellow]")
        preview = {k: v for k, v in kwargs.items() if k != "request_raw_content"}
        if kwargs["request_raw_content"] is not None:
            preview["request_raw_content"] = f"<{len(kwargs['request_raw_content'])} bytes>"
        click.echo(json.dumps(preview, indent=2, ensure_ascii=False, default=str), err=True)
        return

    try:
        new_change_request_log = db_service.add_change_request(submission_id, change_request_enum, **kwargs)
        console_err.print(
            f"[green]Submission '{submission_id}' has undergone a change request of '{new_change_request_log.change.value}'. Log ID: {new_change_request_log.id}[/green]"
        )
        if new_change_request_log.data:
            console_err.print(f"  Data: {new_change_request_log.data}")

    except SubmissionNotFoundError as e:
        _abort_missing_submission(e, submission_id)
    except Exception as e:
        console_err.print(f"[red]An unexpected error occurred: {e}[/red]")
        traceback.print_exc()
        raise click.ClickException(f"Failed to update submission state: {e}") from e


def _research_consented_now(submission: Submission) -> bool | None:
    """Research consent for the submission re-evaluated as of now.

    Unlike the persisted ``consented`` field (evaluated at the submission date),
    this recomputes consent from the stored redacted metadata using the current date.

    :param submission: Submission whose stored metadata to evaluate.
    :returns: ``True``/``False`` for the consent decision now, or ``None`` when
        no metadata is stored (e.g. rows migrated without backpopulated metadata)
        or when the stored metadata cannot be parsed.
    """
    if not submission.submission_metadata:
        return None
    try:
        metadata = GrzSubmissionMetadata.model_validate(submission.submission_metadata)
    except ValidationError:
        log.debug("Could not parse stored metadata for submission %s to evaluate consent now.", submission.id)
        return None
    return metadata.consents_to_research(date=date.today())


def _build_attribute_table(submission: Submission, research_consented_now: bool | None) -> rich.table.Table:
    """Build the attribute table shown by ``submission show``.

    :param submission: Submission to render.
    :param research_consented_now: Research consent re-evaluated as of now
        (see :func:`_research_consented_now`), or ``None`` when unavailable.
    :returns: A populated rich table of submission attributes.
    """
    rows: list[tuple[str, Any]] = []
    for label, attr_name in (
        ("tanG", "tan_g"),
        ("Pseudonym", "pseudonym"),
        ("Submission Uploaded Date", "submission_uploaded_date"),
        ("Submission Size", "submission_size"),
        ("Submission Type", "submission_type"),
        ("Submitter ID", "submitter_id"),
        ("Case ID", "case_id"),
        ("Data Node ID", "data_node_id"),
        ("Disease Type", "disease_type"),
        ("Genomic Study Type", "genomic_study_type"),
        ("Genomic Study Subtype", "genomic_study_subtype"),
        ("Basic QC Passed", "basic_qc_passed"),
        ("Research consent (at submission)", "consented"),
        ("Selected For QC", "selected_for_qc"),
        ("Detailed QC Passed", "detailed_qc_passed"),
    ):
        rows.append((label, getattr(submission, attr_name)))
        if attr_name == "consented":
            # Adjacent row: research consent re-evaluated as of now (recomputed from stored metadata).
            rows.append(("Research consent (now)", research_consented_now))
    return _attribute_table(rows)


@submission.command("show")
@_submission_id_argument
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

    research_consented_now = _research_consented_now(submission)

    if output_json:
        submission_dict = submission.model_dump(mode="json")
        submission_dict["research_consented_now"] = research_consented_now
        submission_dict["states"] = []

        for state_log in sorted(submission.states, key=lambda s: s.timestamp):
            signature_status, verifying_key_comment = _verify_signature(
                ctx.obj["public_keys"], state_log.author_name, state_log
            )
            state_dict = state_log.model_dump(
                mode="json", include={"id", "timestamp", "state", "data", "failure_reason", "grzctl_versions"}
            )

            state_dict["data_steward"] = state_log.author_name
            state_dict["data_steward_signature"] = signature_status
            state_dict["signature_key_comment"] = verifying_key_comment
            submission_dict["states"].append(state_dict)

        json.dump(submission_dict, sys.stdout)
        sys.stdout.write("\n")
        return

    attribute_table = _build_attribute_table(submission, research_consented_now)

    renderables: list[rich.console.RenderableType] = [rich.padding.Padding(attribute_table, (1, 0))]
    if submission.states:
        state_table = rich.table.Table(title="State History", show_header=True)
        state_table.add_column("Log ID", style="dim", width=12)
        state_table.add_column("Timestamp (UTC)", style="yellow")
        state_table.add_column("State", style="green")
        state_table.add_column("Failure Reason", style="red", min_width=15)
        state_table.add_column("Data", style="cyan", overflow="ellipsis")
        state_table.add_column("Dependency Versions", style="blue")
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
                state_log.failure_reason.value if state_log.failure_reason else "",
                data_str,
                json.dumps(state_log.grzctl_versions) if state_log.grzctl_versions else _TEXT_MISSING,
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


def _fetch_metadata_json(s3_client: Any, bucket: str, submission_id: str) -> str | None:
    """Return the raw metadata.json content for *submission_id*, or None when not found.

    Raises for any S3 error that is not a simple 404/NoSuchKey.
    """
    key = f"{submission_id}/metadata/metadata.json"
    try:
        response = s3_client.get_object(Bucket=bucket, Key=key)
        return response["Body"].read().decode("utf-8")
    except botocore.exceptions.ClientError as exc:
        code = exc.response.get("Error", {}).get("Code", "")
        if code in {"404", "NoSuchKey"}:
            return None
        raise


class _BackfillResult(StrEnum):
    UPDATED = "updated"
    UP_TO_DATE = "up_to_date"
    NOT_FOUND = "not_found"
    WOULD_OVERWRITE = "would_overwrite"
    ERROR = "error"
    CONSENT_MISMATCH = "consent_mismatch"


class _BackfillOutcome(NamedTuple):
    """What one submission's backfill did, and whether its case link is still missing.

    The two are independent: a submission whose case key could not be resolved is still
    written, since the reason is rarely about that submission and the rest of the metadata is
    what the Prüfbericht is built from. Reporting the unresolved link *instead* of the outcome
    would count such a submission as not updated when it was.

    :param status: What happened to the submission's fields, donors and case link.
    :param link_unresolved: Whether the case could not be resolved. Re-running once the cause
        is cleared links it.
    """

    status: _BackfillResult
    link_unresolved: bool = False


def _backfill_submission(  # noqa: C901, PLR0911, PLR0913
    current_submission: Submission,
    s3_client: Any,
    bucket: str,
    db_service: SubmissionDb,
    dry_run: bool,
    force: bool,
    ignore_fields: set[str],
    allow_overwrite: frozenset[str] = frozenset(),
) -> _BackfillOutcome:
    """Fetch metadata.json from S3 for one submission and commit a diff to the database.

    Uses the same :func:`SubmissionDb.diff` / :func:`SubmissionDb.commit_changes` path
    as ``grzctl db submission populate`` so that every derived field (not only
    *submission_size* and *submission_metadata*) is kept consistent, donor records are
    synchronised, the case link is resolved and committed (``"case_id"`` in
    *ignore_fields* disables that), and already-up-to-date submissions are detected
    without a write.

    When *force* is False, a destructive diff (existing non-NULL field would change)
    is skipped instead of committed, preserving manually-corrected values. The caller
    is expected to pre-filter already-populated rows so re-runs do not re-pay the S3
    network cost for them.
    """
    submission_id = current_submission.id

    try:
        raw_json = _fetch_metadata_json(s3_client, bucket, submission_id)
    except Exception as exc:
        console_err.print(f"[red]  {submission_id}: S3 error: {exc}[/red]")
        return _BackfillOutcome(_BackfillResult.ERROR)

    if raw_json is None:
        # this is expected for submissions residing in the other consent bucket, so we do not explicitly log that here
        # but still report them in the final stats
        return _BackfillOutcome(_BackfillResult.NOT_FOUND)

    try:
        metadata = GrzSubmissionMetadata.model_validate_json(raw_json)
    except Exception as exc:
        console_err.print(f"[red]  {submission_id}: failed to parse metadata.json: {exc}[/red]")
        return _BackfillOutcome(_BackfillResult.ERROR)

    # The archived copy is redacted; the stored row holds the submitter's values.
    still_redacted = current_submission.restore_redacted_fields(metadata)
    ignore_fields = set(ignore_fields) | still_redacted

    try:
        # Backfill has no authoritative date to offer: the archive does not carry one, and the
        # submitter's declared submission_date is not when the upload finished, which is what this
        # column means and what the reporting windows filter on. populate reads the real one from
        # the S3 object; backfill does not fetch it. Passing None lets diff() supply a value for
        # construction and then exclude the field, so backfill never writes it. Passing the row's
        # own value instead would compare a snapshot taken before the run against the row diff()
        # re-reads, and write the stale one.
        changes = db_service.diff(
            submission_id,
            metadata,
            submission_uploaded_date=None,
            ignore_fields=ignore_fields or None,
        )
    except Exception as exc:
        console_err.print(f"[red]  {submission_id}: diff failed: {exc}[/red]")
        return _BackfillOutcome(_BackfillResult.ERROR)

    # See _BackfillOutcome: only the link is withheld, not the rest of the submission.
    link_unresolved = changes.case_link_error is not None
    if link_unresolved:
        console_err.print(f"[yellow]  {submission_id}: case link unresolved: {changes.case_link_error}[/yellow]")

    if not changes.has_pending:
        console_err.print(f"[dim]  {submission_id}: already up to date, skipping.[/dim]")
        return _BackfillOutcome(_BackfillResult.UP_TO_DATE, link_unresolved)

    # Filling a field that was NULL destroys nothing, so it is always written. Replacing one that
    # already has a value needs saying so: --force permits every such overwrite, --allow-overwrite
    # only the fields it names, and anything else is held back and reported.
    allowed = {diff.key for diff in changes.fields.pending} if force else allow_overwrite
    changes.fields, withheld = changes.fields.withhold_destructive(allowed)
    if withheld:
        console_err.print(
            f"[dim]  {submission_id}: not overwriting {', '.join(diff.key for diff in withheld)} "
            f"(use --force for all, or --allow-overwrite for named fields).[/dim]"
        )

    # --allow-overwrite is built from SubmissionBase, which has no case_id, so it cannot name the
    # case link. Replacing one would undo a deliberate case relink, so an unpermitted change holds
    # the whole submission back.
    if not force and changes.case_link is not None and changes.case_link.state is DiffState.UPDATED:
        console_err.print(
            f"[dim]  {submission_id}: would overwrite the case link (case {changes.case_link.before}), "
            "skipping (use --force to overwrite).[/dim]"
        )
        return _BackfillOutcome(_BackfillResult.WOULD_OVERWRITE)

    if not changes.has_pending:
        return _BackfillOutcome(_BackfillResult.WOULD_OVERWRITE)

    if dry_run:
        console_err.print(
            f"[yellow]  [dry-run] {submission_id}: would update: {', '.join(changes.pending_changes)}[/yellow]"
        )
        return _BackfillOutcome(_BackfillResult.UPDATED, link_unresolved)

    try:
        db_service.commit_changes(submission_id, changes)
        console_err.print(f"[green]  {submission_id}: updated ({', '.join(changes.pending_changes)}).[/green]")
        return _BackfillOutcome(_BackfillResult.UPDATED, link_unresolved)
    except Exception as exc:
        console_err.print(f"[red]  {submission_id}: failed to commit: {exc}[/red]")
        return _BackfillOutcome(_BackfillResult.ERROR)


@db.command("backfill")
@grzctl_configuration
@click.option(
    "--dry-run/--no-dry-run",
    default=False,
    help="Preview which submissions would be updated without writing to the database.",
)
@click.option(
    "--force/--no-force",
    default=False,
    help="Overwrite existing non-NULL fields when the metadata.json value differs (destructive diffs). "
    "Without this flag, such fields are reported and left alone while the rest is still written.",
)
@click.option(
    "--allow-overwrite",
    "allow_overwrite",
    type=click.Choice(list(SubmissionBase.model_fields.keys() - SubmissionBase.immutable_fields), case_sensitive=False),
    multiple=True,
    help="Overwrite only these existing non-NULL fields when the metadata.json value differs "
    "(may be repeated). Other destructive changes are held back and the submission is updated in "
    "part. Mutually exclusive with --force, which permits every overwrite.",
)
@click.option(
    "--submission-id",
    "submission_ids",
    multiple=True,
    metavar="SUBMISSION_ID",
    help="Restrict backfill to these submission IDs (may be repeated). Mutually exclusive with --start-date/--end-date.",
)
@click.option(
    "--start-date",
    type=click.DateTime(formats=["%Y-%m-%d"]),
    default=datetime.min,
    help="Process only submissions processed on or after this date (inclusive). Defaults to the beginning of time.",
)
@click.option(
    "--end-date",
    type=click.DateTime(formats=["%Y-%m-%d"]),
    default=datetime.max,
    help="Process only submissions processed on or before this date (inclusive). Defaults to the end of time.",
)
@_ignore_field_option
@click.pass_context
def backfill(  # noqa: C901, PLR0912, PLR0913, PLR0915
    ctx: click.Context,
    configuration: GrzctlConfig,
    dry_run: bool,
    force: bool,
    allow_overwrite: tuple[str, ...],
    submission_ids: tuple[str, ...],
    start_date: datetime,
    end_date: datetime,
    ignore_field: tuple[str, ...],
    **kwargs,
):
    r"""Backfill submission fields for existing submissions by re-reading metadata.json from S3.

    Uses the same diff/commit path as ``grzctl db submission populate``: only fields
    that are actually missing or changed are written, donor records are synchronised,
    and already-up-to-date submissions are silently skipped.

    A field that is missing in the database is always filled. One that already has a
    different value is only overwritten with --force, or when --allow-overwrite names it;
    any other overwrite is reported and held back, so a submission is updated in part
    rather than skipped entirely.

    Both the consented and non-consented archive buckets are always scanned.
    If a submission's metadata.json is found in both, an error is raised.

    The case link is resolved and written like any other missing value, but --allow-overwrite
    cannot name it, so replacing an existing link holds the whole submission back until
    --force is passed. Pass --ignore-field case_id to skip case linking altogether.

    Candidate selection (mutually exclusive):

    \b
      (default)                all submissions within the date window (defaults to the full historical range)
      --submission-id ...      explicit list of submission IDs
      --start-date/--end-date  narrow the default date window

    This command is idempotent: re-running it is always safe.
    """
    # ── Validate option combinations ────────────────────────────────────────
    if submission_ids and (start_date != datetime.min or end_date != datetime.max):
        raise click.UsageError("--submission-id and --start-date/--end-date are mutually exclusive.")

    if force and allow_overwrite:
        raise click.UsageError("--force and --allow-overwrite are mutually exclusive; --force already permits all.")

    ignore_fields = set(ignore_field)

    # ── Build S3 clients for both archive targets ────────────────────────────
    archive_targets = [
        ("consented", configuration.archives.consented.s3.bucket, init_s3_client(configuration.archives.consented.s3)),
        (
            "non_consented",
            configuration.archives.non_consented.s3.bucket,
            init_s3_client(configuration.archives.non_consented.s3),
        ),
    ]

    db_service = get_submission_db_instance(ctx.obj["db_url"], author=ctx.obj["author"])

    # ── Determine which submissions to process ──────────────────────────────
    if submission_ids:
        candidates: list[Submission] = []
        for sid, sub in zip(submission_ids, db_service.get_submissions(list(submission_ids)), strict=True):
            if sub is None:
                console_err.print(f"[yellow]Warning: submission '{sid}' not found in database, skipping.[/yellow]")
            else:
                candidates.append(sub)
    else:
        candidates = list(db_service.list_processed_between(start_date.date(), end_date.date()))
        console_err.print(
            f"[cyan]Date window: {start_date.date()} to {end_date.date()} ({len(candidates)} submission(s)).[/cyan]"
        )

    counts: Counter[_BackfillResult] = Counter()

    console_err.print(
        f"[cyan]{'[dry-run] ' if dry_run else ''}Processing {len(candidates)} submission(s) "
        f"across consented and non-consented archives…[/cyan]"
    )

    # ── Fetch metadata from S3 and update DB ────────────────────────────────
    consent_mismatches = 0
    expired_consents = 0
    links_unresolved = 0

    for submission in tqdm(candidates):
        # fetch from archive targets
        results: dict[str, _BackfillOutcome] = {}
        for label, bucket, client in archive_targets:
            results[label] = _backfill_submission(
                submission,
                client,
                bucket,
                db_service,
                dry_run,
                force,
                ignore_fields,
                frozenset(allow_overwrite),
            )

        found_in = [label for label, res in results.items() if res.status != _BackfillResult.NOT_FOUND]

        if len(found_in) > 1:
            console_err.print(
                f"[red]  {submission.id}: ERROR: metadata.json found in both consented and non_consented archives[/red]"
            )
            counts[_BackfillResult.ERROR] += 1
            continue

        if not found_in:
            counts[_BackfillResult.NOT_FOUND] += 1
            continue

        actual_archive = found_in[0]  # "consented" or "non_consented"
        counts[results[actual_archive].status] += 1
        links_unresolved += results[actual_archive].link_unresolved

        # check DB `consented` vs actual archive
        if submission.consented is not None:
            expected_archive = "consented" if submission.consented else "non_consented"
            if actual_archive != expected_archive:
                consent_mismatches += 1
                console_err.print(
                    f"[bold red]  CONSENT MISMATCH: {submission.id} has DB consented={submission.consented} "
                    f"(expected '{expected_archive}'), but was found in '{actual_archive}' bucket![/bold red]"
                )

        # check consent expiration (if stored in consented archive)
        # (submission.submission_metadata is populated after commit, or check parsed metadata)
        if actual_archive == "consented" and submission.submission_metadata:
            try:
                meta = GrzSubmissionMetadata.model_validate(submission.submission_metadata)
                if not meta.consents_to_research(date=date.today()):
                    expired_consents += 1
                    console_err.print(
                        f"[yellow]  CONSENT EXPIRED: {submission.id} is in 'consented' archive, "
                        f"but research consent has expired as of today ({date.today()}).[/yellow]"
                    )
            except ValidationError:
                log.warning(f"Error validating submission metadata for {submission.id}: {traceback.format_exc()}")

    # ── Summary ─────────────────────────────────────────────────────────────
    prefix = "[dry-run] " if dry_run else ""
    verb = "Would update" if dry_run else "Updated"
    console_err.print(
        f"\n[cyan]{prefix}Done. {verb}: {counts[_BackfillResult.UPDATED]}\n"
        f"  Up to date: {counts[_BackfillResult.UP_TO_DATE]}\n"
        f"  Not in bucket (split consent): {counts[_BackfillResult.NOT_FOUND]}\n"
        f"  Would overwrite (needs --force): {counts[_BackfillResult.WOULD_OVERWRITE]}\n"
        f"  Case link unresolved: {links_unresolved} (also counted above)\n"
        f"  Errors: {counts[_BackfillResult.ERROR]}[/cyan]"
    )

    if consent_mismatches or expired_consents:
        console_err.print(
            f"\n[bold yellow]Consent issues:\n"
            f"  Bucket ↔ DB consent mismatches: {consent_mismatches}\n"
            f"  Expired consents in consented archive: {expired_consents}[/bold yellow]"
        )
    if counts[_BackfillResult.ERROR]:
        sys.exit(1)


@db.command("sync-from-inbox")
@grzctl_configuration
@inbox_options
@click.pass_context
def sync_from_inbox(
    ctx: click.Context,
    configuration: GrzctlConfig,
    submitter_id: str,
    inbox_name: str,
    **kwargs,
):
    """Synchronize the database with submissions found in the inbox."""
    try:
        s3_options = configuration.resolve_inbox(submitter_id=submitter_id, inbox_name=inbox_name).s3
    except Exception:
        console_err.print(
            f"[red]Error resolving S3 configuration for inbox '{inbox_name}': {traceback.format_exc()}[/red]"
        )
        sys.exit(1)

    db_url = ctx.obj["db_url"]
    author = ctx.obj["author"]
    db_service = get_submission_db_instance(db_url, author=author)

    bucket_name = s3_options.bucket
    inbox_desc = f"'{inbox_name}' (bucket '{bucket_name}')" if inbox_name != bucket_name else f"'{bucket_name}'"

    try:
        console_err.print(f"[cyan]Scanning inbox {inbox_desc}...[/cyan]")
        s3_submissions = query_submissions(s3_options, show_cleaned=False)

        console_err.print(f"[cyan]Synchronizing {len(s3_submissions)} submissions with database...[/cyan]")
        sync_submissions(db_service, s3_submissions, author)

        console_err.print("[green]Synchronization complete.[/green]")

    except Exception:
        console_err.print(f"[red]Error during synchronization: {traceback.format_exc()}[/red]")
        traceback.print_exc()
        sys.exit(1)
