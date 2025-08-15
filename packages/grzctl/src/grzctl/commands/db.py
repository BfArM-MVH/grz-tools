"""Command for managing a submission database"""

import enum
import json
import logging
import sys
import traceback
from collections import namedtuple
from datetime import date
from typing import Any

import click
import rich.console
import rich.padding
import rich.panel
import rich.table
import rich.text
from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PublicKey
from grz_common.cli import config_file, output_json
from grz_common.constants import REDACTED_TAN
from grz_db.errors import (
    DatabaseConfigurationError,
    DuplicateSubmissionError,
    DuplicateTanGError,
    SubmissionNotFoundError,
)
from grz_db.models.author import Author
from grz_db.models.base import VerifiableLog
from grz_db.models.submission import (
    ChangeRequestEnum,
    ChangeRequestLog,
    Submission,
    SubmissionDb,
    SubmissionStateEnum,
    SubmissionStateLog,
)
from grz_pydantic_models.submission.metadata import GrzSubmissionMetadata, LibraryType

from ..models.config import DbConfig
from . import limit
from .pruefbericht import get_pruefbericht_library_type

console = rich.console.Console()
console_err = rich.console.Console(stderr=True)
log = logging.getLogger(__name__)
_TEXT_MISSING = rich.text.Text("missing", style="italic yellow")


def get_submission_db_instance(db_url: str, author: Author | None = None) -> SubmissionDb:
    """Creates and returns an instance of SubmissionDb."""
    return SubmissionDb(db_url=db_url, author=author)


@click.group(help="Database operations")
@config_file
@click.pass_context
def db(ctx: click.Context, config_file: str):
    """Database operations"""
    config = DbConfig.from_path(config_file)
    db_config = config.db
    if not db_config:
        raise ValueError("DB config not found")
    author_name = db_config.author.name

    if path := db_config.author.private_key_path:
        with open(path, "rb") as f:
            private_key_bytes = f.read()
    elif key := db_config.author.private_key:
        private_key_bytes = key.encode("utf-8")
    else:
        raise ValueError("Either private_key or private_key_path must be provided.")

    from cryptography.hazmat.primitives.serialization import load_ssh_public_key

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
    ctx.obj = {"author": author, "public_keys": public_keys, "db_url": db_config.database_url}


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
@output_json
@limit
@click.pass_context
def list_submissions(ctx: click.Context, output_json: bool, limit: int):
    """Lists all submissions in the database with their latest state."""
    db = ctx.obj["db_url"]
    db_service = get_submission_db_instance(db)

    try:
        submissions = db_service.list_submissions(limit=limit)
    except Exception as e:
        raise click.ClickException(str(e)) from e

    if not submissions:
        console_err.print("[yellow]No submissions found in the database.[/yellow]")
        return

    table = rich.table.Table(title="All Submissions")
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
@output_json
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


class SignatureStatus(enum.StrEnum):
    """Enum for signature status."""

    VERIFIED = "Verified"
    FAILED = "Failed"
    ERROR = "Error"
    UNKNOWN = "Unknown"

    def rich_display(self, comment: str | None) -> str:
        """Displays the signature status in rich format."""
        match self:
            case "Verified":
                return "[green]Verified[/green]" if comment is None else f"[green]Verified ({comment})[/green]"
            case "Failed":
                return "[red]Failed[/red]"
            case "Error":
                return "[red]Error[/red]"
            case "Unknown" | _:
                return "[yellow]Unknown Key[/yellow]"


def _verify_signature(
    public_keys: dict[str, Ed25519PublicKey], expected_key_comment: str, verifiable_log: VerifiableLog
) -> tuple[SignatureStatus, str | None]:
    signature_status = SignatureStatus.UNKNOWN
    verifying_key_comment = None
    if public_key := public_keys.get(expected_key_comment):
        try:
            signature_status = SignatureStatus.VERIFIED if verifiable_log.verify(public_key) else SignatureStatus.FAILED
        except Exception as e:
            signature_status = SignatureStatus.ERROR
            log.error(e)
    else:
        log.info("Found no key with matching username in comment, trying all keys")
        for comment, public_key in public_keys.items():
            try:
                if verifiable_log.verify(public_key):
                    signature_status = SignatureStatus.VERIFIED
                    verifying_key_comment = comment
                    # stop trying after first verification success
                    break
            except Exception as e:
                signature_status = SignatureStatus.ERROR
                log.error(e)

    return signature_status, verifying_key_comment


def _build_submission_dict_from(
    log_obj: SubmissionStateLog | ChangeRequestLog | None,
    submission: Submission,
    signature_status: SignatureStatus,
) -> dict[str, Any]:
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
    except (DuplicateSubmissionError, DuplicateTanGError) as e:
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


def _diff_metadata(
    submission: Submission, metadata: GrzSubmissionMetadata, ignore_fields: set[str]
) -> list[tuple[str, Any, Any]]:
    """Given a database submission and a metadata.json file, report changed fields and their before/after values if they are not in ignore_fields."""
    changes = []

    simple_fields = {"tan_g", "submission_date", "submission_type", "submitter_id", "disease_type"}

    for field in simple_fields - ignore_fields:
        if field == "tan_g" and metadata.submission.tan_g == REDACTED_TAN:
            raise ValueError(
                "Refusing to populate a seemingly-redacted TAN (all zeros). "
                "Add 'tan_g' to --ignore-field or use 'grzctl db submission modify' directly."
            )
        submission_attr = getattr(submission, field)
        metadata_attr = getattr(metadata.submission, field)
        if submission_attr != metadata_attr:
            changes.append((field, submission_attr, metadata_attr))

    # pseudonym (TODO: change after phase 0)
    if "pseudonym" not in ignore_fields and (submission.pseudonym != metadata.submission.local_case_id):
        if not metadata.submission.local_case_id:
            raise ValueError(
                "Refusing to populate a seemingly-redacted local case ID (empty). "
                "Add 'pseudonym' to --ignore-field or use 'grzctl db submission modify' directly."
            )
        changes.append(("pseudonym", submission.pseudonym, metadata.submission.local_case_id))

    # data node id
    if "data_node_id" not in ignore_fields and (submission.data_node_id != metadata.submission.genomic_data_center_id):
        changes.append(("data_node_id", submission.data_node_id, metadata.submission.genomic_data_center_id))

    # library type
    metadata_library_type = LibraryType(get_pruefbericht_library_type(metadata))
    if "library_type" not in ignore_fields and (submission.library_type != metadata_library_type):
        changes.append(("library_type", submission.library_type, metadata_library_type))

    # consent state
    consented = metadata.consents_to_research(date=date.today())
    if submission.consented != consented:
        changes.append(("consented", submission.consented, consented))

    return changes


@submission.command()
@click.argument("submission_id", type=str)
@click.argument("metadata_path", metavar="path/to/metadata.json", type=str)
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
def populate(ctx: click.Context, submission_id: str, metadata_path: str, confirm: bool, ignore_field: list[str]):
    """Populate the submission database from a metadata JSON file."""
    log.debug(f"Ignored fields for populate: {ignore_field}")

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

    with open(metadata_path) as metadata_file:
        metadata = GrzSubmissionMetadata.model_validate_json(metadata_file.read())

    changes = _diff_metadata(submission, metadata, set(ignore_field))
    if not changes:
        console_err.print("[green]Database is already up to date with the provided metadata![/green]")
        ctx.exit()

    diff_table = rich.table.Table()
    diff_table.add_column("Key")
    diff_table.add_column("Before")
    diff_table.add_column("After")
    for key, before, after in changes:
        diff_table.add_row(key, str(before) if before is not None else _TEXT_MISSING, str(after))
    panel = rich.panel.Panel.fit(diff_table, title="Pending Changes")
    console.print(panel)

    if not confirm or click.confirm(
        "Are you sure you want to commit these changes to the database?",
        default=False,
        show_default=True,
    ):
        for key, _before, after in changes:
            _ = db_service.modify_submission(submission_id, key=key, value=after)
        console_err.print("[green]Database populated successfully.[/green]")


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
@click.pass_context
def show(ctx: click.Context, submission_id: str):
    """
    Show details of a submission.
    """
    db = ctx.obj["db_url"]
    db_service = get_submission_db_instance(db)
    submission = db_service.get_submission(submission_id)
    if not submission:
        console_err.print(f"[red]Error: Submission with ID '{submission_id}' not found.[/red]")
        raise click.Abort()

    attribute_table = rich.table.Table(box=None)
    attribute_table.add_column("Attribute", justify="right")
    attribute_table.add_column("Value")
    for label, attr_name in (
        ("tanG", "tan_g"),
        ("Pseudonym", "pseudonym"),
        ("Submission Date", "submission_date"),
        ("Submission Type", "submission_type"),
        ("Submitter ID", "submitter_id"),
        ("Data Node ID", "data_node_id"),
        ("Disease Type", "disease_type"),
        ("Library Type", "library_type"),
        ("Basic QC Passed", "basic_qc_passed"),
        ("Consented", "consented"),
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
