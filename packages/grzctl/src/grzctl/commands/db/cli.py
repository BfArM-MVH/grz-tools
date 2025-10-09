"""Command for managing a submission database"""

import csv
import json
import logging
import sys
import traceback
from collections import namedtuple
from datetime import UTC, date, datetime
from enum import StrEnum
from operator import itemgetter
from pathlib import Path
from typing import Any

import click
import rich.console
import rich.padding
import rich.panel
import rich.table
import rich.text
import textual.logging
from grz_common.cli import FILE_R_E, config_file, output_json
from grz_common.logging import LOGGING_DATEFMT, LOGGING_FORMAT
from grz_db.errors import (
    DatabaseConfigurationError,
    DuplicateSubmissionError,
    DuplicateTanGError,
    SubmissionNotFoundError,
)
from grz_db.models.author import Author
from grz_db.models.submission import (
    ChangeRequestEnum,
    ChangeRequestLog,
    DetailedQCResult,
    Donor,
    Submission,
    SubmissionDb,
    SubmissionStateEnum,
    SubmissionStateLog,
)
from grz_pydantic_models.common import StrictBaseModel
from grz_pydantic_models.submission.metadata import (
    REDACTED_TAN,
    GenomicStudySubtype,
    GrzSubmissionMetadata,
    LibraryType,
    Relation,
    SequenceSubtype,
    SequenceType,
)
from pydantic import Field

from ...models.config import DbConfig
from .. import limit
from . import SignatureStatus, _verify_signature
from .tui import DatabaseBrowser

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

    simple_fields = {
        "tan_g",
        "submission_date",
        "submission_type",
        "submitter_id",
        "coverage_type",
        "disease_type",
        "genomic_study_type",
        "genomic_study_subtype",
    }

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

    # consent state
    consented = metadata.consents_to_research(date=date.today())
    if submission.consented != consented:
        changes.append(("consented", submission.consented, consented))

    return changes


def _diff_donors(
    donors_in_db: tuple[Donor, ...], metadata: GrzSubmissionMetadata
) -> tuple[tuple[Donor, ...], tuple[Donor, ...], tuple[rich.console.RenderableType, ...]]:
    pseudonym2before = {donor.pseudonym: donor for donor in donors_in_db}

    updated_donors = []
    pending_pseudonyms = set()
    donor_diff_tables: list[rich.console.RenderableType] = []
    for donor in metadata.donors:
        donor_after = Donor(
            submission_id=metadata.submission_id,
            pseudonym="index" if donor.relation == Relation.index_ else donor.donor_pseudonym,
            relation=donor.relation,
            library_types={datum.library_type for datum in donor.lab_data},
            sequence_types={datum.sequence_type for datum in donor.lab_data},
            sequence_subtypes={datum.sequence_subtype for datum in donor.lab_data},
            mv_consented=donor.consents_to_mv(),
            research_consented=donor.consents_to_research(date=date.today()),
            research_consent_missing_justifications={
                consent.no_scope_justification
                for consent in donor.research_consents
                if consent.no_scope_justification is not None
            }
            if donor.research_consents
            else None,
        )
        donor_before = pseudonym2before.get(donor_after.pseudonym, None)
        pending_pseudonyms.add(donor_after.pseudonym)
        if donor_before == donor_after:
            continue
        updated_donors.append(donor_after)

        table_title = f"[green]Donor '{donor_after.pseudonym}' added/updated[/green]"
        diff_table = rich.table.Table(title=table_title, min_width=len(table_title), title_justify="left")
        diff_table.add_column("Key")
        diff_table.add_column("Before")
        diff_table.add_column("After")
        for field in sorted(donor_after.model_fields.keys() - {"submission_id", "pseudonym"}):
            before = getattr(donor_before, field, None)
            after = getattr(donor_after, field)
            if before != after:
                diff_table.add_row(
                    field, _TEXT_MISSING if before is None else rich.pretty.Pretty(before), rich.pretty.Pretty(after)
                )
        if diff_table.row_count:
            donor_diff_tables.append(diff_table)

    deleted_donors = tuple(filter(lambda donor: donor.pseudonym not in pending_pseudonyms, pseudonym2before.values()))
    for deleted_donor in deleted_donors:
        donor_diff_tables.append(rich.text.Text(f"Donor {deleted_donor.pseudonym} deleted", style="red"))

    return tuple(updated_donors), deleted_donors, tuple(donor_diff_tables)


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
def populate(ctx: click.Context, submission_id: str, metadata_path: str, confirm: bool, ignore_field: list[str]):  # noqa: C901
    """Populate the submission database from a metadata JSON file."""
    log.debug("Ignored fields for populate: %s", ignore_field)

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

    with open(metadata_path, encoding="utf-8") as metadata_file:
        metadata = GrzSubmissionMetadata.model_validate_json(metadata_file.read())

    changes = _diff_metadata(submission, metadata, set(ignore_field))

    # consent records
    updated_donors, deleted_donors, donor_diff_tables = _diff_donors(
        donors_in_db=db_service.get_donors(submission_id=Submission.id), metadata=metadata
    )

    if (not changes) and (not updated_donors) and (not deleted_donors):
        console_err.print("[green]Database is already up to date with the provided metadata![/green]")
        ctx.exit()

    diff_table: rich.console.RenderableType
    if changes:
        diff_table = rich.table.Table(title="Submission Metadata")
        diff_table.add_column("Key")
        diff_table.add_column("Before")
        diff_table.add_column("After")
        for key, before, after in sorted(changes, key=itemgetter(0)):
            diff_table.add_row(key, str(before) if before is not None else _TEXT_MISSING, str(after))
    else:
        diff_table = rich.padding.Padding(rich.text.Text("No changes to submission-level metadata."), pad=(0, 0, 1, 0))

    panel = rich.panel.Panel.fit(rich.console.Group(diff_table, *donor_diff_tables, fit=True), title="Pending Changes")
    console.print(panel)

    if not confirm or click.confirm(
        "Are you sure you want to commit these changes to the database?",
        default=False,
        show_default=True,
    ):
        for key, _before, after in changes:
            _ = db_service.modify_submission(submission_id, key=key, value=after)
        for updated_donor in updated_donors:
            _ = db_service.add_donor(updated_donor)
        for deleted_donor in deleted_donors:
            db_service.delete_donor(deleted_donor)
        console_err.print("[green]Database populated successfully.[/green]")


class QCStatus(StrEnum):
    PASS = "PASS"  # noqa: S105
    FAIL = "FAIL"
    TOO_LOW = "TOO LOW"


class QCReportRow(StrictBaseModel):
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
@click.argument("report_csv_path", metavar="path/to/report.csv", type=FILE_R_E)
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
        ("Genomic Study Type", "genomic_study_type"),
        ("Genomic Study Subtype", "genomic_study_subtype"),
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
