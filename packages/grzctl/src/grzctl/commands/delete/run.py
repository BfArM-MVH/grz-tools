"""``grzctl delete run``: carry out a BfArM deletion request."""

from pathlib import Path

import click
import rich.console
import rich.table
import yaml

from ...models.config import GrzctlConfig
from ...models.deletion_request import DeletionRequest, deletion_request_problems
from .. import grzctl_configuration
from ..change_request import _read_raw_content
from ..db.cli import get_submission_db_instance
from .resolve import (
    DeletionOutcome,
    DeletionRequestRefusedError,
    ResolvedDeletionRequest,
    resolve_deletion_request,
    short_tan_g,
)

console = rich.console.Console()
console_err = rich.console.Console(stderr=True)


def _load_request(data_file: Path) -> DeletionRequest:
    """Read and check a deletion request, reporting every problem at once.

    :raises click.Abort: if the request has a problem.
    """
    try:
        raw = yaml.safe_load(data_file.read_text(encoding="utf-8"))
    except yaml.YAMLError as e:
        console_err.print(f"[red]Error: {data_file} is no valid YAML: {e}[/red]")
        raise click.Abort() from e
    if not isinstance(raw, dict):
        console_err.print(
            f"[red]Error: {data_file} does not hold a deletion request; see grzctl delete template.[/red]"
        )
        raise click.Abort()
    if problems := deletion_request_problems(raw):
        console_err.print("[red]Error: the deletion request needs attention:[/red]")
        for problem in problems:
            console_err.print(f"  [red]•[/red] {problem}")
        raise click.Abort()
    return DeletionRequest.model_validate(raw)


def _print_plan(request: DeletionRequest, resolved: ResolvedDeletionRequest, configuration: GrzctlConfig) -> None:
    submitter = configuration.leistungserbringer.get(request.submitter_id)
    alias = f" ({submitter.alias})" if submitter is not None and submitter.alias != request.submitter_id else ""
    console.print(
        f"Deletion request of {request.submitter_id}{alias}, sent {request.requested_at}, "
        f"instructed {request.instructed_at}"
    )
    table = rich.table.Table(box=None)
    for column in ("tanG", "Submission", "State", "Result"):
        table.add_column(column)
    for r in resolved.tan_gs:
        latest = r.submission.get_latest_state() if r.submission is not None else None
        result = r.outcome.value
        if r.outcome == DeletionOutcome.NOT_PRESENT:
            result = f"{result} at {configuration.identifiers.grz}"
        table.add_row(
            short_tan_g(r.tan_g),
            r.submission.id if r.submission is not None else "—",
            latest.state.value if latest is not None else "—",
            result,
        )
    console.print(table)
    for warning in resolved.case_warnings:
        console_err.print(f"[yellow]⚠ {warning}[/yellow]")
    if any(r.outcome == DeletionOutcome.NOT_PRESENT for r in resolved.tan_gs):
        console_err.print(
            "[yellow]A tanG that is not present has no submission in the database. A submission not yet in the "
            "database is deleted by hand, not by grzctl delete run.[/yellow]"
        )


@click.command("run")
@grzctl_configuration
@click.option(
    "--data-file",
    "data_file",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="The deletion request, filled in from grzctl delete template.",
)
@click.option(
    "--raw-content",
    "raw_content_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="The request form as the LE sent it (PDF), to keep with the deletion.",
)
@click.option("--dry-run", is_flag=True, help="Check the request and show what would be deleted, but delete nothing.")
def run(configuration: GrzctlConfig, data_file: Path, raw_content_path: Path, dry_run: bool):
    """Delete the submissions of a BfArM deletion request."""
    request = _load_request(data_file)
    _read_raw_content(raw_content_path)

    db = get_submission_db_instance(db_url=configuration.db.database_url)
    try:
        resolved = resolve_deletion_request(request, db, configuration.identifiers.grz)
    except DeletionRequestRefusedError as e:
        console_err.print("[red]Error: the deletion request cannot be carried out:[/red]")
        for problem in e.problems:
            console_err.print(f"  [red]•[/red] {problem}")
        raise click.Abort() from e

    _print_plan(request, resolved, configuration)
    if dry_run:
        console_err.print("[green]Dry run: nothing was deleted.[/green]")
        return
    raise click.ClickException("Deleting is not implemented yet; run with --dry-run.")
