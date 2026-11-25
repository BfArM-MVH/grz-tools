"""``grzctl delete run``: carry out a BfArM deletion request."""

from pathlib import Path

import click
import rich.console
import rich.table
import yaml
from grz_db.models.submission import RequestRawContentType, SubmissionDb

from ...models.config import GrzctlConfig
from ...models.deletion_request import DeletionRequest, deletion_request_problems
from .. import grzctl_configuration
from ..change_request import _read_raw_content
from ..db.cli import get_submission_db_instance
from .execute import StorageTargets, delete_submission, storage_targets
from .resolve import (
    DeletionOutcome,
    DeletionRequestRefusedError,
    ResolvedDeletionRequest,
    resolve_deletion_request,
    short_tan_g,
)
from .storage import Holdings, archive_holdings, has_object_lock, inbox_holdings

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


def _size(num_bytes: int) -> str:
    size = float(num_bytes)
    for unit in ("B", "KB", "MB"):
        if size < 1000:
            return f"{size:.0f} {unit}" if unit == "B" else f"{size:.1f} {unit}"
        size /= 1000
    return f"{size:.1f} GB"


def _holdings_text(holdings: Holdings) -> str:
    return f"{holdings.objects} objects, {_size(holdings.bytes)}" if holdings.objects else "nothing"


def _print_storage(targets: dict[str, StorageTargets]) -> None:
    """Show what each bucket holds of each submission to delete."""
    table = rich.table.Table(box=None)
    for column in ("Submission", "Location", "To delete"):
        table.add_column(column)
    for submission_id, submission_targets in targets.items():
        for location, s3_options in submission_targets.archives.items():
            if (holdings := archive_holdings(s3_options, submission_id)).objects:
                table.add_row(submission_id, location, _holdings_text(holdings))
        location, s3_options = submission_targets.inbox
        table.add_row(submission_id, location, _holdings_text(inbox_holdings(s3_options, submission_id)))
    console.print(table)


def _locked_buckets(targets: dict[str, StorageTargets]) -> list[str]:
    """Each bucket among *targets* that has object lock enabled, which keeps its objects from being deleted."""
    buckets = {
        s3_options.bucket: (location, s3_options)
        for submission_targets in targets.values()
        for location, s3_options in submission_targets.buckets.items()
    }
    return [
        f"{location} (bucket {bucket})"
        for bucket, (location, s3_options) in buckets.items()
        if has_object_lock(s3_options)
    ]


def _resolve(request: DeletionRequest, db: SubmissionDb, grz_id: str) -> ResolvedDeletionRequest:
    """Resolve the request, or abort naming why it cannot be carried out."""
    try:
        return resolve_deletion_request(request, db, grz_id)
    except DeletionRequestRefusedError as e:
        console_err.print("[red]Error: the deletion request cannot be carried out:[/red]")
        for problem in e.problems:
            console_err.print(f"  [red]•[/red] {problem}")
        raise click.Abort() from e


def _delete_all(
    configuration: GrzctlConfig,
    resolved: ResolvedDeletionRequest,
    targets: dict[str, StorageTargets],
    request: DeletionRequest,
    request_form: tuple[bytes, RequestRawContentType],
) -> None:
    """Delete each submission to delete, stopping at the first that fails."""
    db = get_submission_db_instance(db_url=configuration.db.database_url, author=configuration.db.signing_author)
    for r in resolved.tan_gs:
        if r.submission is None or r.outcome not in (DeletionOutcome.DELETE, DeletionOutcome.RESUME):
            continue
        try:
            delete_submission(
                db,
                r.submission,
                targets[r.submission.id],
                request,
                request_form,
                record_request=r.outcome == DeletionOutcome.DELETE,
            )
        except Exception as e:
            raise click.ClickException(
                f"Deleting {r.submission.id} failed: {e}. It stays in state Deleting, and the submissions after it "
                "are not deleted yet. Run grzctl delete run again with the same request to continue."
            ) from e
        console_err.print(f"[green]Deleted {r.submission.id} (tanG {short_tan_g(r.tan_g)}).[/green]")


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
@click.option("--yes-i-really-mean-it", is_flag=True, help="Delete without asking for confirmation.")
def run(
    configuration: GrzctlConfig, data_file: Path, raw_content_path: Path, dry_run: bool, yes_i_really_mean_it: bool
):
    """Delete the submissions of a BfArM deletion request.

    Deletes each submission from the archives that hold it, except its log files, and from its
    inbox, and purges its personal data from the database. If deleting a submission fails, it
    stays in state Deleting; run the command again with the same request to continue.
    """
    request = _load_request(data_file)
    request_form, request_form_type = _read_raw_content(raw_content_path)
    if request_form is None or request_form_type is None:
        raise click.UsageError("--raw-content is required.")

    db = get_submission_db_instance(db_url=configuration.db.database_url)
    resolved = _resolve(request, db, configuration.identifiers.grz)
    _print_plan(request, resolved, configuration)

    targets = {s.id: storage_targets(configuration, db, s.id, request.submitter_id) for s in resolved.to_delete}
    if targets:
        _print_storage(targets)
    if locked := _locked_buckets(targets):
        console_err.print(
            f"[red]Error: object lock keeps objects from being deleted in {', '.join(locked)}. Nothing was deleted.[/red]"
        )
        raise click.Abort()

    if dry_run:
        console_err.print("[green]Dry run: nothing was deleted.[/green]")
        return
    if not targets:
        console_err.print("[green]Nothing to delete.[/green]")
        return
    if not yes_i_really_mean_it and not click.confirm(
        f"Delete {len(targets)} submission(s) for good?", default=False, show_default=True
    ):
        console_err.print("Nothing was deleted.")
        raise click.Abort()
    _delete_all(configuration, resolved, targets, request, (request_form, request_form_type))
