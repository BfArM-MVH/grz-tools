"""``grzctl delete purge-records``: delete the records of deletions once they may no longer be kept.

The BfArM SOP on deleting datasets in phase 0 (section IV.g) allows the records of a deletion to
be kept for at most three years after it was completed.
"""

import datetime

import click
import rich.console
import rich.table
from grz_db.models.submission import ChangeRequestEnum, Submission
from grz_pydantic_models.common import as_aware_datetime

from ...models.config import GrzctlConfig
from .. import grzctl_configuration
from ..db.cli import get_submission_db_instance
from .confirmation import RECORD_RETENTION_YEARS, years_later
from .resolve import new_table

console = rich.console.Console()
console_err = rich.console.Console(stderr=True)


def _deleted_on(submission: Submission) -> datetime.date:
    return as_aware_datetime(submission.deleted_at).astimezone(datetime.UTC).date()  # type: ignore[arg-type]


def _request(submission: Submission) -> tuple[str, str]:
    """Date of BfArM's instruction and the LE of the submission's deletion request, to find its e-mails by."""
    requests = [change for change in submission.changes if change.change == ChangeRequestEnum.DELETE]
    data = (requests[-1].data or {}) if requests else {}
    return data.get("instructed_at", "unknown"), data.get("submitter_id") or submission.submitter_id or "unknown"


def records_due(db_url: str, due_before: datetime.date) -> list[Submission]:
    """The deleted submissions whose records must be purged by *due_before*: deleted at least three years before it.

    :param db_url: URL of the submission database.
    :param due_before: The last day the records may be due.
    """
    last_deletion_day = years_later(due_before, -RECORD_RETENTION_YEARS)
    deleted_before = datetime.datetime.combine(
        last_deletion_day + datetime.timedelta(days=1), datetime.time(), datetime.UTC
    )
    return list(get_submission_db_instance(db_url=db_url).list_deletion_records(deleted_before))


@click.command("purge-records")
@grzctl_configuration
@click.option(
    "--due-before",
    required=True,
    type=click.DateTime(formats=["%Y-%m-%d"]),
    help="Purge the records of every deletion whose three years end by this date (YYYY-MM-DD).",
)
@click.option("--dry-run", is_flag=True, help="Show whose records would be purged, but purge nothing.")
@click.option("--yes-i-really-mean-it", is_flag=True, help="Purge without asking for confirmation.")
def purge_records(
    configuration: GrzctlConfig, due_before: datetime.datetime, dry_run: bool, yes_i_really_mean_it: bool
):
    """Purge the records of deletions that may be kept no longer.

    The records of a deletion, its tanG and its change requests with the request form, may be
    kept for three years after it. Purges them for every deletion whose three years end by
    --due-before, and lists the e-mails of each request, to delete from the mailbox too.
    """
    due = due_before.date()
    submissions = records_due(configuration.db.database_url, due)
    if not submissions:
        console_err.print(f"[green]No records are due by {due}.[/green]")
        return

    table = new_table("Submission", "Deleted on", "Due on", "LE", "BfArM instruction of")
    requests = {submission.id: _request(submission) for submission in submissions}
    for submission in submissions:
        instructed_at, submitter_id = requests[submission.id]
        deleted_on = _deleted_on(submission)
        table.add_row(
            submission.id,
            deleted_on.isoformat(),
            years_later(deleted_on, RECORD_RETENTION_YEARS).isoformat(),
            submitter_id,
            instructed_at,
        )
    console.print(table)

    if dry_run:
        console_err.print("[green]Dry run: nothing was purged.[/green]")
        return
    if not yes_i_really_mean_it and not click.confirm(
        f"Purge the records of {len(submissions)} deletion(s) for good?", default=False, show_default=True
    ):
        console_err.print("Nothing was purged.")
        raise click.Abort()

    db = get_submission_db_instance(db_url=configuration.db.database_url, author=configuration.db.signing_author)
    for submission in submissions:
        db.purge_deletion_records(submission.id)
    console_err.print(f"[green]Purged the records of {len(submissions)} deletion(s).[/green]")
    console_err.print("Delete the e-mails of these requests from the mailbox too (BfArM instruction of, LE):")
    for instructed_at, submitter_id in sorted(set(requests.values())):
        console_err.print(f"  {instructed_at}  {submitter_id}")
