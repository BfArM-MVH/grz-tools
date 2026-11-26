"""``grzctl delete mark-done``: record a deletion carried out by hand before ``grzctl delete`` existed."""

import datetime
from dataclasses import asdict

import click
import rich.console
from grz_db.models.submission import ChangeRequestEnum, SubmissionStateEnum

from ... import get_recorded_versions
from ...models.config import GrzctlConfig
from .. import grzctl_configuration
from ..db.cli import get_submission_db_instance
from .execute import describe_holdings, leftovers, storage_targets

console_err = rich.console.Console(stderr=True)


@click.command("mark-done")
@grzctl_configuration
@click.argument("submission_id")
@click.option(
    "--deleted-on",
    required=True,
    type=click.DateTime(formats=["%Y-%m-%d"]),
    help="Date the submission was deleted by hand (YYYY-MM-DD).",
)
def mark_done(configuration: GrzctlConfig, submission_id: str, deleted_on: datetime.datetime):
    """Record that a submission was deleted by hand, before grzctl delete existed.

    The submission must have a Delete change request, and its archives and inbox must hold
    nothing of it any more. Records the states Deleting and Deleted on --deleted-on, so that the
    reports count the deletion in its quarter, and purges the submission's personal data.
    """
    deleted_date = deleted_on.date()
    today = datetime.date.today()
    if deleted_date > today:
        raise click.BadParameter(f"{deleted_date} is in the future.", param_hint="--deleted-on")

    db = get_submission_db_instance(db_url=configuration.db.database_url)
    submission = db.get_submission(submission_id)
    if submission is None:
        raise click.ClickException(f"Submission {submission_id} is not in the database.")
    if submission.deletion_started:
        raise click.ClickException(f"The deletion of {submission_id} is recorded already.")
    changes = next((s.changes for s in db.list_change_requests() if s.id == submission_id), [])
    if not any(change.change == ChangeRequestEnum.DELETE for change in changes):
        raise click.ClickException(
            f"{submission_id} has no Delete change request, so nothing records that its deletion was requested."
        )

    submitter_id = submission.submitter_id or submission_id.split("_", maxsplit=1)[0]
    targets = storage_targets(configuration, db, submission_id, submitter_id)
    if left := leftovers(targets, submission_id):
        raise click.ClickException(f"{describe_holdings(left)}. Delete {submission_id} with grzctl delete run instead.")

    signing_db = get_submission_db_instance(
        db_url=configuration.db.database_url, author=configuration.db.signing_author
    )
    at = datetime.datetime.combine(deleted_date, datetime.time(), datetime.UTC)
    signing_db.update_submission_state(
        submission_id, SubmissionStateEnum.DELETING, timestamp=at, grzctl_versions=get_recorded_versions()
    )
    purge = signing_db.purge_personal_data(submission_id)
    signing_db.update_submission_state(
        submission_id,
        SubmissionStateEnum.DELETED,
        timestamp=at,
        data={"deleted_by_hand": True, "recorded_on": today.isoformat(), "database": asdict(purge)},
        grzctl_versions=get_recorded_versions(),
    )
    console_err.print(f"[green]Recorded {submission_id} as deleted on {deleted_date}.[/green]")
