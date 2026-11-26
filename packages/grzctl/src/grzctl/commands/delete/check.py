"""``grzctl delete check``: find deleted submissions that a restore from a backup brought back.

The BfArM SOP on deleting datasets in phase 0 (section IV.d) asks the data node to delete again,
without delay, what a restore from a backup brings back of a deleted submission.
"""

from dataclasses import dataclass
from enum import StrEnum

import click
import rich.console
import rich.table
from grz_db.models.submission import SubmissionDb, SubmissionStateEnum

from ... import get_recorded_versions
from ...models.config import GrzctlConfig
from .. import grzctl_configuration
from ..db.cli import get_submission_db_instance
from .execute import (
    StorageTargets,
    archive_targets,
    describe_holdings,
    finish_deletion,
    leftovers,
    storage_targets,
)
from .resolve import new_table
from .storage import tombstoned_submissions

console = rich.console.Console()
console_err = rich.console.Console(stderr=True)


class Finding(StrEnum):
    """What ``grzctl delete check`` can find of a deleted submission."""

    STORAGE_RESTORED = "storage restored"
    DATABASE_RESTORED = "database restored"
    INTERRUPTED = "deletion interrupted"
    NOT_IN_DATABASE = "not in the database"
    UNCHECKED = "not checked"


REPAIRABLE = (Finding.STORAGE_RESTORED, Finding.DATABASE_RESTORED)
"""What ``grzctl delete check --repair`` deletes again."""


@dataclass(frozen=True)
class CheckResult:
    """A deleted submission that is not entirely deleted.

    :param submission_id: ID of the submission.
    :param finding: What is wrong.
    :param detail: What was found, in words.
    :param targets: Its buckets, when they are known.
    """

    submission_id: str
    finding: Finding
    detail: str
    targets: StorageTargets | None = None


def _submitter_id(submission_id: str) -> str:
    return submission_id.split("_", maxsplit=1)[0]


def find_restored(configuration: GrzctlConfig, db: SubmissionDb) -> list[CheckResult]:
    """Find what is not entirely deleted of the submissions whose deletion started.

    Checks each submission the database records as deleted for objects in its buckets, and each
    tombstone in the archives for a deletion the database records.

    :param configuration: The grzctl configuration.
    :param db: The submission database.
    :returns: One result for each submission that is not entirely deleted.
    """
    results = []
    deletions = {submission.id: submission for submission in db.list_deletions()}
    for submission_id, submission in deletions.items():
        if not submission.is_deleted:
            results.append(
                CheckResult(
                    submission_id,
                    Finding.INTERRUPTED,
                    "its deletion stopped before it finished; run grzctl delete run with its request again",
                )
            )
            continue
        try:
            targets = storage_targets(configuration, db, submission_id, _submitter_id(submission_id))
        except click.ClickException as e:
            results.append(CheckResult(submission_id, Finding.UNCHECKED, e.message))
            continue
        if held := leftovers(targets, submission_id):
            detail = describe_holdings(held, "the {location} holds {objects} objects again")
            results.append(CheckResult(submission_id, Finding.STORAGE_RESTORED, detail, targets))

    for location, s3_options in archive_targets(configuration).items():
        for submission_id in sorted(tombstoned_submissions(s3_options) - deletions.keys()):
            detail = f"the {location} holds its tombstone, but the database records no deletion"
            if db.get_submission(submission_id) is None:
                results.append(CheckResult(submission_id, Finding.NOT_IN_DATABASE, detail))
                continue
            try:
                targets = storage_targets(configuration, db, submission_id, _submitter_id(submission_id))
            except click.ClickException as e:
                results.append(CheckResult(submission_id, Finding.UNCHECKED, f"{detail}; {e.message}"))
                continue
            results.append(CheckResult(submission_id, Finding.DATABASE_RESTORED, detail, targets))
    return results


def repair(db: SubmissionDb, result: CheckResult) -> None:
    """Delete again what a restore brought back of a deleted submission, and scrub it again.

    :param db: The submission database, with an author to sign the entries.
    :param result: A result of :func:`find_restored` with a finding in :data:`REPAIRABLE`.
    """
    if result.finding not in REPAIRABLE or result.targets is None:
        raise ValueError(f"{result.submission_id}: '{result.finding}' cannot be repaired.")
    if result.finding == Finding.DATABASE_RESTORED:
        db.update_submission_state(
            result.submission_id, SubmissionStateEnum.DELETING, grzctl_versions=get_recorded_versions()
        )
        note = {
            "repeated": "database restored from a backup; the request is on record in the mailbox, "
            "its change request was lost with the restore"
        }
    else:
        note = {"repeated": "storage restored from a backup"}
    finish_deletion(db, result.submission_id, result.targets, note=note)


@click.command("check")
@grzctl_configuration
@click.option("--repair", "repair_found", is_flag=True, help="Delete again what a restore brought back.")
def check(configuration: GrzctlConfig, repair_found: bool):
    """Find deleted submissions that a restore from a backup brought back.

    Run it after every restore of the database or of a bucket. Finds data in the archives or the
    inbox of a submission recorded as deleted, a tombstone in an archive of a submission the
    database records no deletion of, and deletions that stopped before they finished. With
    --repair, deletes the data again and purges the submission's personal data again. Exits with
    an error while anything is left that it did not repair.
    """
    db = get_submission_db_instance(db_url=configuration.db.database_url)
    results = find_restored(configuration, db)
    if not results:
        console_err.print("[green]Every deleted submission is deleted entirely.[/green]")
        return

    table = new_table("Submission", "Finding", "Detail")
    for result in results:
        table.add_row(result.submission_id, result.finding.value, result.detail)
    console.print(table)

    left = results
    if repair_found:
        signing_db = get_submission_db_instance(
            db_url=configuration.db.database_url, author=configuration.db.signing_author
        )
        for result in results:
            if result.finding not in REPAIRABLE:
                continue
            try:
                repair(signing_db, result)
            except Exception as e:
                raise click.ClickException(f"Deleting {result.submission_id} again failed: {e}") from e
            console_err.print(f"[green]Deleted {result.submission_id} again.[/green]")
        left = [result for result in results if result.finding not in REPAIRABLE]
    if left:
        hint = "" if repair_found else " Run grzctl delete check --repair to delete again what a restore brought back."
        raise click.ClickException(f"{len(left)} submission(s) are not entirely deleted.{hint}")
