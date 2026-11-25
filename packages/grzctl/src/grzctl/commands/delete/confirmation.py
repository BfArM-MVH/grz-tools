"""``grzctl delete confirmation``: the GRZ's confirmation of a deletion to BfArM.

The BfArM SOP on deleting datasets in phase 0 (section IV.e) asks the data node to confirm a
deletion to BfArM's function mailbox by S/MIME encrypted e-mail, naming the tanGs. The text is
German, as the rest of that correspondence.
"""

import datetime
from pathlib import Path

import click
from grz_db.models.submission import SubmissionDb
from grz_pydantic_models.common import as_aware_datetime
from grz_pydantic_models.dates import date_to_quarter_year

from ...models.config import GrzctlConfig
from ...models.deletion_request import DeletionRequest
from .. import grzctl_configuration
from ..db.cli import get_submission_db_instance
from .run import _load_request

BACKUP_DELETION_DAYS = 90
"""Days within which copies in backups must be deleted (BfArM SOP, section IV.d)."""

RECORD_RETENTION_YEARS = 3
"""Years for which the records of a deletion, its tanGs among them, may be kept (BfArM SOP, section IV.g)."""


def years_later(date: datetime.date, years: int) -> datetime.date:
    """The same day *years* later, or earlier if negative; 29 February becomes 28 February."""
    try:
        return date.replace(year=date.year + years)
    except ValueError:  # 29 February
        return date.replace(year=date.year + years, day=28)


def _german(date: datetime.date) -> str:
    return date.strftime("%d.%m.%Y")


def build_confirmation(request: DeletionRequest, db: SubmissionDb, grz_id: str, today: datetime.date) -> str:
    """Write the confirmation of a deletion request to BfArM.

    :param request: The deletion request.
    :param db: The submission database.
    :param grz_id: ID of this GRZ.
    :param today: The date to judge whether the request's records may still be kept.
    :returns: The text of the e-mail, the subject line first.
    :raises click.ClickException: if a submission of the request is not deleted yet, or the
        records of the request are purged already.
    """
    if years_later(request.instructed_at, RECORD_RETENTION_YEARS) < today:
        raise click.ClickException(
            f"The request was instructed on {request.instructed_at}, more than {RECORD_RETENTION_YEARS} years ago. "
            "Its records are purged, so its tanGs can no longer be told from tanGs never present."
        )

    deleted_on: dict[str, datetime.date] = {}
    not_deleted = []
    for tan_g in request.tan_g:
        submission = db.get_submission_by_tan_g(tan_g)
        if submission is None:
            continue
        if submission.deleted_at is None:
            not_deleted.append(submission.id)
        else:
            deleted_on[tan_g] = as_aware_datetime(submission.deleted_at).astimezone(datetime.UTC).date()
    if not_deleted:
        raise click.ClickException(
            f"Not all datasets are deleted yet: {', '.join(not_deleted)}. Run grzctl delete run first."
        )

    intro = "hiermit bestätigen wir die Löschung der Datensätze gemäß" if deleted_on else "zu"
    lines = [
        "Betreff: Löschung von Datensätzen im genomDE MV",
        "",
        "Sehr geehrte Damen und Herren,",
        "",
        f"{intro} Ihrer Anweisung vom {_german(request.instructed_at)} (Antrag auf Löschung des "
        f"Leistungserbringers {request.submitter_id} vom {_german(request.requested_at)}):",
        "",
    ]
    for tan_g in request.tan_g:
        result = (
            f"gelöscht am {_german(deleted_on[tan_g])}" if tan_g in deleted_on else f"liegt beim GRZ {grz_id} nicht vor"
        )
        lines.append(f"  {tan_g}  {result}")

    if deleted_on:
        last_deletion = max(deleted_on.values())
        quarters = sorted({date_to_quarter_year(date)[::-1] for date in deleted_on.values()})
        backups_by = last_deletion + datetime.timedelta(days=BACKUP_DELETION_DAYS)
        lines += [
            "",
            f"Gelöscht wurden beim GRZ {grz_id} die genomischen Daten und Metadaten in Archiv und Inbox sowie die "
            "personenbezogenen Daten in der Datenbank. Ausgenommen sind die technischen Log-Dateien der "
            "Datenübermittlung.",
            "",
            "Kopien in Backups werden gemäß unseren Sicherungs- und Aufbewahrungszyklen spätestens bis zum "
            f"{_german(backups_by)} gelöscht. Bis dahin werden sie nicht weiter verarbeitet und nach einer "
            "Wiederherstellung unverzüglich erneut gelöscht.",
            "",
            "Die Löschung wird im Quartalsbericht für "
            + " und ".join(f"Q{quarter}/{year}" for year, quarter in quarters)
            + " in der Spalte number_of_deletions aufgeführt. Die tanGs werden zur Dokumentation bis längstens "
            f"{_german(years_later(last_deletion, RECORD_RETENTION_YEARS))} aufbewahrt.",
        ]
    return "\n".join(lines) + "\n"


@click.command("confirmation")
@grzctl_configuration
@click.option(
    "--data-file",
    "data_file",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="The deletion request, as given to grzctl delete run.",
)
def confirmation(configuration: GrzctlConfig, data_file: Path):
    """Print the confirmation of a deletion request to BfArM.

    Send it S/MIME encrypted to PlattformMvGenomseq@bfarm.de. It names the tanGs in full.
    Refuses while a submission of the request is not deleted yet.
    """
    request = _load_request(data_file)
    db = get_submission_db_instance(db_url=configuration.db.database_url)
    click.echo(build_confirmation(request, db, configuration.identifiers.grz, datetime.date.today()), nl=False)
