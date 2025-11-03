"""Command for determining whether a submission is consented for research."""

import datetime
import json
import logging
import sys
from pathlib import Path

import click
import rich.console
import rich.table
import rich.text
from grz_common.cli import FILE_R_E, output_json, show_details, submission_dir
from grz_common.workers.submission import GrzSubmissionMetadata, SubmissionMetadata

log = logging.getLogger(__name__)


@click.command()
@submission_dir
@click.option(
    "--metadata-file",
    metavar="PATH",
    type=FILE_R_E,
    required=False,
    help="Direct path to the metadata.json file.",
)
@output_json
@show_details
@click.option("--date", help="date for which to check consent validity in ISO format (default: today)")
def consent(submission_dir, metadata_file, output_json, show_details, date):
    """
    Check if a submission is consented for research.

    Returns 'true' if consented, 'false' if not.
    """
    in_legacy_mode = submission_dir is not None
    in_flexible_mode = metadata_file is not None

    if in_legacy_mode and in_flexible_mode:
        raise click.UsageError("'--submission-dir' is mutually exclusive with '--metadata-file'.")

    if in_legacy_mode:
        metadata_path = Path(submission_dir) / "metadata" / "metadata.json"
    elif in_flexible_mode:
        metadata_path = Path(metadata_file)
    else:
        raise click.UsageError("You must specify either '--submission-dir' or '--metadata-file'.")

    metadata = SubmissionMetadata(metadata_path).content

    date = datetime.date.today() if date is None else datetime.date.fromisoformat(date)
    consents = _gather_consent_information(metadata, date)
    overall_consent = metadata.consents_to_research(date)

    match output_json, show_details:
        case True, True:
            json.dump(consents, sys.stdout)
        case True, False:
            json.dump(overall_consent, sys.stdout)
        case False, True:
            _print_rich_table(consents)
        case False, False:
            click.echo(str(overall_consent).lower())


def _print_rich_table(consents: dict[str, bool]):
    console = rich.console.Console()
    table = rich.table.Table()
    table.add_column("Donor", no_wrap=True)
    table.add_column("Research Consent", no_wrap=True)
    for donor_pseudonym, consent_value in consents.items():
        research_consent = rich.text.Text(
            "True" if consent_value else "False",
            style="green" if consent_value else "red",
        )
        table.add_row(
            donor_pseudonym,
            research_consent,
        )
    console.print(table)


def _gather_consent_information(metadata: GrzSubmissionMetadata, date: datetime.date) -> dict[str, bool]:
    consents = {donor.donor_pseudonym: False for donor in metadata.donors}
    for donor in metadata.donors:
        consents[donor.donor_pseudonym] = donor.consents_to_research(date)

    return consents
