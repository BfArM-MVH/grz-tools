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
from grz_common.cli import output_json, show_details, submission_dir
from grz_common.workers.submission import GrzSubmissionMetadata, SubmissionMetadata

log = logging.getLogger(__name__)


@click.command()
@submission_dir
@output_json
@show_details
@click.option("--date", help="date for which to check consent validity in ISO format (default: today)")
def consent(submission_dir, output_json, show_details, date):
    """
    Check if a submission is consented for research.

    Returns 'true' if consented, 'false' if not.
    A submission is considered consented if all donors have consented for research, that is
    the FHIR MII IG Consent profiles all have a "permit" provision for code 2.16.840.1.113883.3.1937.777.24.5.3.8
    """
    metadata = SubmissionMetadata(Path(submission_dir) / "metadata" / "metadata.json").content

    date = datetime.date.today() if date is None else datetime.date.fromisoformat(date)
    consents = _gather_consent_information(metadata, date)
    overall_consent = _submission_has_research_consent(consents)

    match output_json, show_details:
        case True, True:
            json.dump(consents, sys.stdout)
        case True, False:
            json.dump(overall_consent, sys.stdout)
        case False, True:
            _print_rich_table(consents)
        case False, False:
            click.echo(str(overall_consent).lower())


def _submission_has_research_consent(consents):
    return all(consents.values())


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
