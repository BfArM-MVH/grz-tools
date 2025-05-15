"""Command for determining whether a submission is consented for research."""

import json
import logging
import sys
from pathlib import Path

import click
import rich.console
import rich.table
import rich.text

from ..workers.submission import SubmissionMetadata
from .common import output_json, show_details, submission_dir

log = logging.getLogger(__name__)

MDAT_WISSENSCHAFTLICH_NUTZEN_EU_DSGVO_NIVEAU = "2.16.840.1.113883.3.1937.777.24.5.3.8"
FHIR_PROVISION_PERMIT = "permit"
FHIR_PROVISION_DENY = "deny"


@click.command()
@submission_dir
@output_json
@show_details
def consent(submission_dir, output_json, show_details):
    """
    Check if a submission is consented for research.

    Returns 'true' if consented, 'false' if not.
    A submission is considered consented if all donors have consented for research, that is
    the FHIR MII IG Consent profiles all have a "permit" provision for code 2.16.840.1.113883.3.1937.777.24.5.3.8
    """
    metadata = SubmissionMetadata(Path(submission_dir) / "metadata" / "metadata.json").content

    consents = _gather_consent_information(metadata)

    if not show_details:
        if output_json:
            click.echo(json.dumps(all(consents.values())))
        else:
            click.echo(str(all(consents.values())).lower())
    elif output_json:
        json.dump(consents, sys.stdout)
    else:
        _print_rich_table(consents)


def _print_rich_table(consents):
    console = rich.console.Console()
    table = rich.table.Table()
    table.add_column("Donor", no_wrap=True)
    table.add_column("Consent Type", no_wrap=True)
    table.add_column("Consent State", no_wrap=True)
    for donor_pseudonym, consent_value in consents.items():
        status_text = rich.text.Text(
            "True" if consent_value else "False",
            style="green" if consent_value else "red",
        )
        table.add_row(
            donor_pseudonym,
            "research consent",
            status_text,
        )
    console.print(table)


def _gather_consent_information(metadata):
    consents = {donor.donor_pseudonym: False for donor in metadata.donors}
    for donor in metadata.donors:
        for research_consent in donor.research_consents:
            # TODO rename scope to resource
            mii_consent = research_consent.scope
            if isinstance(mii_consent, str):
                mii_consent = json.loads(mii_consent)

            if top_level_provision := mii_consent.get("provision"):
                if top_level_provision.get("type") != FHIR_PROVISION_DENY:
                    sys.exit(
                        f"Top level provision type must be deny, not {top_level_provision.get('type')}. "
                        f"(Opt-in via nested provisions)"
                    )
                else:
                    nested_provisions = top_level_provision.get("provision")
                    consents[donor.donor_pseudonym] = _check_nested_provisions(nested_provisions)
            else:
                consents[donor.donor_pseudonym] = False

    return consents


def _check_nested_provisions(provisions):
    for provision in provisions:
        if provision.get("type") == FHIR_PROVISION_PERMIT:
            for code in provision.get("code"):
                # TODO check system code as well
                # TODO check permit period start and end as well
                for coding in code.get("coding"):
                    code_ = coding.get("code")
                    if code_ == MDAT_WISSENSCHAFTLICH_NUTZEN_EU_DSGVO_NIVEAU:
                        return True
    return False
