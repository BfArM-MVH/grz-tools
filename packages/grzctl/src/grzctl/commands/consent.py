"""Command for determining whether a submission is consented for research."""

import enum
import json
import logging
import sys
import typing
from pathlib import Path
from typing import Any

import click
import rich.console
import rich.table
import rich.text
from grz_common.cli import output_json, show_details, submission_dir
from grz_common.workers.submission import GrzSubmissionMetadata, SubmissionMetadata

log = logging.getLogger(__name__)

MDAT_WISSENSCHAFTLICH_NUTZEN_EU_DSGVO_NIVEAU = "2.16.840.1.113883.3.1937.777.24.5.3.8"


class FhirProvision(enum.StrEnum):
    """Possible FHIR Provision options."""

    PERMIT = "permit"
    DENY = "deny"


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


def _gather_consent_information(metadata: GrzSubmissionMetadata) -> dict[str, bool]:
    consents = {donor.donor_pseudonym: False for donor in metadata.donors}
    for donor in metadata.donors:
        for research_consent in donor.research_consents:
            mii_consent = research_consent.scope
            if isinstance(mii_consent, str):
                mii_consent = json.loads(mii_consent)
            mii_consent = typing.cast(dict[str, Any], mii_consent)

            if top_level_provision := mii_consent.get("provision"):
                if top_level_provision.get("type") != FhirProvision.DENY:
                    sys.exit(
                        f"The root provision type must be deny, not {top_level_provision.get('type')}, "
                        f"since the profile follows an opt-in consent scheme. "
                        f"Explicit opt-in consents must be made via nested provisions."
                    )
                else:
                    nested_provisions = top_level_provision.get("provision")
                    consents[donor.donor_pseudonym] = _check_nested_provisions(nested_provisions)

    return consents


def _check_nested_provisions(provisions: list[dict[str, Any]]) -> bool:
    for provision in provisions:
        if provision.get("type") == FhirProvision.PERMIT:
            for codeable_concept in provision.get("code", []):
                for coding in codeable_concept.get("coding", []):
                    code = coding.get("code")
                    if isinstance(code, str):
                        if code == MDAT_WISSENSCHAFTLICH_NUTZEN_EU_DSGVO_NIVEAU:
                            return True
                    elif isinstance(code, dict):
                        if (value := code.get("value")) and value == MDAT_WISSENSCHAFTLICH_NUTZEN_EU_DSGVO_NIVEAU:
                            return True
                    else:
                        raise ValueError(code, f"Expected str or dict, got {type(code)}")
    return False
