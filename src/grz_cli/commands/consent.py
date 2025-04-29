"""Command for determining whether a submission is consented for research."""

import json
import logging
import sys
from pathlib import Path

import click

from ..workers.submission import SubmissionMetadata
from .common import submission_dir

log = logging.getLogger(__name__)

MDAT_WISSENSCHAFTLICH_NUTZEN_EU_DSGVO_NIVEAU = "2.16.840.1.113883.3.1937.777.24.5.3.8"
FHIR_PROVISION_PERMIT = "permit"
FHIR_PROVISION_DENY = "deny"


@click.command()
@submission_dir
def consent(submission_dir):
    """
    Check if a submission is consented for research.

    Returns 'true' if consented, 'false' if not.
    A submission is considered consented if all donors have consented for research, that is
    the FHIR MII IG Consent profiles all have a "permit" provision for code 2.16.840.1.113883.3.1937.777.24.5.3.8
    """
    metadata = SubmissionMetadata(Path(submission_dir) / "metadata" / "metadata.json").content

    consents = {donor.donor_pseudonym: False for donor in metadata.donors}
    for donor in metadata.donors:
        for research_consent in donor.research_consents:
            # TODO rename scope
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

    if all(consents.values()):
        click.echo("true")
    else:
        click.echo("false")


def _check_nested_provisions(provisions):
    for provision in provisions:
        if provision.get("type") == FHIR_PROVISION_PERMIT:
            for code in provision.get("code"):
                # TODO check system code as well
                # TODO check permit period start and end as well
                if code == MDAT_WISSENSCHAFTLICH_NUTZEN_EU_DSGVO_NIVEAU:
                    return True
    return False
