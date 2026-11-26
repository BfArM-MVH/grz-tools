"""``grzctl delete template``: print the YAML template of a BfArM deletion request."""

import click

from ...models.deletion_request import DATE_PLACEHOLDER, PLACEHOLDER_MARKER

DELETION_REQUEST_TEMPLATE = f"""\
# Antrag auf Löschung von Daten (BfArM request form, version 1.1).
# Copy the values from the form, then check them with grzctl delete run --dry-run.

# Date of BfArM's instruction e-mail from PlattformMvGenomseq@bfarm.de, which is not on the form
instructed_at: {DATE_PLACEHOLDER}
# Meldender LE: Institutionskennzeichen (9 digits)
submitter_id: "{PLACEHOLDER_MARKER} LE-ID>"
# Kontaktperson
requester_name: "{PLACEHOLDER_MARKER} name>"
requester_email: "{PLACEHOLDER_MARKER} e-mail address>"
# Datum des Versands des Antrags
requested_at: {DATE_PLACEHOLDER}
# Begründung für Löschung
reason: |
  {PLACEHOLDER_MARKER} reason>
# Beteiligtes GRZ (Knoten-ID)
data_node_id: "{PLACEHOLDER_MARKER} GRZ-ID>"
# tanG of each genomic dataset to delete, one per line
tan_g:
  - "{PLACEHOLDER_MARKER} tanG>"
# Kommentare (optional)
comment: ""
"""


@click.command("template")
def template():
    """Print the YAML template of a BfArM deletion request.

    Redirect it to a file, fill it in from the request form and check it with
    grzctl delete run --dry-run: grzctl delete template > antrag.yaml
    """
    click.echo(DELETION_REQUEST_TEMPLATE, nl=False)
