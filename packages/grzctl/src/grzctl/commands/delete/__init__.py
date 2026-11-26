"""``grzctl delete``: carry out a BfArM request to delete submissions."""

import click

from .check import check
from .confirmation import confirmation
from .mark_done import mark_done
from .purge_records import purge_records
from .run import run
from .template import template


@click.group()
def delete():
    """Delete submissions on a BfArM deletion request.

    Fill in the request from grzctl delete template, carry out the deletion with
    grzctl delete run, and confirm it to BfArM with the text of grzctl delete confirmation.
    """


delete.add_command(template)
delete.add_command(run)
delete.add_command(confirmation)
delete.add_command(check)
delete.add_command(purge_records)
delete.add_command(mark_done)
