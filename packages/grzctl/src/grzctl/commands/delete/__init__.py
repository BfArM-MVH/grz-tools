"""``grzctl delete``: carry out a BfArM request to delete submissions."""

import click

from .template import template


@click.group()
def delete():
    """Delete submissions on a BfArM deletion request.

    Fill in the request from grzctl delete template, then carry out the deletion with
    grzctl delete run.
    """


delete.add_command(template)
