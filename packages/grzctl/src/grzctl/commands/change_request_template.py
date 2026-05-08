"""YAML template for change-request data.

Provides:

* ``render_yaml_template`` — pure helper that builds a fill-in YAML template
  for the given change-request type. Used by the ``db submission
  change-request`` error paths to print the expected schema when validation
  fails.
* ``change_request_template`` — top-level Click command that prints the same
  template on demand. Lives outside the ``db`` group so it does not require
  ``--config-file`` or DB credentials.
"""

from __future__ import annotations

import click
import rich.console
from grz_db.models.submission import ChangeRequestEnum

console_err = rich.console.Console(stderr=True)


_TYPE_GUIDANCE: dict[ChangeRequestEnum, str] = {
    ChangeRequestEnum.DELETE: "describe what to delete and why",
    ChangeRequestEnum.MODIFY: "describe the field/value to modify",
    ChangeRequestEnum.TRANSFER: "describe the transfer destination",
}


def render_yaml_template(change: ChangeRequestEnum) -> str:
    """Build a fill-in YAML template for a change-request type.

    Each placeholder embeds ``<FILL IN`` so the model validator rejects
    templates submitted unchanged.
    """
    guidance = _TYPE_GUIDANCE.get(change, "describe the requested change")
    return (
        "# Full name of the person requesting the change\n"
        "requester_name: <FILL IN requester name>\n"
        "\n"
        "# Email address of the requester\n"
        "requester_email: <FILL IN email address>\n"
        "\n"
        "# Date the change was requested (YYYY-MM-DD)\n"
        "requested_at: YYYY-MM-DD\n"
        "\n"
        "# Verbatim email text from the requester (multi-line allowed).\n"
        "# Required unless --raw-content (e.g. PDF) is supplied via the CLI.\n"
        "request_email_content: |\n"
        f"  <FILL IN paste verbatim email content here — {guidance}>\n"
        "\n"
        "# Optional type-specific extras (free-form key/value pairs)\n"
        "data: {}\n"
    )


@click.command("change-request-template")
@click.argument(
    "change_str",
    metavar="CHANGE",
    type=click.Choice(ChangeRequestEnum.list(), case_sensitive=False),
)
def change_request_template(change_str: str):
    """Print a YAML data template for the given change-request type."""
    change_request_enum = ChangeRequestEnum(change_str)
    click.echo(render_yaml_template(change_request_enum))
