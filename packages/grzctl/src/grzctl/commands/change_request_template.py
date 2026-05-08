"""YAML template for change-request data.

Provides:

* ``render_yaml_template`` — pure helper that builds a fill-in YAML template
  from a Pydantic model. Used by the ``db submission change-request`` error
  paths to print the expected schema when validation fails.
* ``change_request_template`` — top-level Click command that prints the same
  template on demand. Lives outside the ``db`` group so it does not require
  ``--config-file`` or DB credentials.
"""

from __future__ import annotations

import datetime
from typing import Any

import click
import rich.console
from grz_db.models.submission import CHANGE_REQUEST_DATA_SCHEMAS, ChangeRequestEnum
from pydantic import BaseModel
from pydantic.fields import FieldInfo

console_err = rich.console.Console(stderr=True)


def render_yaml_template(model_cls: type[BaseModel]) -> str:
    """Build a fill-in YAML template from a Pydantic model's fields.

    Each field becomes a ``# <description>`` comment followed by ``name: <placeholder>``.
    Placeholder shape is driven by the field's annotation (e.g. multi-line block
    scalar for ``*_content`` strings, ``YYYY-MM-DD`` for dates).
    """
    lines: list[str] = []
    for name, info in model_cls.model_fields.items():
        if info.description:
            lines.append(f"# {info.description}")
        lines.append(f"{name}: {_placeholder_for(name, info)}")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def _placeholder_for(name: str, info: FieldInfo) -> str:
    annotation: Any = info.annotation
    if annotation is datetime.date:
        # Not a parseable date — fails Pydantic validation if left as-is.
        return "YYYY-MM-DD"
    if annotation is str:
        # Each placeholder embeds "<FILL IN" so the schema's model validator
        # rejects them if the user submits the template unchanged.
        # No colons inside angle brackets — YAML parses '<FILL IN: x>' as a nested mapping.
        if name.endswith("_content"):
            return "|\n  <FILL IN paste verbatim email content here (multi-line allowed)>"
        if "email" in name:
            return "<FILL IN email address>"
        return f"<FILL IN {name.replace('_', ' ')}>"
    return "null"


@click.command("change-request-template")
@click.argument(
    "change_str",
    metavar="CHANGE",
    type=click.Choice(ChangeRequestEnum.list(), case_sensitive=False),
)
def change_request_template(change_str: str):
    """Print a YAML data template for the given change-request type."""
    change_request_enum = ChangeRequestEnum(change_str)
    schema_cls = CHANGE_REQUEST_DATA_SCHEMAS.get(change_request_enum)
    if schema_cls is None:
        console_err.print(
            f"[yellow]No data schema is defined for change type '{change_request_enum.value}'.[/yellow]"
        )
        return
    click.echo(render_yaml_template(schema_cls))
