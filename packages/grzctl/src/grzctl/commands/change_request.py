"""Loading + validation for change-request input, and the offline
``change-request-validate`` command.

The helpers here perform *pure* validation of change-request data supplied on
the CLI (``--data`` / ``--data-file`` / ``--raw-content``): required audit
fields, ``<FILL IN`` placeholder rejection, and model validation. They touch no
database and need no configuration, so they are shared between:

* ``grzctl db ... submission change-request`` — validates, then writes, and
* ``grzctl change-request-validate`` — validates only, never connects to a DB.

Keeping validation config-free means it is impossible to accidentally write a
live change request while only meaning to check an input file.
"""

from __future__ import annotations

import json
from datetime import date
from pathlib import Path

import click
import rich.console
import yaml
from grz_db.models.submission import ChangeRequestEnum, ChangeRequestLogBase, RequestRawContentType
from pydantic import ValidationError

from .change_request_template import render_yaml_template

console_err = rich.console.Console(stderr=True)


_COLUMNAR_KEYS = {"requester_name", "requester_email", "requested_at", "request_email_content"}

# Placeholder tokens emitted by ``render_yaml_template`` (kept in sync with the template
# and the grz_db ``_TEMPLATE_PLACEHOLDER_MARKER`` model validator), so an unedited template
# is reported clearly rather than via pydantic's field-by-field coercion errors.
_PLACEHOLDER_MARKER = "<FILL IN"
_DATE_PLACEHOLDER = "YYYY-MM-DD"

_RAW_CONTENT_TYPE_BY_SUFFIX: dict[str, RequestRawContentType] = {
    ".pdf": RequestRawContentType.PDF,
    ".png": RequestRawContentType.PNG,
}


def _load_change_request_input(data_json: str | None, data_file: Path | None) -> dict | None:
    """Resolve the raw change-request input dict from inline JSON or a JSON/YAML file."""
    if data_json is not None and data_file is not None:
        console_err.print("[red]Error: --data and --data-file are mutually exclusive.[/red]")
        raise click.Abort()
    if data_json is not None:
        try:
            return json.loads(data_json)
        except json.JSONDecodeError as e:
            console_err.print(f"[red]Error: Invalid JSON string for --data: {e}[/red]")
            raise click.Abort() from e
    if data_file is None:
        return None
    text = data_file.read_text()
    # YAML is a superset of JSON, so a single YAML parse accepts .json, .yaml and
    # .yml alike — the extension is only a hint, not a hard requirement. This also
    # means a correctly-formed file loads even if it carries the "wrong" extension.
    try:
        raw = yaml.safe_load(text)
    except yaml.YAMLError as e:
        console_err.print(f"[red]Error: could not parse '{data_file}' as YAML or JSON.[/red]")
        console_err.print(f"[red]  {e}[/red]")
        console_err.print(
            "[yellow]Hint: the file must contain the change-request fields as key: value pairs. "
            "Run `grzctl change-request-template <CHANGE>` for a valid example.[/yellow]"
        )
        raise click.Abort() from e
    if not isinstance(raw, dict):
        got = "an empty file" if raw is None else f"a top-level {type(raw).__name__}"
        console_err.print(
            f"[red]Error: '{data_file}' must contain a mapping of change-request fields at the top level, "
            f"but got {got}.[/red]"
        )
        console_err.print(
            "[yellow]Hint: run `grzctl change-request-template <CHANGE>` to see the expected structure. "
            "If you saved the template, make sure you replaced the placeholders rather than the field names.[/yellow]"
        )
        raise click.Abort()
    return raw


def _read_raw_content(path: Path | None) -> tuple[bytes | None, RequestRawContentType | None]:
    """Read the optional binary attachment and infer its type from the file extension."""
    if path is None:
        return None, None
    content_type = _RAW_CONTENT_TYPE_BY_SUFFIX.get(path.suffix.lower())
    if content_type is None:
        console_err.print(
            f"[red]Error: cannot infer raw-content type from extension '{path.suffix}'. "
            f"Supported extensions: {', '.join(sorted(_RAW_CONTENT_TYPE_BY_SUFFIX))}.[/red]"
        )
        raise click.Abort()
    return path.read_bytes(), content_type


def _build_change_request_kwargs(
    raw_input: dict | None,
    raw_content_bytes: bytes | None,
    raw_content_type: RequestRawContentType | None,
    change: ChangeRequestEnum,
) -> dict:
    """Split loaded input dict into columnar fields + nested ``data`` extras and validate.

    Required-ness for new entries is enforced here (the schema permits NULL so that
    historical rows pre-dating the audit columns can still be loaded).
    """
    raw_input = dict(raw_input or {})
    nested_extras = raw_input.pop("data", None)
    unknown_keys = {k: v for k, v in raw_input.items() if k not in _COLUMNAR_KEYS}
    extras = nested_extras if nested_extras else (unknown_keys or None)

    kwargs = {
        "requester_name": raw_input.get("requester_name"),
        "requester_email": raw_input.get("requester_email"),
        "requested_at": raw_input.get("requested_at"),
        "request_email_content": raw_input.get("request_email_content"),
        "request_raw_content": raw_content_bytes,
        "request_raw_content_type": raw_content_type,
        "data": extras,
    }

    # Collect every input problem in one pass so the user can fix them all at once,
    # rather than discovering them one error (and one re-run) at a time.
    problems = _collect_input_problems(kwargs)
    if problems:
        console_err.print("[red]Error: the change-request input needs attention:[/red]")
        for problem in problems:
            console_err.print(f"  [red]•[/red] {problem}")
        console_err.print(
            f"[yellow]See `grzctl change-request-template {change.value}` for a filled-in example.[/yellow]"
        )
        raise click.Abort()

    # Backstop: let the model catch anything the checks above don't cover.
    try:
        ChangeRequestLogBase(change=change, **kwargs)
    except ValidationError as e:
        console_err.print(f"[red]Error: data failed validation for change type '{change.value}':[/red]")
        click.echo(str(e), err=True)
        raise click.Abort() from e
    return kwargs


def _collect_input_problems(kwargs: dict) -> list[str]:
    """Return human-readable descriptions of every problem with the change-request input.

    A single pass over all fields — required-ness, leftover ``<FILL IN`` / ``YYYY-MM-DD``
    placeholders, date format, and the email-or-attachment rule — so callers can report
    everything at once instead of one error per run.
    """
    problems: list[str] = []

    for field in ("requester_name", "requester_email", "requested_at"):
        value = kwargs[field]
        if value is None or (isinstance(value, str) and not value.strip()):
            problems.append(f"{field}: required, but missing")
        elif isinstance(value, str) and _PLACEHOLDER_MARKER in value:
            problems.append(f"{field}: still has the '{_PLACEHOLDER_MARKER} ...>' placeholder — replace it")

    content = kwargs["request_email_content"]
    if isinstance(content, str) and _PLACEHOLDER_MARKER in content:
        problems.append(f"request_email_content: still has the '{_PLACEHOLDER_MARKER} ...>' placeholder — replace it")

    requested_at = kwargs["requested_at"]
    if isinstance(requested_at, str) and requested_at.strip() and _PLACEHOLDER_MARKER not in requested_at:
        if requested_at.strip().upper() == _DATE_PLACEHOLDER:
            problems.append(
                f"requested_at: still the '{_DATE_PLACEHOLDER}' placeholder — use a real date, e.g. 2026-03-27"
            )
        else:
            try:
                date.fromisoformat(requested_at)
            except ValueError:
                problems.append(
                    f"requested_at: invalid date format '{requested_at}' (expected YYYY-MM-DD, e.g. 2026-03-27)"
                )

    if kwargs["request_email_content"] is None and kwargs["request_raw_content"] is None:
        problems.append("provide either 'request_email_content' (in --data/--data-file) or --raw-content (or both)")

    return problems


def resolve_and_validate_change_request(
    change: ChangeRequestEnum,
    data_json: str | None,
    data_file: Path | None,
    raw_content_path: Path | None,
) -> dict:
    """Load and fully validate change-request input from CLI options.

    Performs no database access. Returns the kwargs for ``add_change_request``.
    On any problem, prints the error (and the expected template) and raises
    :class:`click.Abort`.
    """
    raw_content_bytes, raw_content_type = _read_raw_content(raw_content_path)
    raw_input = _load_change_request_input(data_json, data_file)
    if raw_input is None and raw_content_bytes is None:
        console_err.print(
            f"[red]Error: change request '{change.value}' requires audit data via --data or --data-file "
            f"(and/or --raw-content).[/red]"
        )
        console_err.print("[yellow]Expected fields (YAML template):[/yellow]")
        click.echo(render_yaml_template(change), err=True)
        raise click.Abort()
    return _build_change_request_kwargs(raw_input, raw_content_bytes, raw_content_type, change)


def _preview_kwargs(kwargs: dict) -> dict:
    """Copy of the validated kwargs with binary content elided for display."""
    preview = {k: v for k, v in kwargs.items() if k != "request_raw_content"}
    if kwargs["request_raw_content"] is not None:
        preview["request_raw_content"] = f"<{len(kwargs['request_raw_content'])} bytes>"
    return preview


@click.command("change-request-validate")
@click.argument("change_str", metavar="CHANGE", type=click.Choice(ChangeRequestEnum.list(), case_sensitive=False))
@click.option("--data", "data_json", type=str, default=None, help='Inline JSON data (e.g., \'{"k":"v"}\').')
@click.option(
    "--data-file",
    "data_file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Path to a JSON or YAML file with the change-request fields (see `grzctl change-request-template`).",
)
@click.option(
    "--raw-content",
    "raw_content_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Optional PDF/PNG attachment to include in validation; type is verified by magic bytes.",
)
def change_request_validate(
    change_str: str,
    data_json: str | None,
    data_file: Path | None,
    raw_content_path: Path | None,
):
    """Validate change-request input offline — no DB config, no writes.

    Runs exactly the same checks as ``db ... submission change-request`` (required
    audit fields, ``<FILL IN`` placeholder rejection, model validation) but never
    connects to a database and needs no ``--config-file``. It therefore cannot
    accidentally write a change request, making it the safe way to check a data
    file before the real command.
    """
    change = ChangeRequestEnum(change_str)
    kwargs = resolve_and_validate_change_request(change, data_json, data_file, raw_content_path)
    console_err.print(f"[green]OK: change-request input for '{change.value}' is valid.[/green]")
    console_err.print("[yellow]Validated fields:[/yellow]")
    click.echo(json.dumps(_preview_kwargs(kwargs), indent=2, ensure_ascii=False, default=str), err=True)
