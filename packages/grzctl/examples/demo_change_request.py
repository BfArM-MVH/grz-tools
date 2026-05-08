#!/usr/bin/env python
"""End-to-end demo of the schema-validated DELETE change-request flow.

Spins up a throwaway sqlite database in a temp dir, adds a submission, then
exercises the template, validation, dry-run, and real-write paths via Click's
``CliRunner`` — i.e. the same code paths the unit tests cover, but with the
output printed to your terminal so you can watch the flow.

Run::

    uv run python packages/grzctl/examples/demo_change_request.py
    # add --keep to preserve the temp dir for inspection
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
import tempfile
from pathlib import Path

import click.testing
import cryptography.hazmat.primitives.serialization as cryptser
import grzctl.cli
import yaml
from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PrivateKey

SUB_ID = "260840108_2026-03-27_deadbeef"


def _setup_workspace(workspace: Path) -> Path:
    """Generate keys + config, return the path to the YAML config file."""
    private_key = Ed25519PrivateKey.generate()
    private_key_path = workspace / "alice.sec"
    private_key_path.write_bytes(
        private_key.private_bytes(
            encoding=cryptser.Encoding.PEM,
            format=cryptser.PrivateFormat.OpenSSH,
            encryption_algorithm=cryptser.NoEncryption(),
        )
    )

    public_key_path = workspace / "alice.pub"
    public_key_bytes = private_key.public_key().public_bytes(
        encoding=cryptser.Encoding.OpenSSH,
        format=cryptser.PublicFormat.OpenSSH,
    )
    public_key_path.write_bytes(public_key_bytes + b" alice")

    config_path = workspace / "config.yaml"
    config_path.write_text(
        yaml.dump(
            {
                "db": {
                    "database_url": f"sqlite:///{workspace / 'demo.sqlite'}",
                    "author": {
                        "name": "alice",
                        "private_key_path": str(private_key_path),
                        "private_key_passphrase": "",
                    },
                    "known_public_keys": str(public_key_path),
                }
            }
        )
    )
    return config_path


def _step(title: str) -> None:
    print()
    print(f"=== {title} ===")


def _run(runner: click.testing.CliRunner, cli, args: list[str]) -> click.testing.Result:
    result = runner.invoke(cli, args)
    sys.stdout.write(result.stdout)
    if result.stderr:
        sys.stdout.write(result.stderr)
    sys.stdout.flush()
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument(
        "--keep",
        action="store_true",
        help="Do not delete the temporary workspace at exit (path will be printed).",
    )
    args = parser.parse_args()

    workspace = Path(tempfile.mkdtemp(prefix="grzctl-demo-"))
    try:
        cfg = _setup_workspace(workspace)
        cli = grzctl.cli.build_cli()
        runner = click.testing.CliRunner()

        common = ["db", "--config-file", str(cfg)]
        cr_cmd = [*common, "submission", "change-request", SUB_ID, "Delete"]

        _step(f"0. Init throwaway DB at {workspace} and add submission")
        _run(runner, cli, [*common, "init"])
        _run(runner, cli, [*common, "submission", "add", SUB_ID])

        _step("1. Print the YAML template (no config required)")
        _run(runner, cli, ["change-request-template", "Delete"])

        _step("2. Submitting the unfilled template fails (date placeholder)")
        template = runner.invoke(cli, ["change-request-template", "Delete"]).stdout
        as_is = workspace / "as-is.yaml"
        as_is.write_text(template)
        _run(runner, cli, [*cr_cmd, "--data-file", str(as_is), "--dry-run"])

        _step("3. Filling only the date still fails (placeholder validator catches the rest)")
        partial = workspace / "partial.yaml"
        partial.write_text(template.replace("YYYY-MM-DD", "2026-03-27"))
        _run(runner, cli, [*cr_cmd, "--data-file", str(partial), "--dry-run"])

        _step("4. Author a fully filled-in YAML data file")
        good_data = {
            "requester_name": "Erika Mustermann",
            "requester_email": "demo.requester@example.org",
            "requested_at": "2026-03-27",
            "requester_email_content": (
                "Liebe Kolleginnen und Kollegen vom Beispielklinikum,\n"
                f"bitte den Datensatz mit der Submission-ID {SUB_ID} löschen.\n"
                "Vielen Dank, Erika\n"
            ),
            "internal_note": "demo placeholder note",
        }
        good = workspace / "delete.yaml"
        good.write_text(yaml.safe_dump(good_data, allow_unicode=True))
        print(good.read_text())

        _step("5. Dry-run with valid data — validates and reports, no DB write")
        _run(runner, cli, [*cr_cmd, "--data-file", str(good), "--dry-run"])

        _step("6. Real insert")
        _run(runner, cli, [*cr_cmd, "--data-file", str(good)])

        _step("7. Confirm what landed in the DB")
        result = runner.invoke(cli, [*common, "list-change-requests", "--json"])
        sys.stdout.write(json.dumps(json.loads(result.stdout), indent=2, ensure_ascii=False))
        sys.stdout.write("\n")

        _step("8. Same content via --data (inline JSON) — same constraints apply")
        _run(runner, cli, [*cr_cmd, "--data", json.dumps(good_data), "--dry-run"])

        _step("9. Modify is unschematized — loose pass-through still works")
        _run(runner, cli, [*common, "submission", "change-request", SUB_ID, "Modify"])

        print()
        print("Demo complete.")
        return 0
    finally:
        if args.keep:
            print(f"\n[--keep] workspace preserved at: {workspace}")
        else:
            shutil.rmtree(workspace, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
