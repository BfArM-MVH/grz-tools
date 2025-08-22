"""Tests for grzctl report."""

import datetime

import grzctl.cli
from click.testing import CliRunner


def test_report_processed(temp_db_config_file_path):
    env = {"GRZ_DB__AUTHOR__PRIVATE_KEY_PASSPHRASE": "test"}

    runner = CliRunner(env=env)
    cli = grzctl.cli.build_cli()
    execute = lambda args: runner.invoke(cli, args, catch_exceptions=False)

    args_prefix = ["db", "--config-file", temp_db_config_file_path]
    result = execute([*args_prefix, "init"])
    assert result.exit_code == 0, result.output

    submission_id = "123456789_1970-01-01_a0b1c2d3"
    result = execute([*args_prefix, "submission", "add", submission_id])
    assert result.exit_code == 0, result.output

    result = execute(
        [
            *args_prefix,
            "submission",
            "modify",
            submission_id,
            "basic_qc_passed",
            "True",
        ]
    )
    assert result.exit_code == 0, result.output
    result = execute(
        [
            *args_prefix,
            "submission",
            "modify",
            submission_id,
            "detailed_qc_passed",
            "False",
        ]
    )
    assert result.exit_code == 0, result.output

    today = datetime.date.today()
    report_args = [
        "report",
        "--config-file",
        temp_db_config_file_path,
        "processed",
        "--since",
        (today - datetime.timedelta(days=1)).isoformat(),
        "--until",
        (today + datetime.timedelta(days=1)).isoformat(),
    ]

    # No submissions with state "reported" or "qced" in the db
    result = runner.invoke(cli, report_args, catch_exceptions=False)
    _comment, _header, *lines = result.stdout.splitlines()
    assert not lines

    # Add a submission with state "reported"
    result = execute(
        [
            *args_prefix,
            "submission",
            "update",
            submission_id,
            "reported",
        ]
    )
    assert result.exit_code == 0, result.output

    result = runner.invoke(cli, report_args, catch_exceptions=False)
    _comment, _header, *lines = result.stdout.splitlines()
    assert f"{submission_id}\tyes\tno\t{today.isoformat()}" == lines[0], lines
