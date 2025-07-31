"""Tests for the db subcommand."""

import json
import os

import grzctl.cli
from click.testing import CliRunner


def test_db(
    temp_db_config_file_path,
):
    env = {"GRZ_DB__AUTHOR__PRIVATE_KEY_PASSPHRASE": "test"}
    os.environ.update(env)

    runner = CliRunner(env={"GRZ_DB__AUTHOR__PRIVATE_KEY_PASSPHRASE": "test"})
    cli = grzctl.cli.build_cli()
    execute = lambda args: runner.invoke(cli, args, catch_exceptions=False)

    args_prefix = ["db", "--config-file", temp_db_config_file_path]
    # first initialize an empty DB
    init_args = [*args_prefix, "init"]
    result = execute(init_args)
    assert result.exit_code == 0, result.output

    # then add a submission
    add_args = [*args_prefix, "submission", "add", "S01"]
    result = execute(add_args)
    assert result.exit_code == 0, result.output

    # then update the submission's tanG and/or pseudonym
    modify_args = [
        *args_prefix,
        "submission",
        "modify",
        "S01",
        "tan_g",
        "2c26b46b68ffc68ff99b453c1d30413413422d706483bfa0f98a5e886266e7ae",
    ]
    result = execute(modify_args)
    assert result.exit_code == 0, result.output
    modify_args = [*args_prefix, "submission", "modify", "S01", "pseudonym", "bar"]
    result = execute(modify_args)
    assert result.exit_code == 0, result.output

    # then update a submission
    # … downloading …
    update_args = [*args_prefix, "submission", "update", "S01", "Downloading"]
    result = execute(update_args)
    assert result.exit_code == 0, result.output

    # … downloaded …
    update_args = [*args_prefix, "submission", "update", "S01", "Downloaded"]
    result = execute(update_args)
    assert result.exit_code == 0, result.output

    # then show details for the submission
    show_args = [*args_prefix, "submission", "show", "S01"]
    result = execute(show_args)
    assert result.exit_code == 0, result.output

    # then list all submissions
    list_args = [*args_prefix, "list"]
    result = execute(list_args)
    assert result.exit_code == 0, result.output

    # do the same again with json output
    list_args = [*args_prefix, "list", "--json"]
    result = execute(list_args)
    assert result.exit_code == 0, result.output
    expected_output = [
        {
            "id": "S01",
            "tan_g": "2c26b46b68ffc68ff99b453c1d30413413422d706483bfa0f98a5e886266e7ae",
            "pseudonym": "bar",
            "latest_state": {
                "state": "Downloaded",
                "timestamp": "MASKED",
                "data": None,
                "data_steward": "Alice",
                "data_steward_signature": "Verified",
            },
            "latest_change_request": {},
        }
    ]
    actual_output = json.loads(result.output)
    actual_output[0]["latest_state"]["timestamp"] = "MASKED"
    assert actual_output == expected_output
