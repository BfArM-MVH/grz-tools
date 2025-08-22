"""Tests for the db subcommand."""

import json

import grzctl.cli
from click.testing import CliRunner


def test_db(
    temp_db_config_file_path,
):
    env = {"GRZ_DB__AUTHOR__PRIVATE_KEY_PASSPHRASE": "test"}

    runner = CliRunner(env=env)
    cli = grzctl.cli.build_cli()
    execute = lambda args: runner.invoke(cli, args, catch_exceptions=False)

    args_prefix = ["db", "--config-file", temp_db_config_file_path]
    # first initialize an empty DB
    init_args = [*args_prefix, "init"]
    result = execute(init_args)
    assert result.exit_code == 0, result.output

    submission_id = "123456789_1970-01-01_a0b1c2d3"
    tan_g = "2c26b46b68ffc68ff99b453c1d30413413422d706483bfa0f98a5e886266e7ae"

    # then add a submission
    add_args = [*args_prefix, "submission", "add", submission_id]
    result = execute(add_args)
    assert result.exit_code == 0, result.output

    # then update the submission's tanG and/or pseudonym
    modify_args = [
        *args_prefix,
        "submission",
        "modify",
        submission_id,
        "tan_g",
        tan_g,
    ]
    result = execute(modify_args)
    assert result.exit_code == 0, result.output
    modify_args = [*args_prefix, "submission", "modify", submission_id, "pseudonym", "bar"]
    result = execute(modify_args)
    assert result.exit_code == 0, result.output

    # then update a submission
    # … downloading …
    update_args = [*args_prefix, "submission", "update", submission_id, "Downloading"]
    result = execute(update_args)
    assert result.exit_code == 0, result.output

    # … downloaded …
    update_args = [*args_prefix, "submission", "update", submission_id, "Downloaded"]
    result = execute(update_args)
    assert result.exit_code == 0, result.output

    # then show details for the submission
    show_args = [*args_prefix, "submission", "show", submission_id]
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
            "id": submission_id,
            "tan_g": tan_g,
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
