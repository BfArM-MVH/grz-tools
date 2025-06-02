"""Tests for the db subcommand."""

import os

from click.testing import CliRunner

import grz_cli.cli


def test_db(
    temp_config_file_path,
):
    env = {"GRZ_DB_AUTHOR_PASSPHRASE": "test"}
    os.environ.update(env)

    runner = CliRunner(env={"GRZ_DB_AUTHOR_PASSPHRASE": "test"})
    cli = grz_cli.cli.build_cli(grz_mode=True)
    execute = lambda args: runner.invoke(cli, args, catch_exceptions=False)

    args_prefix = ["db", "--config-file", temp_config_file_path]
    # first initialize an empty DB
    init_args = [*args_prefix, "init"]
    result = execute(init_args)
    assert result.exit_code == 0, result.output

    # then add a submission
    add_args = [*args_prefix, "submission", "add", "S01", "--tan-g", "foo", "--pseudonym", "bar"]
    result = execute(add_args)
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
