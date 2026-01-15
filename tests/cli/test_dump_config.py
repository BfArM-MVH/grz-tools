import json

import yaml
from click.testing import CliRunner
from grz_cli.cli import build_cli
from grz_common.utils.config import merge_config_dicts


def test_dump_config_no_files(caplog):
    """
    GIVEN the grz-cli
    WHEN the dump-config command is called without any config files
    THEN the command should succeed and log that no config files were loaded
    """
    runner = CliRunner()
    cli = build_cli()
    result = runner.invoke(cli, ["dump-config"])

    assert result.exit_code == 0
    if "config.yaml" in caplog.text:
        assert "Configuration files to load" in caplog.text
    else:
        assert "Configuration files to load: []" in caplog.text
        assert "Merged configuration: {}" in caplog.text


def test_dump_config_with_files(
    caplog,
    temp_s3_config_file_path,
    temp_db_config_file_path,
):
    """
    GIVEN the grz-cli
    WHEN the dump-config command is called with multiple config files
    THEN the command should succeed and log the merged configuration
    """
    runner = CliRunner()
    cli = build_cli()
    result = runner.invoke(
        cli,
        ["dump-config", "--config-file", str(temp_s3_config_file_path), "--config-file", str(temp_db_config_file_path)],
    )

    assert result.exit_code == 0

    # Check that the loaded files are logged
    assert str(temp_s3_config_file_path) in caplog.text
    assert str(temp_db_config_file_path) in caplog.text

    # Check that the merged configuration is logged and correct
    merged_config_log = [rec for rec in caplog.records if "Merged configuration" in rec.message]
    assert len(merged_config_log) == 1
    logged_config_str = merged_config_log[0].message.replace("Merged configuration: ", "")
    logged_config = json.loads(logged_config_str)

    with open(temp_s3_config_file_path) as fd:
        s3_config_content = yaml.safe_load(fd)
    with open(temp_db_config_file_path) as fd:
        db_config_content = yaml.safe_load(fd)

    expected_config = merge_config_dicts(s3_config_content, db_config_content)

    assert logged_config == expected_config

def test_dump_config_global_args(
    caplog,
    temp_s3_config_file_path,
    temp_db_config_file_path,
):
    """
    GIVEN the grz-cli
    WHEN the dump-config command is called with multiple config files
    THEN the command should succeed and log the merged configuration
    """
    runner = CliRunner()
    cli = build_cli()
    result = runner.invoke(
        cli,
        ["--config-file", str(temp_s3_config_file_path), "dump-config", "--config-file", str(temp_db_config_file_path)],
    )

    assert result.exit_code == 0

    # Check that the loaded files are logged
    assert str(temp_s3_config_file_path) in caplog.text
    assert str(temp_db_config_file_path) in caplog.text

    # Check that the merged configuration is logged and correct
    merged_config_log = [rec for rec in caplog.records if "Merged configuration" in rec.message]
    assert len(merged_config_log) == 1
    logged_config_str = merged_config_log[0].message.replace("Merged configuration: ", "")
    logged_config = json.loads(logged_config_str)

    with open(temp_s3_config_file_path) as fd:
        s3_config_content = yaml.safe_load(fd)
    with open(temp_db_config_file_path) as fd:
        db_config_content = yaml.safe_load(fd)

    expected_config = merge_config_dicts(s3_config_content, db_config_content)

    assert logged_config == expected_config
