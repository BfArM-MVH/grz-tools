import grzctl.cli
from click.testing import CliRunner

from .common import copy_submission


def test_consent_submission(working_dir_path):
    copy_submission(working_dir_path, "files", "metadata")

    testargs = [
        "consent",
        "--date",
        "2025-06-25",
        "--submission-dir",
        str(working_dir_path),
    ]
    runner = CliRunner()
    cli = grzctl.cli.build_cli()
    result = runner.invoke(cli, testargs, catch_exceptions=False)
    assert result.stdout.strip() == "true"

    testargs += ["--json", "--details"]
    result = runner.invoke(cli, testargs, catch_exceptions=False)
    assert result.stdout.strip() == (
        '{"index": true, "bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111": true}'
    )
