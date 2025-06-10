import shutil
from pathlib import Path

import grzctl.cli
from click.testing import CliRunner


def test_consent_submission(working_dir_path):
    submission_dir = Path("tests/mock_files/submissions/valid_submission")

    shutil.copytree(submission_dir / "files", working_dir_path / "files", dirs_exist_ok=True)
    shutil.copytree(submission_dir / "metadata", working_dir_path / "metadata", dirs_exist_ok=True)

    testargs = [
        "consent",
        "--submission-dir",
        str(submission_dir),
    ]
    runner = CliRunner()
    cli = grzctl.cli.build_cli()
    result = runner.invoke(cli, testargs, catch_exceptions=False)
    assert result.stdout.strip() == "true"

    testargs += ["--json", "--details"]
    result = runner.invoke(cli, testargs, catch_exceptions=False)
    assert result.stdout.strip() == (
        "{"
        '"aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000": true, '
        '"bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111": true'
        "}"
    )
