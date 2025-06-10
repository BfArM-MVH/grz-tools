import shutil
from pathlib import Path

import grz_cli.cli
from click.testing import CliRunner


def test_validate_submission(
    working_dir_path,
):
    submission_dir = Path("tests/mock_files/submissions/valid_submission")

    shutil.copytree(submission_dir / "files", working_dir_path / "files", dirs_exist_ok=True)
    shutil.copytree(submission_dir / "metadata", working_dir_path / "metadata", dirs_exist_ok=True)

    testargs = [
        "validate",
        "--submission-dir",
        str(working_dir_path),
    ]

    runner = CliRunner()
    cli = grz_cli.cli.build_cli()
    result = runner.invoke(cli, testargs, catch_exceptions=False)

    # check if re-validation is skipped
    result = runner.invoke(cli, testargs, catch_exceptions=False)

    # test if command has correctly checked for:
    # - mismatched md5sums
    # - all files existing

    assert result.exit_code == 0, result.output
