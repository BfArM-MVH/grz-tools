import json
import logging
import os
import shutil
from contextlib import contextmanager
from pathlib import Path

import grz_cli.cli
import pytest
from click.testing import CliRunner
from grz_common.workers.submission import SubmissionValidationError


@pytest.fixture
def chdir(tmp_path):
    """A fixture to temporarily change the working directory."""

    @contextmanager
    def _chdir(path):
        current_dir = Path.cwd()
        try:
            os.chdir(path)
            yield
        finally:
            os.chdir(current_dir)

    return _chdir


@pytest.mark.parametrize("grz_check_flag", ["--with-grz-check", "--no-grz-check"])
def test_validate_submission(
    temp_identifiers_config_file_path,
    working_dir_path,
    grz_check_flag,
    caplog,
):
    have_grz_check = shutil.which("grz-check") is not None

    if (grz_check_flag == "--with-grz-check") and not have_grz_check:
        # explicitly note when skipping this when grz-check not available instead of silently falling back
        pytest.skip(reason="grz-check not installed")

    submission_dir = Path("tests/mock_files/submissions/valid_submission")

    shutil.copytree(submission_dir / "files", working_dir_path / "files", dirs_exist_ok=True)
    shutil.copytree(submission_dir / "metadata", working_dir_path / "metadata", dirs_exist_ok=True)

    testargs = [
        "validate",
        "--config-file",
        temp_identifiers_config_file_path,
        "--submission-dir",
        str(working_dir_path),
        grz_check_flag,
    ]

    runner = CliRunner()
    cli = grz_cli.cli.build_cli()
    with caplog.at_level(logging.INFO):
        result = runner.invoke(cli, testargs, catch_exceptions=False)
        assert result.exit_code == 0, result.output

        if grz_check_flag == "--no-grz-check":
            assert "Starting checksum validation (fallback)..." in caplog.text
        else:
            assert "Starting file validation with `grz-check`..." in caplog.text
    caplog.clear()

    # check if re-validation is skipped
    with caplog.at_level(logging.INFO):
        result = runner.invoke(cli, testargs, catch_exceptions=False)

        if grz_check_flag == "--no-grz-check":
            assert "Starting checksum validation (fallback)..." in caplog.text
        else:
            assert "Starting file validation with `grz-check`..." in caplog.text

    # test if command has correctly checked for:
    # - mismatched md5sums
    # - all files existing

    assert result.exit_code == 0, result.output


def test_validate_submission_incorrect_grz_id(
    temp_identifiers_config_file_path,
    working_dir_path,
):
    submission_dir = Path("tests/mock_files/submissions/valid_submission")

    shutil.copytree(submission_dir / "files", working_dir_path / "files", dirs_exist_ok=True)
    shutil.copytree(submission_dir / "metadata", working_dir_path / "metadata", dirs_exist_ok=True)

    with open(working_dir_path / "metadata" / "metadata.json", mode="r+") as metadata_file:
        metadata = json.load(metadata_file)

        # put a GRZ id into the metadata that does not match what's in the config
        metadata["submission"]["genomicDataCenterId"] = "GRZX00000"

        metadata_file.seek(0)
        json.dump(metadata, metadata_file)
        metadata_file.truncate()

    testargs = [
        "validate",
        "--config-file",
        temp_identifiers_config_file_path,
        "--submission-dir",
        str(working_dir_path),
    ]

    runner = CliRunner()
    cli = grz_cli.cli.build_cli()
    # set catch_exceptions to True because we expect this to fail
    result = runner.invoke(cli, testargs, catch_exceptions=True)
    exc_type, exc, *_ = result.exc_info
    assert exc_type == SubmissionValidationError
    assert "does not match genomic data center identifier" in str(exc)

    assert result.exit_code == 1, result.output


@pytest.mark.parametrize("grz_check_flag", ["--with-grz-check", "--no-grz-check"])
def test_validate_submission_with_symlink(
    temp_identifiers_config_file_path, working_dir_path, grz_check_flag, caplog, chdir
):
    """
    Tests that the validation can handle relative symlinked files correctly.
    """
    have_grz_check = shutil.which("grz-check") is not None

    if (grz_check_flag == "--with-grz-check") and not have_grz_check:
        pytest.skip(reason="grz-check not installed")

    # setup dir for original files (with a name other than "files", just to make sure)
    source_data_dir = working_dir_path / "source_data"
    source_data_dir.mkdir()

    submission_dir = Path("tests/mock_files/submissions/valid_submission")
    mock_files_dir = submission_dir / "files"
    submission_files_dir = working_dir_path / "files"

    # copy mock files to _source_ directory
    shutil.copytree(mock_files_dir, source_data_dir, dirs_exist_ok=True)

    # create actual submission dir which should contain symlinks
    submission_files_dir.mkdir()

    # create symlinks (from submission 'files/' to 'source_data/')
    for source_file in source_data_dir.iterdir():
        link_path = submission_files_dir / source_file.name
        relative_target_path = os.path.relpath(source_file, start=submission_files_dir)
        link_path.symlink_to(relative_target_path)

    shutil.copytree(submission_dir / "metadata", working_dir_path / "metadata", dirs_exist_ok=True)

    config_file = Path(temp_identifiers_config_file_path).resolve()

    testargs = [
        "validate",
        "--config-file",
        str(config_file),
        "--submission-dir",
        ".",  # relative submission dir path
        grz_check_flag,
    ]

    runner = CliRunner()
    cli = grz_cli.cli.build_cli()

    # trigger validation from within working_dir_path
    with chdir(working_dir_path):
        with caplog.at_level(logging.INFO):
            result = runner.invoke(cli, testargs, catch_exceptions=False)

    assert result.exit_code == 0, result.output

    if grz_check_flag == "--no-grz-check":
        assert "Starting checksum validation (fallback)..." in caplog.text
    else:
        assert "Starting file validation with `grz-check`..." in caplog.text
        assert "Could not find metadata for file in grz-check report" not in caplog.text


@pytest.mark.parametrize("grz_check_flag", ["--with-grz-check", "--no-grz-check"])
def test_validate_submission_with_broken_symlink(
    grz_check_flag,
    temp_identifiers_config_file_path,
    working_dir_path,
):
    """
    Tests that the validation fails with a broken symlink.
    """
    have_grz_check = shutil.which("grz-check") is not None

    if (grz_check_flag == "--with-grz-check") and not have_grz_check:
        pytest.skip(reason="grz-check not installed")

    submission_dir = Path("tests/mock_files/submissions/valid_submission")
    submission_files_dir = working_dir_path / "files"

    submission_files_dir.mkdir()

    broken_link_path = submission_files_dir / "broken_link.bed"
    broken_link_path.symlink_to("non_existent.bed")

    shutil.copytree(submission_dir / "metadata", working_dir_path / "metadata", dirs_exist_ok=True)

    metadata_path = working_dir_path / "metadata" / "metadata.json"
    with open(metadata_path, "r+") as f:
        metadata = json.load(f)
        first = metadata["donors"][0]["labData"][0]["sequenceData"]["files"][0]
        first["filePath"] = "broken_link.bed"
        f.seek(0)
        json.dump(metadata, f)
        f.truncate()

    testargs = [
        "validate",
        "--config-file",
        temp_identifiers_config_file_path,
        "--submission-dir",
        str(working_dir_path),
        grz_check_flag,
    ]

    runner = CliRunner()
    cli = grz_cli.cli.build_cli()
    result = runner.invoke(cli, testargs, catch_exceptions=True)

    assert result.exit_code == 1
    assert isinstance(result.exception, SubmissionValidationError)
