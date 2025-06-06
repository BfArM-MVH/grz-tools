import filecmp
import os
import shutil
from pathlib import Path
from unittest import mock

import grz_cli.cli
import grzctl.cli
from click.testing import CliRunner


def are_dir_trees_equal(dir1, dir2):
    """
    Compare two directories recursively. Files in each directory are
    assumed to be equal if their names and contents are equal.

    :param dir1: First directory path
    :param dir2: Second directory path
    :return: True if the directory trees are the same and
        there were no errors while accessing the directories or files,
        False otherwise.
    """
    dirs_cmp = filecmp.dircmp(dir1, dir2)
    if len(dirs_cmp.left_only) > 0 or len(dirs_cmp.right_only) > 0 or len(dirs_cmp.funny_files) > 0:
        return False
    (_, mismatch, errors) = filecmp.cmpfiles(dir1, dir2, dirs_cmp.common_files, shallow=False)
    if len(mismatch) > 0 or len(errors) > 0:
        return False
    for common_dir in dirs_cmp.common_dirs:
        new_dir1 = os.path.join(dir1, common_dir)
        new_dir2 = os.path.join(dir2, common_dir)
        if not are_dir_trees_equal(new_dir1, new_dir2):
            return False
    return True


def test_upload_download_submission(
    working_dir_path,
    tmpdir_factory,
    temp_s3_config_file_path,
    remote_bucket,
):
    submission_dir = Path("tests/mock_files/submissions/valid_submission")

    shutil.copytree(
        submission_dir / "encrypted_files",
        working_dir_path / "encrypted_files",
        dirs_exist_ok=True,
    )
    shutil.copytree(submission_dir / "metadata", working_dir_path / "metadata", dirs_exist_ok=True)

    with mock.patch(
        "grz_common.models.s3.S3Options.__getattr__",
        lambda self, name: None if name == "endpoint_url" else AttributeError,
    ):
        # upload encrypted submission
        upload_args = [
            "upload",
            "--submission-dir",
            str(working_dir_path),
            "--config-file",
            temp_s3_config_file_path,
        ]

        runner = CliRunner()
        cli = grz_cli.cli.build_cli()
        result = runner.invoke(cli, upload_args, catch_exceptions=False)

        assert result.exit_code == 0, result.output
        assert len(result.output) != 0, result.stderr

        submission_id = result.stdout.strip()

        objects_in_bucket = list(remote_bucket.objects.all())
        assert len(objects_in_bucket) > 0, "Upload failed: No objects were found in the mock S3 bucket!"

        # download
        download_dir = tmpdir_factory.mktemp("submission_download")
        download_dir_path = Path(download_dir.strpath)

        # download encrypted submission
        download_args = [
            "download",
            "--submission-id",
            submission_id,
            "--output-dir",
            str(download_dir_path),
            "--config-file",
            temp_s3_config_file_path,
        ]
        cli = grzctl.cli.build_cli()
        result = runner.invoke(cli, download_args, catch_exceptions=False)

        assert result.exit_code == 0, result.output

    assert are_dir_trees_equal(
        working_dir_path / "encrypted_files",
        download_dir_path / "encrypted_files",
    ), "Encrypted files are different!"
    assert are_dir_trees_equal(
        working_dir_path / "metadata",
        download_dir_path / "metadata",
    ), "Metadata is different!"
