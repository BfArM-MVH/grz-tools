"""grzctl sets ``local_storage_mode`` on the local directories that hold submission data, whatever the umask."""

import os
import stat
from pathlib import Path

import grzctl.cli
import pytest
import yaml
from click.testing import CliRunner
from grz_common.workers.submission import EncryptedSubmission

from ..conftest import _GRZ_PRIVATE_KEY_PATH, _grzctl_config_dict
from .common import copy_submission

LOCAL_STORAGE_MODE = 0o750
"""The configured mode. Under :func:`strict_umask`, a new directory gets it only from ``chmod``."""

NESTED_FILE = "regions/target_regions.bed"
"""A file path of the submission with a subdirectory."""


def _mode(path: Path) -> int:
    return stat.S_IMODE(path.stat().st_mode)


@pytest.fixture
def strict_umask():
    """A umask that clears the group and other bits, so that ``mkdir`` alone gives no :data:`LOCAL_STORAGE_MODE`."""
    previous_umask = os.umask(0o077)
    yield
    os.umask(previous_umask)


@pytest.fixture
def config_path(tmp_path: Path) -> Path:
    """A grzctl config whose inbox ``testing`` decrypts the example submission, with :data:`LOCAL_STORAGE_MODE`."""
    data = _grzctl_config_dict(
        leistungserbringer={"260914050": {"inbox_buckets": {"testing": {"private_key_path": _GRZ_PRIVATE_KEY_PATH}}}},
    )
    data["local_storage_mode"] = f"{LOCAL_STORAGE_MODE:04o}"
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(data))
    return config_path


@pytest.fixture
def nested_submission_dir(tmp_path: Path) -> Path:
    """The encrypted example submission, with ``target_regions.bed`` moved to :data:`NESTED_FILE`."""
    submission_dir = tmp_path / "nested_submission"
    copy_submission(submission_dir, "encrypted_files", "metadata")
    metadata_path = submission_dir / "metadata" / "metadata.json"
    metadata_path.write_text(metadata_path.read_text().replace('"target_regions.bed"', f'"{NESTED_FILE}"'))
    encrypted_files_dir = submission_dir / "encrypted_files"
    (encrypted_files_dir / NESTED_FILE).parent.mkdir()
    (encrypted_files_dir / "target_regions.bed.c4gh").rename(encrypted_files_dir / f"{NESTED_FILE}.c4gh")
    return submission_dir


@pytest.mark.parametrize("existing_metadata_dir", [False, True])
def test_download_sets_the_local_storage_mode(
    config_path: Path,
    nested_submission_dir: Path,
    remote_bucket_with_version,
    tmp_path: Path,
    strict_umask,
    existing_metadata_dir: bool,
):
    encrypted_submission = EncryptedSubmission(
        nested_submission_dir / "metadata", nested_submission_dir / "encrypted_files"
    )
    metadata_path, metadata_key = encrypted_submission.get_metadata_file_path_and_object_id()
    remote_bucket_with_version.upload_file(metadata_path, metadata_key)
    for local_file_path, s3_key in encrypted_submission.get_encrypted_files_and_object_id().items():
        remote_bucket_with_version.upload_file(local_file_path, s3_key)
    output_dir = tmp_path / "download"
    if existing_metadata_dir:
        (output_dir / "metadata").mkdir(parents=True)
        output_dir.chmod(0o755)
        (output_dir / "metadata").chmod(0o755)
    parent_mode = _mode(tmp_path)

    result = CliRunner().invoke(
        grzctl.cli.build_cli(),
        [
            "--config",
            str(config_path),
            "download",
            "--submission-id",
            encrypted_submission.submission_id,
            "--output-dir",
            str(output_dir),
            "--no-update-db",
            "--no-populate",
            "--inbox",
            "testing",
        ],
    )

    assert result.exit_code == 0, result.output
    assert (output_dir / "encrypted_files" / f"{NESTED_FILE}.c4gh").is_file()
    for directory in (
        output_dir,
        output_dir / "metadata",
        output_dir / "encrypted_files",
        output_dir / "encrypted_files" / "regions",
        output_dir / "logs",
    ):
        assert _mode(directory) == LOCAL_STORAGE_MODE, directory
    assert _mode(tmp_path) == parent_mode


def test_decrypt_sets_the_local_storage_mode(config_path: Path, nested_submission_dir: Path, strict_umask):
    files_dir = nested_submission_dir / "files"
    files_dir.mkdir()
    files_dir.chmod(0o755)
    submission_dir_mode = _mode(nested_submission_dir)

    result = CliRunner().invoke(
        grzctl.cli.build_cli(),
        [
            "--config",
            str(config_path),
            "decrypt",
            "--submission-dir",
            str(nested_submission_dir),
            "--inbox",
            "testing",
            "--no-update-db",
        ],
    )

    assert result.exit_code == 0, result.output
    assert (files_dir / NESTED_FILE).is_file()
    for directory in (files_dir, files_dir / "regions", nested_submission_dir / "logs"):
        assert _mode(directory) == LOCAL_STORAGE_MODE, directory
    assert _mode(nested_submission_dir) == submission_dir_mode
