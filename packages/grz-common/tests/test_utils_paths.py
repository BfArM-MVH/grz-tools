"""``ensure_directory_mode`` sets the mode of the directories below a base directory, and of nothing else."""

import os
import stat
from pathlib import Path

import pytest
from grz_common.utils.paths import ensure_directory_mode


def _mode(path: Path) -> int:
    return stat.S_IMODE(path.stat().st_mode)


def test_ensure_directory_mode_sets_the_mode_below_base_only(tmp_path: Path):
    base = tmp_path / "base"
    base.mkdir()
    base.chmod(0o755)
    existing = base / "existing"
    existing.mkdir()
    existing.chmod(0o755)

    ensure_directory_mode(base, existing / "new" / "leaf", 0o750)

    assert _mode(tmp_path / "base") == 0o755
    assert _mode(existing) == 0o750
    assert _mode(existing / "new") == 0o750
    assert _mode(existing / "new" / "leaf") == 0o750


def test_ensure_directory_mode_ignores_the_umask(tmp_path: Path):
    """The umask would clear the group bits of a new directory, but ``chmod`` sets them again."""
    previous_umask = os.umask(0o077)
    try:
        ensure_directory_mode(tmp_path, tmp_path / "submission", 0o770)
    finally:
        os.umask(previous_umask)

    assert _mode(tmp_path / "submission") == 0o770


def test_ensure_directory_mode_without_a_mode_leaves_it_to_the_umask(tmp_path: Path):
    existing = tmp_path / "existing"
    existing.mkdir()
    existing.chmod(0o755)

    previous_umask = os.umask(0o022)
    try:
        ensure_directory_mode(tmp_path, existing / "new", None)
    finally:
        os.umask(previous_umask)

    assert _mode(existing) == 0o755
    assert _mode(existing / "new") == 0o750


def test_ensure_directory_mode_refuses_a_directory_outside_base(tmp_path: Path):
    base = tmp_path / "base"
    base.mkdir()

    with pytest.raises(ValueError):
        ensure_directory_mode(base, tmp_path / "elsewhere", 0o770)

    assert not (tmp_path / "elsewhere").exists()
