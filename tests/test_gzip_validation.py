"""Tests for gzip magic bytes validation."""

import tempfile
from pathlib import Path

import pytest
from grz_common.validation import check_gzip_magic_bytes


def test_valid_gzip_file():
    """Test that a valid gzip file passes validation."""
    with tempfile.NamedTemporaryFile(suffix=".fastq.gz", delete=False) as tmp:
        # Write gzip magic bytes followed by some data
        tmp.write(b"\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\x03")
        tmp.flush()
        tmp_path = Path(tmp.name)

    try:
        error = check_gzip_magic_bytes(tmp_path)
        assert error is None
    finally:
        tmp_path.unlink()


def test_invalid_gzip_magic_bytes():
    """Test that a file with .gz extension but invalid magic bytes fails validation."""
    with tempfile.NamedTemporaryFile(suffix=".fastq.gz", delete=False) as tmp:
        tmp.write(b"This is not a gzip file")
        tmp.flush()
        tmp_path = Path(tmp.name)

    try:
        error = check_gzip_magic_bytes(tmp_path)
        assert error is not None
        assert "does not have gzip magic bytes" in error
        assert "0x54 0x68" in error  # "Th" in hex
    finally:
        tmp_path.unlink()


def test_empty_gz_file():
    """Test that an empty file with .gz extension fails validation."""
    with tempfile.NamedTemporaryFile(suffix=".fastq.gz", delete=False) as tmp:
        # Create an empty file
        tmp_path = Path(tmp.name)

    try:
        error = check_gzip_magic_bytes(tmp_path)
        assert error is not None
        assert "too small to contain gzip header" in error
    finally:
        tmp_path.unlink()


def test_single_byte_gz_file():
    """Test that a 1-byte file with .gz extension fails validation."""
    with tempfile.NamedTemporaryFile(suffix=".fastq.gz", delete=False) as tmp:
        tmp.write(b"\x1f")
        tmp.flush()
        tmp_path = Path(tmp.name)

    try:
        error = check_gzip_magic_bytes(tmp_path)
        assert error is not None
        assert "too small to contain gzip header" in error
    finally:
        tmp_path.unlink()


def test_non_gz_file_not_checked():
    """Test that files without .gz extension are not checked."""
    with tempfile.NamedTemporaryFile(suffix=".txt", delete=False) as tmp:
        tmp.write(b"This is a text file")
        tmp.flush()
        tmp_path = Path(tmp.name)

    try:
        error = check_gzip_magic_bytes(tmp_path)
        assert error is None  # Should not validate non-gz files
    finally:
        tmp_path.unlink()


def test_non_gz_fastq_file():
    """Test that uncompressed FASTQ files are not checked."""
    with tempfile.NamedTemporaryFile(suffix=".fastq", delete=False) as tmp:
        tmp.write(b"@SEQ1\nACGT\n+\nFFFF\n")
        tmp.flush()
        tmp_path = Path(tmp.name)

    try:
        error = check_gzip_magic_bytes(tmp_path)
        assert error is None  # Should not validate non-gz files
    finally:
        tmp_path.unlink()
