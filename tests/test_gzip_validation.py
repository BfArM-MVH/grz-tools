"""Tests for gzip magic bytes validation."""

from pathlib import Path

import pytest
from grz_common.validation import check_gzip_magic_bytes


def test_valid_gzip_file(tmp_path):
    """Test that a valid gzip file passes validation."""
    gz_file = tmp_path / "test.fastq.gz"
    # Write gzip magic bytes followed by some data
    gz_file.write_bytes(b"\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\x03")
    
    error = check_gzip_magic_bytes(gz_file)
    assert error is None


def test_invalid_gzip_magic_bytes(tmp_path):
    """Test that a file with .gz extension but invalid magic bytes fails validation."""
    gz_file = tmp_path / "fake.fastq.gz"
    gz_file.write_bytes(b"This is not a gzip file")
    
    error = check_gzip_magic_bytes(gz_file)
    assert error is not None
    assert "does not have gzip magic bytes" in error
    assert "0x54 0x68" in error  # "Th" in hex


def test_empty_gz_file(tmp_path):
    """Test that an empty file with .gz extension fails validation."""
    gz_file = tmp_path / "empty.fastq.gz"
    gz_file.write_bytes(b"")
    
    error = check_gzip_magic_bytes(gz_file)
    assert error is not None
    assert "too small to contain gzip header" in error


def test_single_byte_gz_file(tmp_path):
    """Test that a 1-byte file with .gz extension fails validation."""
    gz_file = tmp_path / "single.fastq.gz"
    gz_file.write_bytes(b"\x1f")
    
    error = check_gzip_magic_bytes(gz_file)
    assert error is not None
    assert "too small to contain gzip header" in error


def test_non_gz_file_not_checked(tmp_path):
    """Test that files without .gz extension are not checked."""
    txt_file = tmp_path / "test.txt"
    txt_file.write_bytes(b"This is a text file")
    
    error = check_gzip_magic_bytes(txt_file)
    assert error is None  # Should not validate non-gz files


def test_non_gz_fastq_file(tmp_path):
    """Test that uncompressed FASTQ files are not checked."""
    fastq_file = tmp_path / "test.fastq"
    fastq_file.write_bytes(b"@SEQ1\nACGT\n+\nFFFF\n")
    
    error = check_gzip_magic_bytes(fastq_file)
    assert error is None  # Should not validate non-gz files
