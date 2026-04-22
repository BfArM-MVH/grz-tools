"""Smoke and parity tests for the grz_check Python bindings."""

from __future__ import annotations

import gzip
import hashlib
import io
from pathlib import Path

import grz_check
import pytest

REPO_ROOT = Path(__file__).resolve().parents[3]
BAM = REPO_ROOT / "tests" / "resources" / "reads" / "valid_HiFi.bam"

FASTQ_BYTES = b"@read1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n@read2\nTTGGCCAATTGG\n+\nJJJJJJJJJJJJ\n"


@pytest.fixture
def fastq_file(tmp_path: Path) -> Path:
    path = tmp_path / "reads.fastq"
    path.write_bytes(FASTQ_BYTES)
    return path


def test_checksum_path_and_stream_agree(tmp_path: Path) -> None:
    data = b"hello grz-check\n"
    path = tmp_path / "x.bin"
    path.write_bytes(data)

    expected = hashlib.sha256(data).hexdigest()
    assert grz_check.calculate_checksum(path) == expected
    assert grz_check.calculate_checksum(io.BytesIO(data)) == expected


def test_validate_fastq_path(fastq_file: Path) -> None:
    report = grz_check.validate_fastq(fastq_file)
    assert report.is_valid
    assert report.num_records == 2
    assert report.mean_read_length == pytest.approx(12.0)


def test_validate_fastq_stream_matches_path(fastq_file: Path) -> None:
    path_report = grz_check.validate_fastq(fastq_file)
    stream_report = grz_check.validate_fastq(io.BytesIO(FASTQ_BYTES))
    assert path_report.is_valid == stream_report.is_valid
    assert path_report.num_records == stream_report.num_records
    assert path_report.mean_read_length == stream_report.mean_read_length


def test_validate_fastq_gzip_stream() -> None:
    gz_bytes = gzip.compress(FASTQ_BYTES)
    with gzip.GzipFile(fileobj=io.BytesIO(gz_bytes), mode="rb") as reader:
        report = grz_check.validate_fastq(reader)
    assert report.is_valid
    assert report.num_records == 2


def test_validate_fastq_gzip_raw_bytes() -> None:
    """Niffler on the Rust side decompresses gzip transparently — no gzip.open() needed."""
    gz_bytes = gzip.compress(FASTQ_BYTES)
    report = grz_check.validate_fastq(io.BytesIO(gz_bytes))
    assert report.is_valid
    assert report.num_records == 2


def test_validate_fastq_gzip_filelike_raw() -> None:
    """Passing a raw binary file handle to a .gz file — niffler decompresses transparently."""
    gz_bytes = gzip.compress(FASTQ_BYTES)
    report = grz_check.validate_fastq(io.BytesIO(gz_bytes))
    assert report.is_valid
    assert report.num_records == 2


def test_validate_fastq_buffer_protocol_bytes() -> None:
    """bytes objects go through the PyBuffer fast path (no .read() loop)."""
    report = grz_check.validate_fastq(FASTQ_BYTES)
    assert report.is_valid
    assert report.num_records == 2


def test_validate_fastq_buffer_protocol_bytearray() -> None:
    """bytearray goes through PyBuffer (buffer protocol)."""
    report = grz_check.validate_fastq(bytearray(FASTQ_BYTES))
    assert report.is_valid
    assert report.num_records == 2


def test_validate_fastq_buffer_protocol_memoryview() -> None:
    """memoryview goes through PyBuffer."""
    report = grz_check.validate_fastq(memoryview(FASTQ_BYTES))
    assert report.is_valid
    assert report.num_records == 2


def test_validate_fastq_buffer_protocol_gzipped_bytes() -> None:
    """Gzipped bytes via buffer protocol: niffler auto-detects and decompresses."""
    gz_bytes = gzip.compress(FASTQ_BYTES)
    report = grz_check.validate_fastq(gz_bytes)
    assert report.is_valid
    assert report.num_records == 2


def test_checksum_buffer_protocol_bytes() -> None:
    """Checksum over raw bytes via buffer protocol."""
    data = b"hello grz-check\n"
    expected = hashlib.sha256(data).hexdigest()
    assert grz_check.calculate_checksum(data) == expected


def test_validate_fastq_min_length_pass(fastq_file: Path) -> None:
    report = grz_check.validate_fastq(fastq_file, min_mean_read_length=5)
    assert report.is_valid


def test_validate_fastq_min_length_fail(fastq_file: Path) -> None:
    report = grz_check.validate_fastq(fastq_file, min_mean_read_length=100)
    assert not report.is_valid
    assert report.errors


def test_validate_fastq_paired_path_and_stream_agree(tmp_path: Path) -> None:
    r1 = tmp_path / "r1.fastq"
    r2 = tmp_path / "r2.fastq"
    r1.write_bytes(FASTQ_BYTES)
    r2.write_bytes(FASTQ_BYTES)

    p1, p2 = grz_check.validate_fastq_paired(r1, r2)
    s1, s2 = grz_check.validate_fastq_paired(io.BytesIO(FASTQ_BYTES), io.BytesIO(FASTQ_BYTES))
    assert p1.num_records == s1.num_records == 2
    assert p2.num_records == s2.num_records == 2
    assert p1.is_valid == s1.is_valid
    assert p2.is_valid == s2.is_valid


@pytest.mark.skipif(not BAM.exists(), reason="BAM fixture not materialized (git lfs pull)")
def test_validate_bam_path_valid() -> None:
    report = grz_check.validate_bam(BAM)
    assert report.is_valid
    assert report.num_records is not None and report.num_records > 0


@pytest.mark.skipif(not BAM.exists(), reason="BAM fixture not materialized (git lfs pull)")
def test_validate_bam_stream_matches_path() -> None:
    path_report = grz_check.validate_bam(BAM)
    with open(BAM, "rb") as handle:
        stream_report = grz_check.validate_bam(handle)
    assert path_report.is_valid == stream_report.is_valid
    assert path_report.num_records == stream_report.num_records


def test_text_mode_file_rejected(fastq_file: Path) -> None:
    with open(fastq_file) as text_handle:
        with pytest.raises(ValueError, match="binary mode"):
            grz_check.validate_fastq(text_handle)


def test_non_filelike_rejected() -> None:
    with pytest.raises(TypeError):
        grz_check.validate_fastq(42)
