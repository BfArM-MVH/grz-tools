"""Tests for the bam_validation module."""

import importlib.resources
import logging

import grz_check
from grz_common.pipeline.components import ReadStream
from grz_common.pipeline.components.validation import BamValidator

from . import resources


def test_valid_hifi_bam():
    """Valid HiFi BAM files should pass validation"""
    bam_ptr = importlib.resources.files(resources).joinpath("reads", "valid_HiFi.bam")
    with importlib.resources.as_file(bam_ptr) as bam_path:
        report = grz_check.validate_bam(bam_path)

    assert report.is_valid
    assert len(report.errors) == 0


def test_valid_hifi_bam_stream():
    """Valid HiFi BAM files should pass validation"""
    bam_ptr = importlib.resources.files(resources).joinpath("reads", "valid_HiFi.bam")

    with importlib.resources.as_file(bam_ptr) as bam_path:
        validator = BamValidator()
        with open(bam_path, "rb") as f:
            ReadStream(f) >> validator


def test_hard_clipped_primary():
    """A warning should be present if hard-clipped bases are detected in a primary alignment."""
    bam_ptr = importlib.resources.files(resources).joinpath("reads", "hard_clipped_primary.bam")
    with importlib.resources.as_file(bam_ptr) as bam_path:
        report = grz_check.validate_bam(bam_path)

    assert report.is_valid
    assert len(report.errors) == 0

    # Check the native warnings list instead of intercepting python logs
    assert any("primary alignment(s) with hard-clipped bases" in w for w in report.warnings)


def test_hard_clipped_primary_stream(caplog):
    """A warning should be logged if hard-clipped bases are detected in a primary alignment."""
    bam_ptr = importlib.resources.files(resources).joinpath("reads", "hard_clipped_primary.bam")

    with importlib.resources.as_file(bam_ptr) as bam_path, caplog.at_level(logging.WARNING):
        validator = BamValidator()
        with open(bam_path, "rb") as f:
            ReadStream(f) >> validator

    assert "primary alignment(s) with hard-clipped bases" in caplog.text


def test_secondary():
    """A warning should be present if a secondary alignment is detected."""
    bam_ptr = importlib.resources.files(resources).joinpath("reads", "secondary.bam")
    with importlib.resources.as_file(bam_ptr) as bam_path:
        report = grz_check.validate_bam(bam_path)

    assert report.is_valid
    assert len(report.errors) == 0

    # Check for the secondary alignment warning
    assert any("secondary alignment(s)" in w for w in report.warnings)

    # Hard-clipped bases are fine in secondaries, so ensure the primary warning didn't trigger
    assert not any("primary alignment(s) with hard-clipped bases" in w for w in report.warnings)


def test_secondary_stream(caplog):
    """A warning should be logged if a secondary alignment is detected."""
    bam_ptr = importlib.resources.files(resources).joinpath("reads", "secondary.bam")

    with importlib.resources.as_file(bam_ptr) as bam_path, caplog.at_level(logging.WARNING):
        validator = BamValidator()
        with open(bam_path, "rb") as f:
            ReadStream(f) >> validator
    assert "File contains 1 secondary alignment" in caplog.text
    # hard-clipped bases are fine in secondaries
    assert "primary alignment(s) with hard-clipped bases" not in caplog.text
