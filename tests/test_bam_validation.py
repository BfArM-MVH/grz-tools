"""Tests for the bam_validation module using ValidateOperation."""

import importlib.resources
import logging
from pathlib import Path

from grz_common.pipeline.operations import ValidateOperation

from . import resources


def test_valid_hifi_bam():
    """Valid HiFi BAM files should pass validation"""
    bam_ptr = importlib.resources.files(resources).joinpath("reads", "valid_HiFi.bam")
    validator = ValidateOperation()

    with importlib.resources.as_file(bam_ptr) as bam_path:
        valid, errors = validator.validate_bam_file(Path(bam_path))

    assert valid
    assert len(errors) == 0


def test_hard_clipped_primary(caplog):
    """A warning should be logged if hard-clipped bases are detected in a primary alignment."""
    bam_ptr = importlib.resources.files(resources).joinpath("reads", "hard_clipped_primary.bam")
    validator = ValidateOperation()

    with importlib.resources.as_file(bam_ptr) as bam_path, caplog.at_level(logging.WARNING):
        valid, errors = validator.validate_bam_file(Path(bam_path))

    # Check for the new log message format from BamValidator
    assert "Detected hard-clipped bases in primary alignment" in caplog.text
    assert valid  # It is valid, just has warnings
    assert len(errors) == 0


def test_secondary(caplog):
    """A warning should be logged if a secondary alignment is detected."""
    bam_ptr = importlib.resources.files(resources).joinpath("reads", "secondary.bam")
    validator = ValidateOperation()

    with importlib.resources.as_file(bam_ptr) as bam_path, caplog.at_level(logging.WARNING):
        valid, errors = validator.validate_bam_file(Path(bam_path))

    # Check for the new log message format from BamValidator
    assert "Detected secondary alignment in BAM" in caplog.text
    # hard-clipped bases are fine in secondaries
    assert "Detected hard-clipped bases in primary alignment" not in caplog.text
    assert valid
    assert len(errors) == 0
