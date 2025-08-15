"""Tests for the bam_validation module."""

import importlib.resources
import logging

from grz_common.validation.bam import validate_bam

from . import resources


def test_valid_hifi_bam():
    """Valid HiFi BAM files should pass validation"""
    bam_ptr = importlib.resources.files(resources).joinpath("reads", "valid_HiFi.bam")
    with importlib.resources.as_file(bam_ptr) as bam_path:
        errors = list(validate_bam(bam_path))
    assert len(errors) == 0


def test_hard_clipped_primary(caplog):
    """A warning should be logged if hard-clipped bases are detected in a primary alignment."""
    bam_ptr = importlib.resources.files(resources).joinpath("reads", "hard_clipped_primary.bam")
    with importlib.resources.as_file(bam_ptr) as bam_path, caplog.at_level(logging.WARNING):
        errors = list(validate_bam(bam_path))
    assert "Detected hard-clipped base(s) in a primary alignment" in caplog.text
    assert len(errors) == 0


def test_secondary(caplog):
    """A warning should be logged if a secondary alignment is detected."""
    bam_ptr = importlib.resources.files(resources).joinpath("reads", "secondary.bam")
    with importlib.resources.as_file(bam_ptr) as bam_path, caplog.at_level(logging.WARNING):
        errors = list(validate_bam(bam_path))
    assert "Detected a secondary alignment in BAM" in caplog.text
    # hard-clipped bases are fine in secondaries
    assert "Detected hard-clipped base(s) in a primary alignment" not in caplog.text
    assert len(errors) == 0
