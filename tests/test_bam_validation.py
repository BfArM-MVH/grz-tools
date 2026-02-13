"""Tests for the bam_validation module using BamValidator."""

import importlib.resources
import logging

from grz_common.pipeline.components.validation import BamValidator

from . import resources


def test_valid_hifi_bam():
    """Valid HiFi BAM files should pass validation"""
    bam_ptr = importlib.resources.files(resources).joinpath("reads", "valid_HiFi.bam")

    with importlib.resources.as_file(bam_ptr) as bam_path:
        with open(bam_path, "rb") as f:
            with BamValidator(f) as validator:
                while validator.read(1024 * 1024):
                    pass
                # validate() is called on __exit__


def test_hard_clipped_primary(caplog):
    """A warning should be logged if hard-clipped bases are detected in a primary alignment."""
    bam_ptr = importlib.resources.files(resources).joinpath("reads", "hard_clipped_primary.bam")

    with importlib.resources.as_file(bam_ptr) as bam_path, caplog.at_level(logging.WARNING):
        with open(bam_path, "rb") as f:
            with BamValidator(f) as validator:
                while validator.read(1024 * 1024):
                    pass

    assert "Detected hard-clipped bases in primary alignment" in caplog.text


def test_secondary(caplog):
    """A warning should be logged if a secondary alignment is detected."""
    bam_ptr = importlib.resources.files(resources).joinpath("reads", "secondary.bam")

    with importlib.resources.as_file(bam_ptr) as bam_path, caplog.at_level(logging.WARNING):
        with open(bam_path, "rb") as f:
            with BamValidator(f) as validator:
                while validator.read(1024 * 1024):
                    pass

    assert "Detected secondary alignment in BAM" in caplog.text
    # hard-clipped bases are fine in secondaries
    assert "Detected hard-clipped bases in primary alignment" not in caplog.text
