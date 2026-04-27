"""Tests for the bam_validation module using BamValidator."""

import importlib.resources
import logging

from grz_common.pipeline.components import ReadStream
from grz_common.pipeline.components.validation import BamValidator

from . import resources


def test_valid_hifi_bam():
    """Valid HiFi BAM files should pass validation"""
    bam_ptr = importlib.resources.files(resources).joinpath("reads", "valid_HiFi.bam")

    with importlib.resources.as_file(bam_ptr) as bam_path:
        validator = BamValidator()
        with open(bam_path, "rb") as f:
            ReadStream(f) >> validator


def test_hard_clipped_primary(caplog):
    """A warning should be logged if hard-clipped bases are detected in a primary alignment."""
    bam_ptr = importlib.resources.files(resources).joinpath("reads", "hard_clipped_primary.bam")

    with importlib.resources.as_file(bam_ptr) as bam_path, caplog.at_level(logging.WARNING):
        validator = BamValidator()
        with open(bam_path, "rb") as f:
            ReadStream(f) >> validator

    assert "primary alignment(s) with hard-clipped bases" in caplog.text


def test_secondary(caplog):
    """A warning should be logged if a secondary alignment is detected."""
    bam_ptr = importlib.resources.files(resources).joinpath("reads", "secondary.bam")

    with importlib.resources.as_file(bam_ptr) as bam_path, caplog.at_level(logging.WARNING):
        validator = BamValidator()
        with open(bam_path, "rb") as f:
            ReadStream(f) >> validator
    assert "File contains 1 secondary alignment" in caplog.text
    # hard-clipped bases are fine in secondaries
    assert "primary alignment(s) with hard-clipped bases" not in caplog.text
