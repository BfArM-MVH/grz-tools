"""Tests for the bam_validation module."""

import importlib.resources

from grz_cli.bam_validation import validate_bam

from . import resources


def test_valid_hifi_bam():
    """Valid HiFi BAM files should pass validation"""
    bam_ptr = importlib.resources.files(resources).joinpath("reads", "valid_HiFi.bam")
    with importlib.resources.as_file(bam_ptr) as bam_path:
        errors = list(validate_bam(bam_path))
    assert len(errors) == 0
