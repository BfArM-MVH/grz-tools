"""Validation package for grz-cli."""

# ruff: noqa: F401
from .bam import validate_bam
from .fastq import validate_paired_end_reads, validate_single_end_reads
