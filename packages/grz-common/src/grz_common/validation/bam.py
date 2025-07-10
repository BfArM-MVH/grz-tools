"""
Validation methods for BAM files

e.g. warn on headers present, which can contain personally identifying
information
"""

import logging
from collections.abc import Generator
from os import PathLike

import pysam

log = logging.getLogger(__name__)


def validate_bam(bam_path: str | PathLike) -> Generator[str]:
    """Check if a BAM file is valid for submission"""
    errors: list[str] = []

    # disable pysam's SQ header presence enforcement (check_sq)
    #  these are used to optionally sort alignments by a custom contig order,
    #  which is unnecessary in unmapped BAM files
    bam_file = pysam.AlignmentFile(str(bam_path), mode="rb", check_sq=False)

    header = bam_file.header.to_dict()
    # ignore the HD key as it is almost always present and extremely
    # unlikely to accidentally contain identifying information
    concerning_keys = header.keys() - {"HD"}
    if concerning_keys:
        log.warning(f"Detected a header in BAM file '{bam_path}', ensure it contains no private information!")

    yield from errors
