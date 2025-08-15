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

    secondary_warned = False
    hard_clipped_warned = False
    # need until_eof since BAM won't have an index
    for read in bam_file.fetch(until_eof=True):
        if not secondary_warned and read.is_secondary:
            log.warning(
                f"Detected a secondary alignment in BAM file '{bam_path}'. "
                "Please consider filtering them to save upload bandwidth and storage space."
            )
            secondary_warned = True

        if not hard_clipped_warned and not read.is_secondary:
            hard_clipped_count = read.get_cigar_stats()[0][5]
            if hard_clipped_count:
                log.warning(
                    f"Detected hard-clipped base(s) in a primary alignment in BAM file '{bam_path}'. "
                    "This is a loss of information from the raw reads."
                )
                hard_clipped_warned = True

    yield from errors
