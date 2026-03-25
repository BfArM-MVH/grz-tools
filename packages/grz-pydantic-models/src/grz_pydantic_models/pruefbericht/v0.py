"""
First draft of Prüfbericht schema
"""

import datetime
import enum

from pydantic import Field

from ..common import StrictBaseModel
from ..submission.metadata import (
    ClinicalDataNodeId,
    CoverageType,
    DiseaseType,
    GenomicDataCenterId,
    SubmissionType,
    SubmitterId,
    Tan,
)


class DataCategory(enum.StrEnum):
    """Type of submission."""

    clinical = "clinical"
    genomic = "genomic"


import enum


class LibraryType(enum.StrEnum):
    """Sequencing method, if applicable."""

    panel = "panel"
    wes = "wes"
    wgs = "wgs"
    wgs_lr = "wgs_lr"
    none = "none"

    @property
    def reimbursement_priority(self) -> int:
        """Reimbursement priority — higher value = more expensive."""
        priorities = {
            LibraryType.none: 0,
            LibraryType.panel: 1,
            LibraryType.wes: 2,
            LibraryType.wgs: 3,
            LibraryType.wgs_lr: 4,
        }
        return priorities[self]

    @classmethod
    def most_expensive(cls, library_types: set[str]) -> "LibraryType":
        """
        Return the LibraryType with the highest reimbursement value from a set of strings.

        Args:
            library_types: Set of library type strings.

        Returns:
            The LibraryType with the highest reimbursement priority.

        Raises:
            ValueError: If no valid LibraryType values are found in the input set.
        """
        valid_types = {cls(lt) for lt in library_types if lt in cls._value2member_map_}

        if not valid_types:
            raise ValueError(
                f"Submission contained ONLY library types ({', '.join(library_types)}) that cannot be "
                f"submitted in the Prüfbericht. Valid types are: {', '.join(cls._value2member_map_)}."
            )

        return max(valid_types, key=lambda lt: lt.reimbursement_priority)


class SubmittedCase(StrictBaseModel):
    """A single submission to a GRZ.

    For a description of fields common to the submission metadata, see the
    model for the submission metadata.
    """

    submission_date: datetime.date

    submission_type: SubmissionType

    tan: Tan
    """
    T-VNk or T-VNg.
    """

    submitter_id: SubmitterId

    data_node_id: GenomicDataCenterId | ClinicalDataNodeId

    disease_type: DiseaseType

    data_category: DataCategory

    library_type: LibraryType

    coverage_type: CoverageType

    data_quality_check_passed: bool
    """
    Whether the quality check at the data hub passed according to the standards.
    """


class Pruefbericht(StrictBaseModel):
    """Quality control report submitted to BfArM for each submission to a KDK/GRZ."""

    _version: str = "0.4"

    submitted_case: SubmittedCase = Field(alias="SubmittedCase")
