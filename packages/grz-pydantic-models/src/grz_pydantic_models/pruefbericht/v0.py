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


class LibraryType(enum.StrEnum):
    """Sequencing method, if applicable."""

    panel = "panel"
    wes = "wes"
    wgs = "wgs"
    wgs_lr = "wgs_lr"
    none = "none"


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

    _version: str = "0.3"

    submitted_case: SubmittedCase = Field(alias="SubmittedCase")
