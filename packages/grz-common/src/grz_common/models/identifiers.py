from grz_pydantic_models.submission.metadata import GenomicDataCenterId, SubmitterId

from .base import IgnoringBaseModel, IgnoringBaseSettings


class IdentifiersModel(IgnoringBaseModel):
    grz: GenomicDataCenterId
    """Id of the GRZ."""

    le: SubmitterId
    """Id of the Leistungserbringer."""


class IdentifiersConfigModel(IgnoringBaseSettings):
    identifiers: IdentifiersModel
