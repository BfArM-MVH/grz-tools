from __future__ import annotations

from typing import Annotated

from pydantic import Field

from ...common import StrictBaseModel

__all__ = [
    "PercentBasesAboveQualityThreshold",
    "TargetedRegionsAboveMinCoverage",
    "Thresholds",
]


class PercentBasesAboveQualityThreshold(StrictBaseModel):
    """Fraction of bases above a given base quality threshold."""

    quality_threshold: Annotated[float, Field(strict=True, ge=0.0)]
    percent_bases_above: Annotated[float, Field(strict=True, ge=0.0, le=100.0)]


class TargetedRegionsAboveMinCoverage(StrictBaseModel):
    """Fraction of targeted regions above a given minimum coverage."""

    min_coverage: Annotated[float, Field(strict=True, ge=0.0)]
    fraction_above: Annotated[float, Field(strict=True, ge=0.0, le=1.0)]


class Thresholds(StrictBaseModel):
    """Coverage and quality thresholds for a sequencing configuration."""

    mean_depth_of_coverage: Annotated[float, Field(strict=True, ge=0.0)]
    mean_read_length: Annotated[int, Field(strict=True, ge=0)] = 0
    percent_bases_above_quality_threshold: PercentBasesAboveQualityThreshold
    targeted_regions_above_min_coverage: TargetedRegionsAboveMinCoverage
