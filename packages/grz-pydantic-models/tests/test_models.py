import pytest
from grz_pydantic_models.submission.metadata import File, FileType
from grz_pydantic_models.submission.thresholds.v1 import (
    PercentBasesAboveQualityThreshold,
    TargetedRegionsAboveMinCoverage,
    Thresholds,
)
from pydantic import ValidationError


def test_file_paths():
    with pytest.raises(ValidationError, match="File paths must be normalized"):
        File(
            filePath="../test.bed",
            fileType=FileType.bed,
            fileChecksum="01ba4719c80b6fe911b091a7c05124b64eeece964e09c058ef8f9805daca546b",
            fileSizeInBytes=0,
        )

    with pytest.raises(ValidationError, match="File paths must be normalized"):
        File(
            filePath="files/./test.bed",
            fileType=FileType.bed,
            fileChecksum="01ba4719c80b6fe911b091a7c05124b64eeece964e09c058ef8f9805daca546b",
            fileSizeInBytes=0,
        )

    with pytest.raises(ValidationError, match="File paths must be relative"):
        File(
            filePath="/data/sensitive/target.bed",
            fileType=FileType.bed,
            fileChecksum="01ba4719c80b6fe911b091a7c05124b64eeece964e09c058ef8f9805daca546b",
            fileSizeInBytes=0,
        )


def test_thresholds_mean_read_length_default():
    """Test that meanReadLength defaults to -1 to indicate disabled/skip for grz-check."""
    threshold = Thresholds(
        meanDepthOfCoverage=100.0,
        percentBasesAboveQualityThreshold=PercentBasesAboveQualityThreshold(
            qualityThreshold=30.0,
            percentBasesAbove=85.0,
        ),
        targetedRegionsAboveMinCoverage=TargetedRegionsAboveMinCoverage(
            minCoverage=30.0,
            fractionAbove=0.8,
        ),
    )
    assert threshold.mean_read_length == -1


def test_thresholds_mean_read_length_positive():
    """Test that positive meanReadLength values work correctly."""
    threshold = Thresholds(
        meanDepthOfCoverage=100.0,
        meanReadLength=70,
        percentBasesAboveQualityThreshold=PercentBasesAboveQualityThreshold(
            qualityThreshold=30.0,
            percentBasesAbove=85.0,
        ),
        targetedRegionsAboveMinCoverage=TargetedRegionsAboveMinCoverage(
            minCoverage=30.0,
            fractionAbove=0.8,
        ),
    )
    assert threshold.mean_read_length == 70


def test_thresholds_mean_read_length_negative():
    """Test that negative meanReadLength values work correctly to indicate skip."""
    threshold = Thresholds(
        meanDepthOfCoverage=100.0,
        meanReadLength=-1,
        percentBasesAboveQualityThreshold=PercentBasesAboveQualityThreshold(
            qualityThreshold=30.0,
            percentBasesAbove=85.0,
        ),
        targetedRegionsAboveMinCoverage=TargetedRegionsAboveMinCoverage(
            minCoverage=30.0,
            fractionAbove=0.8,
        ),
    )
    assert threshold.mean_read_length == -1
