import importlib.resources
import itertools
import json

import pytest
from grz_pydantic_models.submission.metadata.v1 import GrzSubmissionMetadata
from pydantic import ValidationError

from . import resources


@pytest.mark.parametrize(
    "dataset,version", itertools.product(["panel", "wes_tumor_germline", "wgs_tumor_germline"], ["1.1.1", "1.1.4"])
)
def test_examples(dataset: str, version: str):
    metadata_str = (
        importlib.resources.files(resources).joinpath("example_metadata", dataset, f"v{version}.json").read_text()
    )
    GrzSubmissionMetadata.model_validate_json(metadata_str)


def test_example_wgs_lr():
    metadata_str = (
        importlib.resources.files(resources).joinpath("example_metadata", "wgs_lr", "v1.1.4.json").read_text()
    )
    GrzSubmissionMetadata.model_validate_json(metadata_str)


def test_invalid_short_read_submission_with_bam():
    """BAM files should only be allowed in *_lr lab data"""
    metadata = json.loads(
        importlib.resources.files(resources)
        .joinpath("example_metadata", "wgs_tumor_germline", "v1.1.4.json")
        .read_text()
    )
    # add a BAM file
    metadata["donors"][0]["labData"][0]["sequenceData"]["files"].append(
        {
            "filePath": "donor_001/HV5TMDSX7-1-IDUDI0034_S1_L001_R1_001.bam",
            "fileType": "bam",
            "checksumType": "sha256",
            "fileChecksum": "9e87eabc18b726a94a3ffbd8d84df662388bec07b8e3d501ee6a43309c6d43fd",
            "fileSizeInBytes": 129174728987,
            "readLength": 151,
        }
    )

    with pytest.raises(ValidationError):
        GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))
