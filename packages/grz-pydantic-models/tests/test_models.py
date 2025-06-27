import importlib.resources

import pytest
from grz_pydantic_models.submission.metadata import File, FileType, GrzSubmissionMetadata
from pydantic import ValidationError

from . import resources


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


def test_lab_datum():
    metadata = GrzSubmissionMetadata.model_validate_json(
        importlib.resources.files(resources)
        .joinpath("example_metadata", "wes_tumor_germline", "v1.1.4.json")
        .read_text()
    )
    with pytest.raises(ValueError, match="Long read libraries can't be paired-end."):
        metadata.donors[0].lab_data[0].library_type = "wes_lr"
