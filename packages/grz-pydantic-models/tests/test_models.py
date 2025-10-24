import pytest
from grz_pydantic_models.submission.metadata import File, FileType
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
