"""Classes for parsing and validating submission metadata and files."""

from __future__ import annotations

import json
import logging
from collections.abc import Generator
from os import PathLike
from pathlib import Path

from .file_operations import Crypt4GH, calculate_sha256
from .models.v1_0_0.metadata import File as SubmissionFileMetadata
from .models.v1_0_0.metadata import GrzSubmissionMetadata

log = logging.getLogger(__name__)


class SubmissionMetadata:
    """Class for reading and validating submission metadata"""

    __log = log.getChild("SubmissionMetadata")

    def __init__(self, metadata_file):
        """
        Load, parse and validate the metadata file.

        :param metadata_file: path to the metadata.json file
        :raises json.JSONDecodeError: if failed to read the metadata.json file
        :raises jsonschema.exceptions.ValidationError: if metadata does not match expected schema
        """
        self.file_path = metadata_file
        self.content = self._read_metadata(self.file_path)
        self._checksum = calculate_sha256(self.file_path, progress=False)

        self._files = None

    @classmethod
    def _read_metadata(cls, file_path: Path) -> GrzSubmissionMetadata:
        """
        Load and parse the metadata file in JSON format.

        :param file_path: Path to the metadata JSON file
        :return: Parsed metadata as a dictionary
        :raises json.JSONDecodeError: if failed to read the metadata.json file
        """
        try:
            with open(file_path, encoding="utf-8") as jsonfile:
                metadata = json.load(jsonfile)
                metadata_model = GrzSubmissionMetadata(**metadata)
                return metadata_model
        except json.JSONDecodeError as e:
            cls.__log.error("Invalid JSON format in metadata file: %s", file_path)
            raise e

    @property
    def index_case_id(self) -> str:
        """
        The index case ID of this submission
        """
        return self.content.submission.tan_g

    @property
    def files(self) -> dict[Path, SubmissionFileMetadata]:
        """
        The files liked in the metadata.

        :return: Dictionary of `file_path` -> `SubmissionFileMetadata` pairs.
            Each `file_path` refers to the relative file path from the metadata.
        """
        if self._files is not None:
            return self._files

        submission_files = {}
        for donor in self.content.donors:
            for lab_data in donor.lab_data:
                for sequence_data in lab_data.sequence_data:
                    for file_data in sequence_data.files:
                        file_path = Path(file_data.file_path)
                        if file_path.is_symlink():
                            raise ValueError(
                                f"Provided path is a symlink which is not accepted: {file_path}"
                            )
                        else:
                            submission_files[file_path] = file_data

        self._files = submission_files
        return self._files

    def validate(self) -> Generator[str]:
        """
        Validates this submission's metadata (content).

        :return: Generator of errors
        """
        submission_files: dict[str | PathLike, SubmissionFileMetadata] = {}
        for donor in self.content.donors:
            for lab_data in donor.lab_data:
                for sequence_data in lab_data.sequence_data:
                    for file_data in sequence_data.files:
                        # check if file is already registered
                        file_path = Path(file_data.file_path)
                        if other_metadata := submission_files.get(file_path):
                            # check if metadata matches
                            if file_data != other_metadata:
                                yield f"{file_data.file_path}: Different metadata for the same path observed!"

                            # check if FASTQ data was already linked in another submission
                            if file_data.file_type == "fastq":
                                yield f"{file_data.file_path}: FASTQ file already linked in another submission!"
                            if file_data.file_type == "bam":
                                yield f"{file_data.file_path}: BAM file already linked in another submission!"
                        else:
                            submission_files[file_path] = file_data

    @property
    def checksum(self) -> str:
        """
        Checksum of the metadata file
        """
        return self._checksum


class Submission:
    """Class for handling submission data"""

    __log = log.getChild("Submission")

    def __init__(self, metadata_dir: str | PathLike, files_dir: str | PathLike):
        """
        Initialize the submission object.

        :param metadata_dir: Path to the metadata directory
        :param files_dir: Path to the files directory
        """
        self.metadata_dir = Path(metadata_dir)
        self.files_dir = Path(files_dir)

        self.metadata = SubmissionMetadata(self.metadata_dir / "metadata.json")

    @property
    def files(self) -> dict[Path, SubmissionFileMetadata]:
        """
        The files liked in the metadata.

        :return: Dictionary of `local_file_path` -> `SubmissionFileMetadata` pairs.
        """
        retval = {}
        for file_path, file_metadata in self.metadata.files.items():
            local_file_path = self.files_dir / file_path

            retval[local_file_path] = file_metadata

        return retval

    def validate_checksums(self, progress_log_file: str | PathLike) -> Generator[str]:
        """
        Validates the checksum of the files against the metadata and prints the errors.

        :return: Generator of errors
        """
        from .progress_logging import FileProgressLogger

        progress_logger = FileProgressLogger(log_file_path=progress_log_file)
        # cleanup log file and keep only files listed here
        progress_logger.cleanup(
            keep=[
                (file_path, file_metadata)
                for file_path, file_metadata in self.files.items()
            ]
        )
        # fields:
        # - "errors": List[str]
        # - "validation_passed": bool

        for local_file_path, file_metadata in self.files.items():
            logged_state = progress_logger.get_state(local_file_path, file_metadata)

            # determine if we can skip the verification
            if logged_state is None:
                self.__log.debug("State for %s not calculated yet", local_file_path)
            elif not logged_state.get("validation_passed", False):
                errors = logged_state.get("errors", [])
                yield from errors

                # skip re-verification
                continue
            else:
                self.__log.debug(
                    "Validation for %s already passed, skipping...",
                    str(local_file_path),
                )

                # skip re-verification
                continue

            self.__log.debug("Validating '%s'...", str(local_file_path))
            # validate the file
            errors = list(file_metadata.validate_data(local_file_path))
            validation_passed = len(errors) == 0

            # log state
            progress_logger.set_state(
                local_file_path,
                file_metadata,
                state={
                    "errors": errors,
                    "validation_passed": validation_passed,
                },
            )

            yield from errors


class EncryptedSubmission:
    """The encrypted counterpart to `Submission`. Handles encrypted submission data."""

    __log = log.getChild("EncryptedSubmission")

    def __init__(
        self, metadata_dir: str | PathLike, encrypted_files_dir: str | PathLike
    ):
        """
        Initialize the encrypted submission object.

        :param metadata_dir: Path to the metadata directory
        :param encrypted_files_dir: Path to the encrypted files directory
        """
        self.metadata_dir = Path(metadata_dir)
        self.encrypted_files_dir = Path(encrypted_files_dir)

        self.metadata = SubmissionMetadata(self.metadata_dir / "metadata.json")

    @property
    def encrypted_files(self):
        """
        The encrypted files liked in the metadata.

        :return: Dictionary of `local_file_path` -> `SubmissionFileMetadata` pairs.
        """
        retval = {}
        for file_path, file_metadata in self.metadata.files.items():
            encrypted_file_path = self.get_encrypted_file_path(
                self.encrypted_files_dir / file_path
            )

            retval[encrypted_file_path] = file_metadata

        return retval

    @staticmethod
    def get_encrypted_file_path(file_path: str | PathLike) -> Path:
        """
        Return the path to the encrypted file based on the original file path,
        with additional extension'.c4gh'.
        """
        p = Path(file_path)
        return p.with_suffix(p.suffix + ".c4gh")

    @staticmethod
    def get_encryption_header_path(file_path: str | PathLike) -> Path:
        """
        Return the path to the encryption header file based on the original file path,
        with additional extension'.c4gh_header'.
        """
        p = Path(file_path)
        return p.with_suffix(p.suffix + ".c4gh_header")

    def decrypt(
        self,
        files_dir: str | PathLike,
        progress_log_file: str | PathLike,
        recipient_private_key_path: str | PathLike,
    ) -> Submission:
        """
        Decrypt this encrypted submission with a private key using Crypt4Gh

        :param files_dir: Output directory of the decrypted files
        :param progress_log_file: Path to a log file to store the progress of the decryption process
        :param recipient_private_key_path: Path to the private key file which will be used for decryption
        :return: Submission instance
        """
        files_dir = Path(files_dir)

        if not files_dir.is_dir():
            self.__log.debug(
                "Creating decrypted submission files directory: %s...",
                files_dir,
            )
            files_dir.mkdir(mode=0o770, parents=False, exist_ok=False)

        from .progress_logging import FileProgressLogger

        progress_logger = FileProgressLogger(log_file_path=progress_log_file)

        try:
            private_key = Crypt4GH.retrieve_private_key(recipient_private_key_path)
        except Exception as e:
            self.__log.error(f"Error preparing private key: {e}")
            raise e

        for encrypted_file_path, file_metadata in self.encrypted_files.items():
            logged_state = progress_logger.get_state(encrypted_file_path, file_metadata)
            self.__log.debug("state for %s: %s", encrypted_file_path, logged_state)

            decrypted_file_path = files_dir / file_metadata.file_path

            if (
                (logged_state is None)
                or not logged_state.get("decryption_successful", False)
                or not decrypted_file_path.is_file()
            ):
                self.__log.info(
                    "Decrypting file: '%s' -> '%s'",
                    str(encrypted_file_path),
                    str(decrypted_file_path),
                )

                try:
                    Crypt4GH.decrypt_file(
                        encrypted_file_path, decrypted_file_path, private_key
                    )

                    self.__log.info(
                        f"Decryption complete for {str(encrypted_file_path)}. "
                    )
                    progress_logger.set_state(
                        encrypted_file_path,
                        file_metadata,
                        state={"decryption_successful": True},
                    )
                except Exception as e:
                    self.__log.error(
                        "Decryption failed for '%s'", str(encrypted_file_path)
                    )

                    progress_logger.set_state(
                        encrypted_file_path,
                        file_metadata,
                        state={"decryption_successful": False, "error": str(e)},
                    )

                    raise e
            else:
                self.__log.info(
                    "File '%s' already decrypted in '%s'",
                    str(encrypted_file_path),
                    str(decrypted_file_path),
                )

        self.__log.info("File decryption completed.")

        return Submission(
            metadata_dir=self.metadata_dir,
            files_dir=files_dir,
        )


class SubmissionValidationError(Exception):
    """Exception raised when validation of a submission fails"""

    pass
