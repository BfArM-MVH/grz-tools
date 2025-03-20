"""Classes for parsing and validating submission metadata and files."""

from __future__ import annotations

import json
import logging
from collections.abc import Generator
from itertools import groupby
from os import PathLike
from pathlib import Path

from grz_pydantic_models.v1_1_1.metadata import (
    ChecksumType,
    File,
    FileType,
    GrzSubmissionMetadata,
    ReadOrder,
    SequenceData,
    SequencingLayout,
)
from grz_pydantic_models.v1_1_1.metadata import File as SubmissionFileMetadata
from pydantic import ValidationError

from .download import S3BotoDownloadWorker
from .fastq_validation import validate_paired_end_reads, validate_single_end_reads
from .file_operations import Crypt4GH, calculate_sha256
from .models.config import ConfigModel
from .progress_logging import FileProgressLogger
from .states import DecryptionState, EncryptionState, ValidationState
from .upload import S3BotoUploadWorker

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
                try:
                    metadata_model = GrzSubmissionMetadata(**metadata)
                except ValidationError as ve:
                    cls.__log.error("Invalid metadata format in metadata file: %s", file_path)
                    raise SystemExit(ve) from ve
                return metadata_model
        except json.JSONDecodeError as e:
            cls.__log.error("Invalid JSON format in metadata file: %s", file_path)
            raise e

    @property
    def transaction_id(self) -> str:
        """
        The index case ID of this submission
        """
        return self.content.submission.tan_g

    @property
    def files(self) -> dict[Path, SubmissionFileMetadata]:
        """
        The files linked in the metadata.

        :return: Dictionary of `file_path` -> `SubmissionFileMetadata` pairs.
            Each `file_path` refers to the relative file path from the metadata.
        """
        if self._files is not None:
            return self._files

        submission_files = {}
        for donor in self.content.donors:
            for lab_data in donor.lab_data:
                for file_data in lab_data.sequence_data.files:
                    file_path = Path(file_data.file_path)
                    if file_path.is_symlink():
                        raise ValueError(f"Provided path is a symlink which is not accepted: {file_path}")
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
                for file_data in lab_data.sequence_data.files:
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

    @staticmethod
    def validate_file_data(metadata: File, local_file_path: Path) -> Generator[str]:
        """
        Validates whether the provided file matches this metadata.

        :param metadata: Metadata model object
        :param local_file_path: Path to the actual file (resolved if symlinked)
        :return: Generator of errors
        """
        # Resolve file path
        local_file_path = local_file_path.resolve()

        # Check if path exists
        if not local_file_path.exists():
            yield f"{str(metadata.file_path)} does not exist!"
            # Return here as following tests cannot work
            return

        # Check if path is a file
        if not local_file_path.is_file():
            yield f"{str(metadata.file_path)} is not a file!"
            # Return here as following tests cannot work
            return

        # Check if the checksum is correct
        if metadata.checksum_type == "sha256":
            calculated_checksum = calculate_sha256(local_file_path)
            if metadata.file_checksum != calculated_checksum:
                yield (
                    f"{str(metadata.file_path)}: Checksum mismatch! "
                    f"Expected: '{metadata.file_checksum}', calculated: '{calculated_checksum}'."
                )
        else:
            yield (
                f"{str(metadata.file_path)}: Unsupported checksum type: {metadata.checksum_type}. "
                f"Supported types: {[e.value for e in ChecksumType]}"
            )

        # Check file size
        if metadata.file_size_in_bytes != local_file_path.stat().st_size:
            yield (
                f"{str(metadata.file_path)}: File size mismatch! "
                f"Expected: '{metadata.file_size_in_bytes}', observed: '{local_file_path.stat().st_size}'."
            )

    def validate_checksums(self, progress_log_file: str | PathLike) -> Generator[str]:
        """
        Validates the checksum of the files against the metadata and prints the errors.

        :return: Generator of errors
        """
        progress_logger = FileProgressLogger[ValidationState](log_file_path=progress_log_file)
        # cleanup log file and keep only files listed here
        progress_logger.cleanup(keep=[(file_path, file_metadata) for file_path, file_metadata in self.files.items()])
        # fields:
        # - "errors": List[str]
        # - "validation_passed": bool

        def validate_file(local_file_path, file_metadata):
            self.__log.debug("Validating '%s'...", str(local_file_path))

            # validate the file
            errors = list(self.validate_file_data(file_metadata, local_file_path))
            validation_passed = len(errors) == 0

            # return log state
            return ValidationState(errors=errors, validation_passed=validation_passed)

        for local_file_path, file_metadata in self.files.items():
            logged_state = progress_logger.get_state(
                local_file_path,
                file_metadata,
                default=validate_file,  # validate the file if the state was not calculated yet
            )

            if logged_state:
                yield from logged_state["errors"]

    def validate_sequencing_data(self, progress_log_file: str | PathLike) -> Generator[str]:
        """
        Quick-validates sequencing data linked in this submission.

        :return: Generator of errors
        """
        from .progress_logging import FileProgressLogger

        progress_logger = FileProgressLogger[ValidationState](log_file_path=progress_log_file)
        # cleanup log file and keep only files listed here
        progress_logger.cleanup(keep=[(file_path, file_metadata) for file_path, file_metadata in self.files.items()])
        # fields:
        # - "errors": List[str]
        # - "validation_passed": bool

        def find_fastq_files(sequence_data: SequenceData) -> list[File]:
            return [f for f in sequence_data.files if f.file_type == FileType.fastq]

        for donor in self.metadata.content.donors:
            for lab_data in donor.lab_data:
                sequencing_layout = lab_data.sequencing_layout
                sequence_data = lab_data.sequence_data
                # find all FASTQ files
                fastq_files = find_fastq_files(sequence_data)

                match sequencing_layout:
                    case SequencingLayout.single_end | SequencingLayout.reverse | SequencingLayout.other:
                        yield from self._validate_single_end(fastq_files, progress_logger)

                    case SequencingLayout.paired_end:
                        yield from self._validate_paired_end(fastq_files, progress_logger)

    def _validate_single_end(
        self,
        fastq_files: list[File],
        progress_logger: FileProgressLogger[ValidationState],
    ) -> Generator[str, None, None]:
        def validate_file(local_file_path, _file_metadata) -> ValidationState:
            self.__log.debug("Validating '%s'...", str(local_file_path))

            # validate the file
            errors = list(validate_single_end_reads(local_file_path))
            validation_passed = len(errors) == 0

            # return log state
            return ValidationState(
                errors=errors,
                validation_passed=validation_passed,
            )

        for fastq_file in fastq_files:
            logged_state = progress_logger.get_state(
                self.files_dir / fastq_file.file_path,
                fastq_file,
                default=validate_file,  # validate the file if the state was not calculated yet
            )
            if logged_state:
                yield from logged_state["errors"]

    def _validate_paired_end(
        self,
        fastq_files: list[File],
        progress_logger: FileProgressLogger[ValidationState],
    ) -> Generator[str, None, None]:
        key = lambda f: (f.flowcell_id, f.lane_id)
        fastq_files.sort(key=key)
        for _key, group in groupby(fastq_files, key):
            files = list(group)

            # separate R1 and R2 files
            fastq_r1_files = [f for f in files if f.read_order == ReadOrder.r1]
            fastq_r2_files = [f for f in files if f.read_order == ReadOrder.r2]

            for fastq_r1, fastq_r2 in zip(fastq_r1_files, fastq_r2_files, strict=True):
                local_fastq_r1_path = self.files_dir / fastq_r1.file_path
                local_fastq_r2_path = self.files_dir / fastq_r2.file_path

                # get saved state
                logged_state_r1 = progress_logger.get_state(
                    local_fastq_r1_path,
                    fastq_r1,
                )
                logged_state_r2 = progress_logger.get_state(
                    local_fastq_r2_path,
                    fastq_r2,
                )
                if logged_state_r1 is None or logged_state_r2 is None or logged_state_r1 != logged_state_r2:
                    # calculate state
                    errors = list(
                        validate_paired_end_reads(
                            local_fastq_r1_path,  # fastq R1
                            local_fastq_r2_path,  # fastq R2
                        )
                    )
                    validation_passed = len(errors) == 0

                    state = ValidationState(
                        errors=errors,
                        validation_passed=validation_passed,
                    )
                    # update state for both files
                    progress_logger.set_state(  # fastq R1
                        local_fastq_r1_path, fastq_r1, state
                    )
                    progress_logger.set_state(  # fastq R2
                        local_fastq_r2_path, fastq_r2, state
                    )

                    yield from state["errors"]
                else:
                    # both fastq states are equal, so simply yield one of them
                    yield from logged_state_r1["errors"]

    def encrypt(
        self,
        encrypted_files_dir: str | PathLike,
        progress_log_file: str | PathLike,
        recipient_public_key_path: str | PathLike,
        submitter_private_key_path: str | PathLike | None = None,
    ) -> EncryptedSubmission:
        """
        Encrypt this submission with a public key using Crypt4Gh

        :param encrypted_files_dir: Output directory of the encrypted files
        :param progress_log_file: Path to a log file to store the progress of the encryption process
        :param recipient_public_key_path: Path to the public key file which will be used for encryption
        :param submitter_private_key_path: Path to the private key file which will be used to sign the encryption
        :return: EncryptedSubmission instance
        """
        encrypted_files_dir = Path(encrypted_files_dir)

        if not Path(recipient_public_key_path).expanduser().is_file():
            msg = f"Public key file does not exist: {recipient_public_key_path}"
            self.__log.error(msg)
            raise FileNotFoundError(msg)
        if not submitter_private_key_path:
            self.__log.warning("No submitter private key provided, skipping signing.")
        elif not Path(submitter_private_key_path).expanduser().is_file():
            msg = f"Private key file does not exist: {submitter_private_key_path}"
            self.__log.error(msg)
            raise FileNotFoundError(msg)

        if not encrypted_files_dir.is_dir():
            self.__log.debug(
                "Creating encrypted submission files directory: %s...",
                encrypted_files_dir,
            )
            encrypted_files_dir.mkdir(mode=0o770, parents=False, exist_ok=False)

        from .progress_logging import FileProgressLogger

        progress_logger = FileProgressLogger[EncryptionState](log_file_path=progress_log_file)

        try:
            public_keys = Crypt4GH.prepare_c4gh_keys(recipient_public_key_path)
        except Exception as e:
            self.__log.error(f"Error preparing public keys: {e}")
            raise e

        for file_path, file_metadata in self.files.items():
            # encryption_successful = True
            logged_state = progress_logger.get_state(file_path, file_metadata)
            self.__log.debug("state for %s: %s", file_path, logged_state)

            encrypted_file_path = encrypted_files_dir / EncryptedSubmission.get_encrypted_file_path(
                file_metadata.file_path
            )
            encrypted_file_path.parent.mkdir(mode=0o770, parents=True, exist_ok=True)

            if (
                (logged_state is None)
                or not logged_state.get("encryption_successful", False)
                or not encrypted_file_path.is_file()
            ):
                self.__log.info(
                    "Encrypting file: '%s' -> '%s'",
                    str(file_path),
                    str(encrypted_file_path),
                )

                try:
                    Crypt4GH.encrypt_file(file_path, encrypted_file_path, public_keys)

                    self.__log.info(f"Encryption complete for {str(file_path)}. ")
                    progress_logger.set_state(
                        file_path,
                        file_metadata,
                        state=EncryptionState(encryption_successful=True),
                    )
                except Exception as e:
                    self.__log.error("Encryption failed for '%s'", str(file_path))

                    progress_logger.set_state(
                        file_path,
                        file_metadata,
                        state=EncryptionState(encryption_successful=False, errors=[str(e)]),
                    )

                    raise e
            else:
                self.__log.info(
                    "File '%s' already encrypted in '%s'",
                    str(file_path),
                    str(encrypted_file_path),
                )

        self.__log.info("File encryption completed.")

        return EncryptedSubmission(
            metadata_dir=self.metadata_dir,
            encrypted_files_dir=encrypted_files_dir,
        )


class EncryptedSubmission:
    """The encrypted counterpart to `Submission`. Handles encrypted submission data."""

    __log = log.getChild("EncryptedSubmission")

    def __init__(self, metadata_dir: str | PathLike, encrypted_files_dir: str | PathLike):
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
            encrypted_file_path = self.get_encrypted_file_path(self.encrypted_files_dir / file_path)

            retval[encrypted_file_path] = file_metadata

        return retval

    def get_metadata_file_path_and_object_id(self) -> tuple[Path, str]:
        """
        :return: tuple with the `local_file_path` and s3_object_id of the metadata file
        """
        return Path(self.metadata.file_path), str(
            Path(self.metadata.transaction_id) / "metadata" / self.metadata.file_path.name
        )

    def get_encrypted_files_and_object_id(self) -> dict[Path, str]:
        """
        :return Dictionary of `local_file_path` -> s3_object_id
        """
        retval = {}
        for local_file_path, file_metadata in self.encrypted_files.items():
            retval[local_file_path] = str(
                Path(self.metadata.transaction_id) / "files" / self.get_encrypted_file_path(file_metadata.file_path)
            )
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

        progress_logger = FileProgressLogger[DecryptionState](log_file_path=progress_log_file)

        try:
            private_key = Crypt4GH.retrieve_private_key(recipient_private_key_path)
        except Exception as e:
            self.__log.error(f"Error preparing private key: {e}")
            raise e

        for encrypted_file_path, file_metadata in self.encrypted_files.items():
            logged_state = progress_logger.get_state(encrypted_file_path, file_metadata)
            self.__log.debug("state for %s: %s", encrypted_file_path, logged_state)

            decrypted_file_path = files_dir / file_metadata.file_path
            if not decrypted_file_path.parent.is_dir():
                decrypted_file_path.parent.mkdir(mode=0o770, parents=True, exist_ok=False)

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
                    Crypt4GH.decrypt_file(encrypted_file_path, decrypted_file_path, private_key)

                    self.__log.info(f"Decryption complete for {str(encrypted_file_path)}. ")
                    progress_logger.set_state(
                        encrypted_file_path,
                        file_metadata,
                        state=DecryptionState(decryption_successful=True),
                    )
                except Exception as e:
                    self.__log.error("Decryption failed for '%s'", str(encrypted_file_path))

                    progress_logger.set_state(
                        encrypted_file_path,
                        file_metadata,
                        state=DecryptionState(decryption_successful=False, errors=[str(e)]),
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


class Worker:
    """Worker class for handling submission processing"""

    __log = log.getChild("Worker")

    def __init__(
        self,
        metadata_dir: str | PathLike,
        files_dir: str | PathLike,
        log_dir: str | PathLike,
        encrypted_files_dir: str | PathLike,
        threads: int = 1,
    ):
        """
        Initialize the worker object.

        :param metadata_dir: Path to the metadata directory
        :param files_dir: Path to the files directory
        :param log_dir: Path to the log directory
        :param encrypted_files_dir: Path to the encrypted files directory
        :param threads: Number of threads to use
        """
        self._threads = threads
        self.__log.debug("Threads: %s", self._threads)

        # metadata dir
        self.metadata_dir = Path(metadata_dir)
        self.__log.debug("Metadata directory: %s", self.metadata_dir)

        # files dir
        self.files_dir = Path(files_dir)
        self.__log.debug("Files directory: %s", self.files_dir)

        # encrypted files dir
        self.encrypted_files_dir = Path(encrypted_files_dir) if encrypted_files_dir is not None else Path()

        self.__log.info("Encrypted files directory: %s", self.encrypted_files_dir)

        # log dir
        self.log_dir = Path(log_dir)
        self.__log.info("Log directory: %s", self.log_dir)

        # create log dir if non-existent
        if not self.log_dir.is_dir():
            self.__log.debug("Creating log directory...")
            self.log_dir.mkdir(mode=0o770, parents=False, exist_ok=False)

        self.progress_file_checksum_validation = self.log_dir / "progress_checksum_validation.cjson"
        self.progress_file_sequencing_data_validation = self.log_dir / "progress_sequencing_data_validation.cjson"
        self.progress_file_encrypt = self.log_dir / "progress_encrypt.cjson"
        self.progress_file_decrypt = self.log_dir / "progress_decrypt.cjson"
        self.progress_file_upload = self.log_dir / "progress_upload.cjson"
        self.progress_file_download = self.log_dir / "progress_download.cjson"

    def parse_submission(self) -> Submission:
        """
        Reads the submission metadata and returns a Submission instance
        """
        submission = Submission(
            metadata_dir=self.metadata_dir,
            files_dir=self.files_dir,
        )
        return submission

    def parse_encrypted_submission(self) -> EncryptedSubmission:
        """
        Reads the submission metadata and returns an EncryptedSubmission instance
        """
        encrypted_submission = EncryptedSubmission(
            metadata_dir=self.metadata_dir,
            encrypted_files_dir=str(self.encrypted_files_dir),
        )
        return encrypted_submission

    def validate(self, force=False):
        """
        Validate this submission

        :param force: Force validation of already validated files
        :raises SubmissionValidationError: if the validation fails
        """
        submission = self.parse_submission()

        self.__log.info("Starting metadata validation...")
        if errors := list(submission.metadata.validate()):
            error_msg = "\n".join(["Metadata validation failed! Errors:", *errors])
            self.__log.error(error_msg)

            raise SubmissionValidationError(error_msg)
        else:
            self.__log.info("Metadata validation successful!")

        if force:
            # delete the log file
            self.progress_file_checksum_validation.unlink()

        self.__log.info("Starting checksum validation...")
        if errors := list(submission.validate_checksums(progress_log_file=self.progress_file_checksum_validation)):
            error_msg = "\n".join(["Checksum validation failed! Errors:", *errors])
            self.__log.error(error_msg)

            raise SubmissionValidationError(error_msg)
        else:
            self.__log.info("Checksum validation successful!")

        self.__log.info("Starting sequencing data validation...")
        if errors := list(
            submission.validate_sequencing_data(progress_log_file=self.progress_file_sequencing_data_validation)
        ):
            error_msg = "\n".join(["Sequencing data validation failed! Errors:", *errors])
            self.__log.error(error_msg)

            raise SubmissionValidationError(error_msg)
        else:
            self.__log.info("Sequencing data validation successful!")

    def encrypt(
        self,
        recipient_public_key_path: str | PathLike,
        submitter_private_key_path: str | PathLike | None = None,
        force=False,
    ) -> EncryptedSubmission:
        """
        Encrypt this submission with a public key using Crypt4Gh.
        :param recipient_public_key_path: Path to the public key file of the recipient.
        :param submitter_private_key_path: Path to the private key file of the submitter.
        :param force: Force encryption of already encrypted files
        :return: EncryptedSubmission instance
        """
        submission = self.parse_submission()

        if force:
            # delete the log file
            self.progress_file_encrypt.unlink()

        encrypted_submission = submission.encrypt(
            encrypted_files_dir=str(self.encrypted_files_dir),
            progress_log_file=self.progress_file_encrypt,
            recipient_public_key_path=recipient_public_key_path,
            submitter_private_key_path=submitter_private_key_path,
        )

        return encrypted_submission

    def decrypt(self, recipient_private_key_path: str | PathLike, force=False) -> Submission:
        """
        Encrypt this submission with a public key using Crypt4Gh.
        :param recipient_public_key_path: Path to the public key file of the recipient.
        :param submitter_private_key_path: Path to the private key file of the submitter.
        :param force: Force decryption of already decrypted files
        :return: EncryptedSubmission instance
        """
        encrypted_submission = self.parse_encrypted_submission()

        if force:
            # delete the log file
            self.progress_file_decrypt.unlink()

        submission = encrypted_submission.decrypt(
            files_dir=self.files_dir,
            progress_log_file=self.progress_file_decrypt,
            recipient_private_key_path=recipient_private_key_path,
        )

        return submission

    def upload(self, config: ConfigModel):
        """
        Upload an encrypted submission

        """
        upload_worker = S3BotoUploadWorker(config, status_file_path=self.progress_file_upload)

        encrypted_submission = self.parse_encrypted_submission()

        upload_worker.upload(encrypted_submission)

    def download(self, config: ConfigModel, submission_id: str):
        """
        Download an encrypted submission
        """
        download_worker = S3BotoDownloadWorker(config, status_file_path=self.progress_file_upload)

        self.__log.info("Preparing output directories...")
        download_worker.prepare_download(self.metadata_dir, self.encrypted_files_dir, self.log_dir)

        self.__log.info("Downloading metadata...")
        download_worker.download_metadata(submission_id, self.metadata_dir, metadata_file_name="metadata.json")

        self.__log.info("Downloading encrypted files...")
        download_worker.download(
            submission_id,
            EncryptedSubmission(self.metadata_dir, self.encrypted_files_dir),
        )
