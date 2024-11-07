"""Classes for parsing and validating submission metadata and files."""

from __future__ import annotations

import json
import logging
from collections.abc import Generator
from functools import partial
from os import PathLike
from pathlib import Path

from .download import S3BotoDownloadWorker
from .file_operations import Crypt4GH, calculate_sha256
from .models.config import Backend, ConfigModel
from .models.v1_0_0.metadata import File as SubmissionFileMetadata
from .models.v1_0_0.metadata import GrzSubmissionMetadata
from .progress_logging import FileProgressLogger
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

        progress_logger = FileProgressLogger(log_file_path=progress_log_file)

        try:
            public_keys = Crypt4GH.prepare_c4gh_keys(recipient_public_key_path)
        except Exception as e:
            self.__log.error(f"Error preparing public keys: {e}")
            raise e

        self._parallel_encrypt(public_keys, encrypted_files_dir, progress_logger)

        self.__log.info("File encryption completed.")

        return EncryptedSubmission(
            metadata_dir=self.metadata_dir,
            encrypted_files_dir=encrypted_files_dir,
        )

    def _parallel_encrypt(
        self,
        public_keys,
        encrypted_files_dir: Path,
        progress_logger: FileProgressLogger,
    ):
        files_to_encrypt = list(
            enumerate(
                self._determine_files_to_encrypt(encrypted_files_dir, progress_logger)
            )
        )

        from tqdm.contrib.concurrent import process_map

        for file_path, file_metadata, state in process_map(
            partial(_process, public_keys), files_to_encrypt
        ):
            progress_logger.set_state(file_path, file_metadata, state=state)

    def _determine_files_to_encrypt(
        self, encrypted_files_dir, progress_logger
    ) -> Generator[ProcessItem, None, None]:
        for file_path, file_metadata in self.files.items():
            logged_state = progress_logger.get_state(file_path, file_metadata)
            self.__log.debug("state for %s: %s", file_path, logged_state)

            encrypted_file_path = (
                encrypted_files_dir
                / EncryptedSubmission.get_encrypted_file_path(file_metadata.file_path)
            )
            encrypted_file_path.parent.mkdir(mode=0o770, parents=True, exist_ok=True)

            if (
                (logged_state is None)
                or not logged_state.get("encryption_successful", False)
                or not encrypted_file_path.is_file()
            ):
                yield (
                    file_path,
                    encrypted_file_path,
                    file_metadata,
                    logged_state,
                    self.__log,
                )
            else:
                self.__log.info(
                    "File '%s' already encrypted in '%s'",
                    str(file_path),
                    str(encrypted_file_path),
                )


type ProcessItem = tuple[Path, Path, SubmissionFileMetadata, dict, logging.Logger]
type ProcessResult = tuple[Path, SubmissionFileMetadata, dict]


def _process(public_keys, item: tuple[int, ProcessItem]) -> ProcessResult:
    idx, (file_path, encrypted_file_path, file_metadata, logged_state, logger) = item
    logger.debug("state for %s: %s", file_path, logged_state)

    logger.info(
        "Encrypting file: '%s' -> '%s'",
        str(file_path),
        str(encrypted_file_path),
    )

    try:
        Crypt4GH.encrypt_file(
            file_path,
            encrypted_file_path,
            public_keys,
            tqdm_kwargs=dict(position=idx, leave=True),
        )

        logger.info(f"Encryption complete for {str(file_path)}. ")
        return file_path, file_metadata, {"encryption_successful": True}
    except Exception as e:
        logger.error("Encryption failed for '%s'", str(file_path))

        return (
            file_path,
            file_metadata,
            {"encryption_successful": False, "error": str(e)},
        )


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


class Worker:
    """Worker class for handling submission processing"""

    __log = log.getChild("Worker")

    def __init__(
        self,
        working_dir: str | PathLike | None = None,
        metadata_dir: str | PathLike | None = None,
        files_dir: str | PathLike | None = None,
        encrypted_files_dir: str | PathLike | None = None,
        log_dir: str | PathLike | None = None,
    ):
        """
        Initialize the worker object.

        :param working_dir: Path to the working directory
        :param metadata_dir: Path to the metadata directory
        :param files_dir: Path to the files directory
        :param encrypted_files_dir: Path to the encrypted files directory
        :param log_dir: Path to the log directory
        """
        self.working_dir = Path(working_dir) if working_dir is not None else Path.cwd()

        self.__log.debug("Working directory: %s", self.working_dir)

        # metadata dir
        self.metadata_dir = (
            Path(metadata_dir)
            if metadata_dir is not None
            else self.working_dir / "metadata"
        )
        self.__log.debug("Metadata directory: %s", self.metadata_dir)

        # files dir
        self.files_dir = (
            Path(files_dir) if files_dir is not None else self.working_dir / "files"
        )
        self.__log.debug("Files directory: %s", self.files_dir)

        # encrypted files dir
        self.encrypted_files_dir = (
            Path(encrypted_files_dir)
            if encrypted_files_dir is not None
            else self.working_dir / "encrypted_files"
        )
        self.__log.info("Encrypted files directory: %s", self.encrypted_files_dir)

        # log dir
        self.log_dir = (
            Path(log_dir) if log_dir is not None else self.working_dir / "logs"
        )
        self.__log.info("Log directory: %s", self.log_dir)

        # create log dir if non-existent
        if not self.log_dir.is_dir():
            self.__log.debug("Creating log directory...")
            self.log_dir.mkdir(mode=0o770, parents=False, exist_ok=False)

        self.progress_file_checksum = self.log_dir / "progress_checksum.cjson"
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
            encrypted_files_dir=self.encrypted_files_dir,
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
            self.progress_file_checksum.unlink()
        self.__log.info("Starting checksum validation...")
        if errors := list(
            submission.validate_checksums(progress_log_file=self.progress_file_checksum)
        ):
            error_msg = "\n".join(["Checksum validation failed! Errors:", *errors])
            self.__log.error(error_msg)

            raise SubmissionValidationError(error_msg)
        else:
            self.__log.info("Checksum validation successful!")

        # TODO: validate FASTQ

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
            encrypted_files_dir=self.encrypted_files_dir,
            progress_log_file=self.progress_file_encrypt,
            recipient_public_key_path=recipient_public_key_path,
            submitter_private_key_path=submitter_private_key_path,
        )

        return encrypted_submission

    def decrypt(
        self, recipient_private_key_path: str | PathLike, force=False
    ) -> Submission:
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
        if config.s3_options.backend == Backend.s3cmd:
            raise NotImplementedError()
        else:
            upload_worker = S3BotoUploadWorker(
                config, status_file_path=self.progress_file_upload
            )

        encrypted_submission = self.parse_encrypted_submission()

        upload_worker.upload(encrypted_submission)

    def download(self, config: ConfigModel):
        """
        Download an encrypted submission

        """
        if config.s3_options.backend == Backend.s3cmd:
            raise NotImplementedError()
        else:
            download_worker = S3BotoDownloadWorker(
                config, status_file_path=self.progress_file_upload
            )

        submission_id = self.metadata_dir.parent.name
        submission_dir = download_worker.prepare_download(self.metadata_dir)

        log.info("Prepared submission directory: %s", submission_dir)

        download_worker.download(
            submission_id,
            self.metadata_dir,
            self.encrypted_files_dir,
            metadata_file_name="metadata.json",
        )
