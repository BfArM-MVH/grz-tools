"""Worker class for handling submission processing."""

from __future__ import annotations

import logging
import shutil
import json
from dataclasses import dataclass
from datetime import date, datetime
from decimal import Decimal, ROUND_HALF_UP
from os import PathLike
from pathlib import Path
from typing import Any


from grz_db.errors import SubmissionNotFoundError
from grz_db.models.author import Author
from grz_db.models.submission import Donor
from grz_pydantic_models.submission.metadata import (
    GrzSubmissionMetadata,
    Relation,
    REDACTED_TAN
)
from grzctl.commands.db.cli import (
    get_submission_db_instance,
    redact_metadata,
    diff_metadata,
    diff_donor,
    SubmissionData,
    DonorDiff,
    DonorData
)

from grzctl.models.config import (
    DbConfig,
    ListConfig,
    DownloadConfig
)

from ..models.identifiers import IdentifiersModel
from ..models.s3 import S3Options
from ..progress import EncryptionState, FileProgressLogger, ValidationState
from ..transfer import init_s3_client
from ..validation import UserInterruptException
from .download import S3BotoDownloadWorker
from .submission import EncryptedSubmission, Submission, SubmissionValidationError
from .upload import S3BotoUploadWorker

log = logging.getLogger(__name__)

# TODO: move to grz_db.models.author?
def _initialise_author_from_dbconfig(db_config: DbConfig) -> Author:
    author = None
    if db_config.author:
        key_path = Path(db_config.author.private_key_path)
        if not key_path.exists():
            raise FileNotFoundError(f"Author private key not found at: {key_path}")

        author = Author(
            name=db_config.author.name,
            private_key_bytes=key_path.read_bytes(),
            private_key_passphrase=db_config.author.private_key_passphrase,
        )
    return Author

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
            metadata_dir=self.metadata_dir, encrypted_files_dir=str(self.encrypted_files_dir), log_dir=self.log_dir
        )
        return encrypted_submission

    def validate(self, identifiers: IdentifiersModel, force=False, with_grz_check=True):
        """
        Validate this submission

        :param identifiers: IdentifiersModel containing GRZ and LE identifiers.
        :param force: Force validation of already validated files
        :param with_grz_check: If True, use the grz-check tool for validation.
        :raises SubmissionValidationError: if the validation fails
        """
        submission = self.parse_submission()

        self.__log.info("Starting metadata validation...")
        if errors := list(submission.metadata.validate(identifiers)):
            error_msg = "\n".join(["Metadata validation failed! Errors:", *errors])
            self.__log.error(error_msg)

            raise SubmissionValidationError(error_msg)
        else:
            self.__log.info("Metadata validation successful!")

        if force:
            # delete the log files if they exist
            self.progress_file_checksum_validation.unlink(missing_ok=True)
            self.progress_file_sequencing_data_validation.unlink(missing_ok=True)

        have_grz_check = shutil.which("grz-check") is not None
        if with_grz_check and have_grz_check:
            try:
                self.__log.info("Starting file validation with `grz-check`...")
                errors = list(
                    submission.validate_files_with_grz_check(
                        checksum_progress_file=self.progress_file_checksum_validation,
                        seq_data_progress_file=self.progress_file_sequencing_data_validation,
                        threads=self._threads,
                    )
                )
                if errors:
                    error_msg = "\n".join(["File validation failed! Errors:", *errors])
                    self.__log.error(error_msg)
                    raise SubmissionValidationError(error_msg)
                else:
                    self.__log.info("File validation successful!")
                return
            except UserInterruptException as e:
                error_msg = "Validation was cancelled by the user and is incomplete."
                self.__log.error(error_msg)
                raise SubmissionValidationError(error_msg) from e

        # Fallback validation
        self.__log.info("Starting checksum validation (fallback)...")
        if errors := list(
            submission._validate_checksums_fallback(progress_log_file=self.progress_file_checksum_validation)
        ):
            error_msg = "\n".join(["Checksum validation failed! Errors:", *errors])
            self.__log.error(error_msg)
            raise SubmissionValidationError(error_msg)
        else:
            self.__log.info("Checksum validation successful!")

        self.__log.info("Starting sequencing data validation (fallback)...")
        if errors := list(
            submission._validate_sequencing_data_fallback(
                progress_log_file=self.progress_file_sequencing_data_validation
            )
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
        force: bool = False,
        check_validation_logs: bool = True,
    ) -> EncryptedSubmission:
        """
        Encrypt this submission with a public key using Crypt4Gh.
        :param recipient_public_key_path: Path to the public key file of the recipient.
        :param submitter_private_key_path: Path to the private key file of the submitter.
        :param force: Force encryption of already encrypted files
        :param check_validation_logs: Check validation logs before encrypting.
        :return: EncryptedSubmission instance
        """
        submission = self.parse_submission()

        if check_validation_logs:
            checksum_progress_logger = FileProgressLogger[ValidationState](self.progress_file_checksum_validation)
            seq_data_progress_logger = FileProgressLogger[ValidationState](
                self.progress_file_sequencing_data_validation
            )
            unvalidated_files = []

            self.__log.info("Verifying validation status of all submission files…")
            for file_path, file_metadata in submission.files.items():
                checksum_state = checksum_progress_logger.get_state(file_path, file_metadata)
                checksum_passed = checksum_state and checksum_state.get("validation_passed", False)

                seq_data_passed = True  # assume true for non-sequence files
                if file_metadata.file_type in {"fastq", "bam"}:
                    seq_data_state = seq_data_progress_logger.get_state(file_path, file_metadata)
                    seq_data_passed = seq_data_state is not None and seq_data_state.get("validation_passed", False)

                if not (checksum_passed and seq_data_passed):
                    unvalidated_files.append(str(file_path))

            if unvalidated_files:
                failed_files = "\n - ".join(unvalidated_files)
                error_msg = (
                    "Will not encrypt, as the following files were not successfully validated:\n"
                    f"{failed_files}\n"
                    "Please re-run the 'validate' command and try again."
                )
                self.__log.error(error_msg)
                raise SubmissionValidationError(error_msg)

            self.__log.info("All files verified as successfully validated.")

        if force:
            # delete the log file if it exists
            self.progress_file_encrypt.unlink(missing_ok=True)

        encrypted_submission = submission.encrypt(
            encrypted_files_dir=str(self.encrypted_files_dir),
            progress_log_file=self.progress_file_encrypt,
            recipient_public_key_path=recipient_public_key_path,
            submitter_private_key_path=submitter_private_key_path,
            force=force,
        )

        return encrypted_submission

    def decrypt(self, recipient_private_key_path: str | PathLike, force: bool = False) -> Submission:
        """
        Encrypt this submission with a public key using Crypt4Gh.
        :param recipient_private_key_path: Path to the private key file of the recipient.
        :param force: Force decryption of already decrypted files
        :return: EncryptedSubmission instance
        """
        encrypted_submission = self.parse_encrypted_submission()

        if force:
            # delete the log file if it exists
            self.progress_file_decrypt.unlink(missing_ok=True)

        submission = encrypted_submission.decrypt(
            files_dir=self.files_dir,
            progress_log_file=self.progress_file_decrypt,
            recipient_private_key_path=recipient_private_key_path,
        )

        return submission

    def upload(self, s3_options: S3Options) -> str:
        """
        Upload an encrypted submission and return the generated submission ID.

        Verifies that all files were successfully encrypted.
        """
        submission = self.parse_submission()
        encryption_progress_logger = FileProgressLogger[EncryptionState](self.progress_file_encrypt)

        incompletely_encrypted_files = []

        self.__log.info("Verifying encryption status of all submission files…")
        for file_path, file_metadata in submission.files.items():
            state = encryption_progress_logger.get_state(file_path, file_metadata)
            if not state or not state.get("encryption_successful", False):
                incompletely_encrypted_files.append(str(file_path))

        if incompletely_encrypted_files:
            failed_files = "\n - ".join(incompletely_encrypted_files)
            error_msg = (
                "Will not upload, as the following files were not successfully encrypted:\n"
                f"{failed_files}\n"
                "Please re-run the 'encrypt' command and try again."
            )
            self.__log.error(error_msg)
            raise SubmissionValidationError(error_msg)

        self.__log.info("All files verified as successfully encrypted.")

        upload_worker = S3BotoUploadWorker(
            s3_options, status_file_path=self.progress_file_upload, threads=self._threads
        )

        encrypted_submission = self.parse_encrypted_submission()

        upload_worker.upload(encrypted_submission)

        return encrypted_submission.submission_id

    def archive(self, s3_options: S3Options):
        """
        Archive an encrypted submission at a GRZ.
        """
        upload_worker = S3BotoUploadWorker(
            s3_options, status_file_path=self.progress_file_upload, threads=self._threads
        )

        encrypted_submission = self.parse_encrypted_submission()

        upload_worker.archive(encrypted_submission)

    def download(self, s3_options: S3Options, submission_id: str, force: bool = False):
        """
        Download an encrypted submission
        """
        if force:
            # delete the log file if it exists
            self.progress_file_download.unlink(missing_ok=True)

        download_worker = S3BotoDownloadWorker(
            s3_options, status_file_path=self.progress_file_download, threads=self._threads
        )

        self.__log.info("Preparing output directories...")
        download_worker.prepare_download(self.metadata_dir, self.encrypted_files_dir, self.log_dir)

        self.__log.info("Downloading metadata...")
        download_worker.download_metadata(submission_id, self.metadata_dir, metadata_file_name="metadata.json")

        self.__log.info("Downloading encrypted files...")
        download_worker.download(submission_id, EncryptedSubmission(self.metadata_dir, self.encrypted_files_dir))


    def populate(self, configuration,  submission_id: str, force: bool):
        updates = False
        # establish connection to s3; no check because connection established during download
        s3_config = DownloadConfig.model_validate(configuration).s3
        s3_client = init_s3_client(s3_config)

        # get the metadata file object and extract the date
        object_key = f"{submission_id}/metadata/metadata.json"
        response = s3_client.head_object(Bucket= s3_config.bucket, Key= object_key)
        submission_date = response["LastModified"].date()

        # establish connection to database; no check because connection established via DbContext
        db_config = DbConfig.model_validate(configuration).db
        author = _initialise_author_from_dbconfig(db_config)
        db = get_submission_db_instance(db_config.database_url, author=author)
        submission = db.get_submission(submission_id)

        # store entries of donors in database
        donors_in_db_submission = {donor.pseudonym : donor for donor in db.get_donors(submission_id=submission_id)}

        # parse metadata
        metadata_file_path = self.metadata_dir / "metadata.json"
        with open(metadata_file_path, encoding="utf-8") as metadata_file:
            metadata = GrzSubmissionMetadata.model_validate_json(metadata_file.read())

            # metadata_file.seek(0)
            # metadata_content = json.load(metadata_file)
            # redact tanG, local case id and potential tanG in the patient information
            # metadata_content = redact_metadata(metadata_content)
            # metadata_string = json.dumps(metadata_content)

            metadata_content = metadata.to_redacted_dict(False)
            metadata_string = json.dumps(metadata_content)
            print(metadata_string)
            exit()

        # add the filesize; would be better with SubmissionMetadata but then without validation
        # submission_size: int = sum([i.file_size_in_bytes for i in submission_metadata.files.values()])
        submission_size: int = 0
        for donor in metadata.donors:
            for lab_datum in donor.lab_data:
                if sequence_data := lab_datum.sequence_data:
                    for file in sequence_data.files:
                        submission_size += file.file_size_in_bytes

        # populate submission metadata
        submission_information = diff_metadata(submission, metadata, submission_size, submission_date, metadata_string)
        if submission_information.change: updates = True

        donors_in_metadata: list[str] = []
        donor_diff = DonorDiff([], [], [], [], [])
        for donor in metadata.donors:
            donor_data = diff_donor(donor, donors_in_db_submission, submission_id)
            donors_in_metadata.append(donor_data.name)
            if donor_data.status != "unchanged":
                updates = True
                if donor_data.status == "new":
                    donor_diff.added.append(donor_data.donor)
                else:
                    donor_diff.updated.append(donor_data.donor)
            else:
                donor_diff.unchanged.append(donor_data.donor)

        donor_diff.deleted = [db_donor for db_donor in donors_in_db_submission.values() if db_donor.pseudonym not in donors_in_metadata]

        if updates and not force:
            raise RuntimeError(f"Changes in {object_key} detected. Re-run with --populate to commit them.")

        if updates:
            # goes below
            if submission_information.update:
                self.__log.info(f"Submission: {submission_id} - Updating fields: {', '.join(submission_information.update)} in database")
            if submission_information.unchanged:
                self.__log.info(f"Submission: {submission_id} - Not updating fields: {', '.join(submission_information.unchanged)} in database")
            if donor_diff.unchanged:
                self.__log.info(f"Submission: {submission_id} - Keep existing donor(s) in database: {','.join([i.pseudonym for i in donor_diff.unchanged])}")
            if donor_diff.added:
                self.__log.info(f"Submission: {submission_id} - Adding new donor(s) to database: {','.join([i.pseudonym for i in donor_diff.added])} from metadata.json")
            if donor_diff.updated:
                self.__log.info(f"Submission: {submission_id} - Modify existing donor(s) in database: {','.join([i.pseudonym for i in donor_diff.updated])} from metadata.json")
            if donor_diff.deleted:
                self.__log.info(f"Submission: {submission_id} - Drop existing donor(s) in db: {','.join([i.pseudonym for i in donor_diff.deleted])}")

            # push to database
            for key, unused, value in submission_information.db_ingest:
                db.modify_submission(submission_id, key, value)
            for donor in donor_diff.added:
                db.add_donor(donor)
            for donor in donor_diff.updated:
                db.update_donor(donor)
            for donor in donor_diff.deleted:
                db.delete_donor(donor)
        else:
            self.__log.info(f"Submission: {submission_id} - No updates from metadata.json necessary")


