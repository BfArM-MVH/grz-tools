"""Worker class for handling submission processing."""

from __future__ import annotations

import logging
from os import PathLike
from pathlib import Path

from ..models.identifiers import IdentifiersModel
from ..models.s3 import S3Options
from .download import S3BotoDownloadWorker
from .submission import EncryptedSubmission, Submission, SubmissionValidationError

log = logging.getLogger(__name__)


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

    def validate(self, identifiers: IdentifiersModel, force=False):
        """
        Validate this submission

        :param identifiers: IdentifiersModel containing GRZ and LE identifiers.
        :param force: Force validation of already validated files
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
        force: bool = False,
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
        :param recipient_public_key_path: Path to the public key file of the recipient.
        :param submitter_private_key_path: Path to the private key file of the submitter.
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
        Upload an encrypted submission and return the generated submission ID
        """
        from .upload import S3BotoUploadWorker

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
        from .upload import S3BotoUploadWorker

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
