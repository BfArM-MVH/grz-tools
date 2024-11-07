"""
Module handling the submission processing:
- Validation
- Encryption
- Decryption
"""

import logging
from os import PathLike
from pathlib import Path

from .download import S3BotoDownloadWorker
from .file_operations import Crypt4GH
from .models.config import Backend, ConfigModel
from .operations.encrypt import _parallel_encrypt
from .progress_logging import FileProgressLogger
from .submission import (
    EncryptedSubmission,
    Submission,
    SubmissionValidationError,
    log,
)
from .upload import S3BotoUploadWorker


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
        Initialize the operations object.

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
        if force:
            # delete the log file
            self.progress_file_encrypt.unlink()

        encrypted_submission = self._encrypt(
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

    def _encrypt(
        self,
        recipient_public_key_path: str | PathLike,
        submitter_private_key_path: str | PathLike | None = None,
        logger: logging.Logger = logging.getLogger(__name__),
    ) -> EncryptedSubmission:
        """
        Encrypt a submission with a public key using Crypt4Gh

        :param recipient_public_key_path: Path to the public key file which will be used for encryption
        :param submitter_private_key_path: Path to the private key file which will be used to sign the encryption
        :return: EncryptedSubmission instance
        """
        submission = self.parse_submission()

        encrypted_files_dir = Path(self.encrypted_files_dir)

        if not Path(recipient_public_key_path).expanduser().is_file():
            msg = f"Public key file does not exist: {recipient_public_key_path}"
            logger.error(msg)
            raise FileNotFoundError(msg)
        if not submitter_private_key_path:
            logger.warning("No submitter private key provided, skipping signing.")
        elif not Path(submitter_private_key_path).expanduser().is_file():
            msg = f"Private key file does not exist: {submitter_private_key_path}"
            logger.error(msg)
            raise FileNotFoundError(msg)

        if not encrypted_files_dir.is_dir():
            logger.debug(
                "Creating encrypted submission files directory: %s...",
                encrypted_files_dir,
            )
            encrypted_files_dir.mkdir(mode=0o770, parents=False, exist_ok=False)

        progress_logger = FileProgressLogger(log_file_path=self.progress_file_encrypt)

        try:
            public_keys = Crypt4GH.prepare_c4gh_keys(recipient_public_key_path)
        except Exception as e:
            logger.error(f"Error preparing public keys: {e}")
            raise e

        _parallel_encrypt(submission, public_keys, encrypted_files_dir, progress_logger)

        logger.info("File encryption completed.")

        return EncryptedSubmission(
            metadata_dir=submission.metadata_dir,
            encrypted_files_dir=encrypted_files_dir,
        )
