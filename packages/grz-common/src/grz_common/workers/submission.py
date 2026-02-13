"""Classes for parsing and validating submission metadata and files."""

from __future__ import annotations

import hashlib
import json
import logging
import os
import shutil
import subprocess
import typing
from collections.abc import Generator
from itertools import groupby
from os import PathLike
from pathlib import Path
from typing import Any

from grz_pydantic_models.submission.metadata import get_accepted_versions
from grz_pydantic_models.submission.metadata.v1 import (
    Donor,
    FileType,
    GrzSubmissionMetadata,
    LabDatum,
    ReadOrder,
    SequencingLayout,
)
from grz_pydantic_models.submission.metadata.v1 import File as SubmissionFileMetadata
from grz_pydantic_models.submission.thresholds import Thresholds
from pydantic import ValidationError
from tqdm.auto import tqdm

from ..constants import TQDM_DEFAULTS
from ..models.identifiers import IdentifiersModel
from ..pipeline.components import TqdmObserver
from ..pipeline.components.crypt4gh import Crypt4GHDecryptor, Crypt4GHEncryptor
from ..pipeline.components.validation import ValidatorObserver
from ..progress import DecryptionState, EncryptionState, FileProgressLogger, ValidationState
from ..utils.crypt import Crypt4GH
from ..validation import UserInterruptException, run_grz_check

log = logging.getLogger(__name__)

S3_MAX_KEY_LENGTH = 1024
# length of uploaded prefix before LE-specified filepath
# e.g. len("123456789_2025-04-01_a70eb6ce/files/") == 36
UPLOADED_FILE_PREFIX_LENGTH = 36


class SubmissionMetadata:
    """Class for reading and validating submission metadata"""

    __log = log.getChild("SubmissionMetadata")

    def __init__(self, metadata_file: Path):
        """
        Load, parse and validate the metadata file.

        :param metadata_file: path to the metadata.json file
        :raises json.JSONDecodeError: if failed to read the metadata.json file
        :raises jsonschema.exceptions.ValidationError: if metadata does not match expected schema
        """
        self.file_path = metadata_file
        self.content = self._read_metadata(self.file_path)
        self._checksum = self._calculate_metadata_checksum(self.file_path)

        self._files: dict | None = None

    def _calculate_metadata_checksum(self, file_path: Path) -> str:
        """Calculate SHA256 checksum of the metadata file."""
        return hashlib.sha256(open(file_path, "rb").read(), usedforsecurity=False).hexdigest()

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
                if not lab_data.sequence_data:
                    continue
                for file_data in lab_data.sequence_data.files:
                    file_path = Path(file_data.file_path)
                    submission_files[file_path] = file_data

        self._files = submission_files
        return self._files

    @staticmethod
    def pair_files(
        files: list[SubmissionFileMetadata],
    ) -> Generator[tuple[SubmissionFileMetadata, SubmissionFileMetadata], None, None]:
        """
        Groups a list of FASTQ files by Flowcell and Lane, yielding (R1, R2) pairs.
        """
        key_func = lambda f: (f.flowcell_id, f.lane_id)
        files.sort(key=key_func)

        for _, group in groupby(files, key_func):
            group_files = list(group)

            r1_files = [f for f in group_files if f.read_order == ReadOrder.r1]
            r2_files = [f for f in group_files if f.read_order == ReadOrder.r2]

            if r1_files and r2_files:
                yield from zip(r1_files, r2_files, strict=True)

    def iter_paired_end_fastqs(
        self,
    ) -> Generator[
        tuple[Donor, LabDatum, list[tuple[SubmissionFileMetadata, SubmissionFileMetadata]], Thresholds], None, None
    ]:
        """
        Yields (Donor, LabDatum, List of R1/R2 Pairs, Threshold) for every Paired-End unit.
        """
        for donor in self.content.donors:
            for lab_datum in donor.lab_data:
                if lab_datum.sequencing_layout != SequencingLayout.paired_end:
                    continue
                if not lab_datum.sequence_data:
                    continue

                files = [f for f in lab_datum.sequence_data.files if f.file_type == FileType.fastq]
                thresholds = self.content.determine_thresholds_for(donor, lab_datum)
                pairs = list(self.pair_files(files))

                if pairs:
                    yield donor, lab_datum, pairs, thresholds

    def iter_single_end_fastqs(
        self,
    ) -> Generator[tuple[Donor, LabDatum, list[SubmissionFileMetadata], Thresholds], None, None]:
        """
        Yields (Donor, LabDatum, List of Files, Threshold) for every Single-End unit.
        """
        for donor in self.content.donors:
            for lab_datum in donor.lab_data:
                # Handle Single, Reverse, Other, etc.
                if lab_datum.sequencing_layout == SequencingLayout.paired_end:
                    continue
                if not lab_datum.sequence_data:
                    continue

                files = [f for f in lab_datum.sequence_data.files if f.file_type == FileType.fastq]
                thresholds = self.content.determine_thresholds_for(donor, lab_datum)

                if files:
                    yield donor, lab_datum, files, thresholds

    def iter_bams(self) -> Generator[tuple[Donor, LabDatum, SubmissionFileMetadata], None, None]:
        """
        Yields every BAM file in the submission.
        """
        for donor in self.content.donors:
            for lab_datum in donor.lab_data:
                if not lab_datum.sequence_data:
                    continue
                for f in lab_datum.sequence_data.files:
                    if f.file_type == FileType.bam:
                        yield donor, lab_datum, f

    def validate(self, identifiers: IdentifiersModel) -> Generator[str]:  # noqa: C901, PLR0912
        """
        Validates this submission's metadata (content).

        :return: Generator of errors
        """
        metadata_schema_version = self.content.get_schema_version()
        accepted_versions = get_accepted_versions()
        if metadata_schema_version not in accepted_versions:
            yield f"Metadata schema version {metadata_schema_version} is outdated. Currently accepting the following versions: {', '.join(accepted_versions)}"

        expected_grz_id, expected_le_id = identifiers.grz, identifiers.le
        if (submitted_grz_id := self.content.submission.genomic_data_center_id) != expected_grz_id:
            yield (
                f"Genomic data center identifier specified in the metadata.json ({submitted_grz_id}) "
                f"does not match genomic data center identifier in config ({expected_grz_id})"
            )

        if (submitted_le_id := self.content.submission.submitter_id) != expected_le_id:
            yield (
                f"Submitter (LE) identifier specified in the metadata.json ({submitted_le_id}) "
                f"does not match submitter (LE) identifier in config ({expected_le_id})"
            )

        submission_files: dict[str | PathLike, SubmissionFileMetadata] = {}
        for donor in self.content.donors:
            for lab_data in donor.lab_data:
                if not lab_data.sequence_data:
                    log.info(f"Skipping validation of empty sequence data for donor {donor}.")
                    continue
                for file_data in lab_data.sequence_data.files:
                    # check if file is already registered
                    file_path = Path(file_data.file_path)
                    if len(str(file_path)) > (S3_MAX_KEY_LENGTH - UPLOADED_FILE_PREFIX_LENGTH):
                        yield f"{file_data.file_path}: File path is too long for the inbox!"
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

    def validate_files_with_grz_check(  # noqa: C901, PLR0912
        self, checksum_progress_file: str | PathLike, seq_data_progress_file: str | PathLike, threads: int | None
    ) -> Generator[str, None, None]:
        """
        Validates submission files using `grz-check` and populates both progress logs.
        """
        checksum_progress_logger = FileProgressLogger[ValidationState](log_file_path=checksum_progress_file)
        checksum_progress_logger.cleanup(keep=[(fp, fm) for fp, fm in self.files.items()])

        seq_data_progress_logger = FileProgressLogger[ValidationState](log_file_path=seq_data_progress_file)
        seq_data_progress_logger.cleanup(keep=[(fp, fm) for fp, fm in self.files.items()])

        log_dir = Path(checksum_progress_file).parent

        grz_check_args = []
        checked_files = set()

        def should_check_file(file_path: Path, file_metadata: SubmissionFileMetadata) -> bool:
            # Check against both logs. If either is missing a "pass", re-check.
            checksum_state = checksum_progress_logger.get_state(file_path, file_metadata)
            seq_data_state = seq_data_progress_logger.get_state(file_path, file_metadata)
            is_seq_file = file_metadata.file_type in ("fastq", "bam")

            checksum_passed = checksum_state and checksum_state.get("validation_passed")
            seq_data_passed = seq_data_state and seq_data_state.get("validation_passed")

            if is_seq_file:
                return not (checksum_passed and seq_data_passed)
            return not checksum_passed

        for _donor, _lab_datum, pairs, thresholds in self.metadata.iter_paired_end_fastqs():
            mean_read_length_threshold = thresholds.mean_read_length
            for r1_meta, r2_meta in pairs:
                r1_path = self.files_dir / r1_meta.file_path
                r2_path = self.files_dir / r2_meta.file_path

                if should_check_file(r1_path, r1_meta) or should_check_file(r2_path, r2_meta):
                    grz_check_args.extend(
                        ["--fastq-paired", str(r1_path), str(r2_path), str(mean_read_length_threshold)]
                    )
                checked_files.update({r1_path, r2_path})

        for _donor, _lab_datum, fastq_files, thresholds in self.metadata.iter_single_end_fastqs():
            mean_read_length_threshold = thresholds.mean_read_length
            for f_meta in fastq_files:
                f_path = self.files_dir / f_meta.file_path
                if f_path not in checked_files:
                    if should_check_file(f_path, f_meta):
                        grz_check_args.extend(["--fastq-single", str(f_path), str(mean_read_length_threshold)])
                    checked_files.add(f_path)

        for _donor, _lab_datum, bam_meta in self.metadata.iter_bams():
            bam_path = self.files_dir / bam_meta.file_path
            if bam_path not in checked_files:
                if should_check_file(bam_path, bam_meta):
                    grz_check_args.extend(["--bam", str(bam_path)])
                checked_files.add(bam_path)

        # Handle any other files with --raw for calculating checksums
        for file_path, file_metadata in self.files.items():
            if file_path not in checked_files and should_check_file(file_path, file_metadata):
                grz_check_args.extend(["--raw", str(file_path)])

        if not grz_check_args:
            self.__log.info("All files are already validated. Skipping `grz-check`.")
        else:
            yield from self._run_grz_check_command(
                grz_check_args, threads, log_dir, checksum_progress_logger, seq_data_progress_logger
            )

        yield from self._aggregate_validation_errors(checksum_progress_logger, seq_data_progress_logger)

    def _run_grz_check_command(
        self,
        grz_check_args: list[str],
        threads: int | None,
        log_dir: Path,
        checksum_logger: FileProgressLogger,
        seq_logger: FileProgressLogger,
    ) -> Generator[str, None, None]:
        """Helper to encapsulate the subprocess execution and report processing."""
        temp_report_path = log_dir / "grz-check.report.jsonl"
        temp_report_path.unlink(missing_ok=True)

        command_args = ["--output", str(temp_report_path), *grz_check_args]
        if threads:
            command_args.extend(["--threads", str(threads)])

        try:
            run_grz_check(command_args)
        except UserInterruptException:
            self.__log.warning("Validation cancelled by user. Processing partial results...")
            raise
        except subprocess.CalledProcessError as e:
            self.__log.error(f"`grz-check` failed with exit code {e.returncode}")
            yield "`grz-check` execution failed. See logs for details."
        finally:
            if temp_report_path.is_file():
                with temp_report_path.open("r") as f:
                    self._process_grz_check_report(f, checksum_logger, seq_logger)
                temp_report_path.unlink()

    def _aggregate_validation_errors(
        self, checksum_progress_logger: FileProgressLogger, seq_data_progress_logger: FileProgressLogger
    ) -> Generator[str, None, None]:
        """Aggregates all errors from both progress loggers into a flat generator."""
        all_errors = set()
        for local_file_path, file_metadata in self.files.items():
            checksum_state = checksum_progress_logger.get_state(local_file_path, file_metadata)
            if checksum_state and not checksum_state.get("validation_passed"):
                for error in checksum_state.get("errors", []):
                    all_errors.add(f"{local_file_path.relative_to(self.files_dir)}: {error}")

            if file_metadata.file_type in ("fastq", "bam"):
                seq_data_state = seq_data_progress_logger.get_state(local_file_path, file_metadata)
                if seq_data_state and not seq_data_state.get("validation_passed"):
                    for error in seq_data_state.get("errors", []):
                        all_errors.add(f"{local_file_path.relative_to(self.files_dir)}: {error}")
        yield from all_errors

    def _process_grz_check_report(  # noqa: C901
        self,
        report_file: typing.TextIO,
        checksum_progress_logger: FileProgressLogger[ValidationState],
        seq_data_progress_logger: FileProgressLogger[ValidationState],
    ):
        """
        Parses the JSONL report from `grz-check` and updates both progress loggers.
        """
        for line in report_file:
            try:
                report_entry = json.loads(line)
                data = report_entry.get("data", {})
                file_path_str = data.get("path")
                if not file_path_str:
                    continue

                file_path = Path(file_path_str).resolve()
                file_metadata = self.files.get(file_path)

                if not file_metadata:
                    self.__log.warning(f"Could not find metadata for file in grz-check report: {file_path_str}")
                    continue

                status = data.get("status")
                errors = data.get("errors", [])
                warnings = data.get("warnings", [])
                checksum = data.get("checksum")

                if warnings:
                    for w in warnings:
                        self.__log.warning(w)

                checksum_issues = []
                if (
                    checksum
                    and (file_metadata.checksum_type or "").lower() == "sha256"
                    and file_metadata.file_checksum != checksum
                ):
                    checksum_issues.append(
                        f"Checksum mismatch! Expected: '{file_metadata.file_checksum}', calculated: '{checksum}'"
                    )

                if file_path.exists() and file_path.is_file():
                    if file_metadata.file_size_in_bytes != file_path.stat().st_size:
                        checksum_issues.append(
                            f"File size mismatch! Expected: '{file_metadata.file_size_in_bytes}', observed: '{file_path.stat().st_size}'."
                        )
                else:
                    checksum_issues.append("File not found for size check.")

                checksum_passed = not checksum_issues
                checksum_state = ValidationState(errors=checksum_issues, validation_passed=checksum_passed)
                checksum_progress_logger.set_state(file_path, file_metadata, checksum_state)

                if file_metadata.file_type in ("fastq", "bam"):
                    integrity_passed = status == "OK"
                    integrity_errors = errors if not integrity_passed else []

                    seq_data_state = ValidationState(errors=integrity_errors, validation_passed=integrity_passed)
                    seq_data_progress_logger.set_state(file_path, file_metadata, seq_data_state)

            except json.JSONDecodeError:
                self.__log.warning(f"Could not parse line in grz-check report: {line.strip()}")
            except Exception as e:
                self.__log.error(f"Error processing grz-check report entry: {line.strip()}. Error: {e}")

    def validate_files_fallback(
        self, checksum_progress_file: str | PathLike, seq_data_progress_file: str | PathLike
    ) -> Generator[str, None, None]:
        """
        Validates files (checksum and content) using streaming validation.
        Replaces `grz-check` usage when external tools are unavailable.
        """
        checksum_logger = FileProgressLogger[ValidationState](log_file_path=checksum_progress_file)
        checksum_logger.cleanup(keep=[(fp, fm) for fp, fm in self.files.items()])

        sequence_logger = FileProgressLogger[ValidationState](log_file_path=seq_data_progress_file)
        sequence_logger.cleanup(keep=[(fp, fm) for fp, fm in self.files.items()])

        # get thresholds for fastq validation
        threshold_map = {}
        for _, _, pairs, thresholds in self.metadata.iter_paired_end_fastqs():
            for f1, f2 in pairs:
                threshold_map[f1.file_path] = thresholds.mean_read_length
                threshold_map[f2.file_path] = thresholds.mean_read_length
        for _, _, files, thresholds in self.metadata.iter_single_end_fastqs():
            for f in files:
                threshold_map[f.file_path] = thresholds.mean_read_length

        # validate files, including fastq read counts
        read_counts = {}
        self._validate_files_fallback(checksum_logger, read_counts, sequence_logger, threshold_map)
        self._validate_paired_end_readcounts(read_counts, sequence_logger)

        yield from self._aggregate_validation_errors(checksum_logger, sequence_logger)

    def _validate_paired_end_readcounts(
        self, read_counts: dict[Any, Any], sequence_logger: FileProgressLogger[ValidationState]
    ):
        """
        Ensure paired end read counts match
        """
        for _, _, pairs, _ in self.metadata.iter_paired_end_fastqs():
            for f1, f2 in pairs:
                c1 = read_counts.get(f1.file_path)
                c2 = read_counts.get(f2.file_path)

                if c1 is not None and c2 is not None and c1 != c2:
                    msg = f"Read count mismatch: {f1.file_path}({c1}) != {f2.file_path}({c2})"

                    for meta in (f1, f2):
                        path = self.files_dir / meta.file_path
                        state = sequence_logger.get_state(path, meta) or ValidationState(
                            errors=[], validation_passed=False
                        )
                        if msg not in state.get("errors", []):
                            state["errors"].append(msg)
                            state["validation_passed"] = False
                            sequence_logger.set_state(path, meta, state)

    def _validate_files_fallback(
        self,
        checksum_logger: FileProgressLogger[ValidationState],
        read_counts: dict[Any, Any],
        sequence_logger: FileProgressLogger[ValidationState],
        threshold_map: dict[Any, Any],
    ):
        for local_path, meta in self.files.items():
            # check existing state from progress logs
            checksum_state = checksum_logger.get_state(local_path, meta)
            sequence_state = sequence_logger.get_state(local_path, meta)

            checksum_passed = checksum_state and checksum_state.get("validation_passed")
            sequence_passed = sequence_state and sequence_state.get("validation_passed")
            is_fastq_or_bam = meta.file_type in {FileType.fastq, FileType.bam}

            # if it's not a fastq/bam file and checksums are fine, skip
            if checksum_passed and (not is_fastq_or_bam or sequence_passed):
                continue

            self.__log.debug("Validating '%s'...", str(local_path))
            errors = []

            try:
                if not local_path.is_file():
                    errors.append(f"File not found: {local_path}")
                else:
                    with open(local_path, "rb") as f:
                        validator = ValidatorObserver(
                            f,
                            file_type=meta.file_type,
                            expected_checksum=meta.file_checksum,
                            mean_read_length_threshold=threshold_map.get(meta.file_path),
                        )
                        while validator.read(1024**2):
                            pass

                        validator.validate()
                        if "read_count" in validator.metrics:
                            read_counts[meta.file_path] = validator.metrics["read_count"]

            except Exception as e:
                errors.append(str(e))

            # update progress logs
            state = ValidationState(errors=errors, validation_passed=not errors)
            checksum_logger.set_state(local_path, meta, state)
            if is_fastq_or_bam:
                sequence_logger.set_state(local_path, meta, state)

    def encrypt(
        self,
        encrypted_files_dir: str | PathLike,
        progress_log_file: str | PathLike,
        recipient_public_key_path: str | PathLike,
        submitter_private_key_path: str | PathLike | None = None,
        force: bool = False,
    ) -> EncryptedSubmission:
        """
        Encrypt this submission with a public key using Crypt4Gh

        :param encrypted_files_dir: Output directory of the encrypted files
        :param progress_log_file: Path to a log file to store the progress of the encryption process
        :param recipient_public_key_path: Path to the public key file which will be used for encryption
        :param submitter_private_key_path: Path to the private key file which will be used to sign the encryption
        :param force: Force encryption even if target files already exist
        :return: EncryptedSubmission instance
        """
        # Import here to avoid circular import issues
        from ..progress import FileProgressLogger  # noqa: PLC0415

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
        # TODO parse and use submitter_private_key_path if available

        if not encrypted_files_dir.is_dir():
            self.__log.debug(
                "Creating encrypted submission files directory: %s...",
                encrypted_files_dir,
            )
            encrypted_files_dir.mkdir(mode=0o770, parents=False, exist_ok=False)

        progress_logger = FileProgressLogger[EncryptionState](log_file_path=progress_log_file)

        try:
            public_keys = Crypt4GH.prepare_c4gh_keys(recipient_public_key_path)
            recipient_public_key = public_keys[0][2]
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

                if encrypted_file_path.exists() and not force:
                    raise RuntimeError(
                        f"'{encrypted_file_path}' already exists. Delete it or use --force to overwrite it."
                    )

                try:
                    with (
                        open(file_path, "rb") as src,
                        open(encrypted_file_path, "wb") as f,
                        tqdm(
                            total=os.stat(file_path).st_size,
                            desc="ENCRYPT ",
                            postfix={"file": Path(file_path).name},
                            leave=False,
                            **TQDM_DEFAULTS,
                        ) as pbar,
                        TqdmObserver(src, pbar=pbar) as monitored,
                        Crypt4GHEncryptor(monitored, recipient_public_key, None) as encryptor,
                    ):
                        shutil.copyfileobj(encryptor, f)

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

    def __init__(
        self, metadata_dir: str | PathLike, encrypted_files_dir: str | PathLike, log_dir: str | PathLike | None = None
    ):
        """
        Initialize the encrypted submission object.

        :param metadata_dir: Path to the metadata directory
        :param encrypted_files_dir: Path to the encrypted files directory
        """
        self.metadata_dir = Path(metadata_dir)
        self.encrypted_files_dir = Path(encrypted_files_dir)
        self.log_dir = Path(log_dir) if log_dir is not None else None

        self.metadata = SubmissionMetadata(self.metadata_dir / "metadata.json")

    @property
    def encrypted_files(self) -> dict[Path, SubmissionFileMetadata]:
        """
        The encrypted files linked in the metadata.

        :return: Dictionary of `local_file_path` -> `SubmissionFileMetadata` pairs.
        """
        retval = {}
        for file_path, file_metadata in self.metadata.files.items():
            encrypted_file_path = self.get_encrypted_file_path(self.encrypted_files_dir / file_path)

            retval[encrypted_file_path] = file_metadata

        return retval

    @property
    def submission_id(self) -> str:
        return self.metadata.content.submission_id

    def get_metadata_file_path_and_object_id(self) -> tuple[Path, str]:
        """
        :return: tuple with the `local_file_path` and s3_object_id of the metadata file
        """
        return Path(self.metadata.file_path), str(Path(self.submission_id) / "metadata" / self.metadata.file_path.name)

    def get_encrypted_files_and_object_id(self) -> dict[Path, str]:
        """
        :return Dictionary of `local_file_path` -> s3_object_id
        """
        retval = {}
        for local_file_path, file_metadata in self.encrypted_files.items():
            retval[local_file_path] = str(
                Path(self.submission_id) / "files" / self.get_encrypted_file_path(file_metadata.file_path)
            )
        return retval

    def get_log_files_and_object_id(self) -> dict[Path, str]:
        """
        :return Dictionary of `local_file_path` -> s3_object_id
        """
        retval = {}
        if self.log_dir is not None:
            log_dir = self.log_dir
            for dirpath, _dirnames, filenames in log_dir.walk():
                for filename in filenames:
                    local_file_path = dirpath / filename
                    retval[local_file_path] = str(
                        Path(self.submission_id) / "logs" / local_file_path.relative_to(log_dir)
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
        # Import here to avoid circular import issues
        from ..progress import FileProgressLogger  # noqa: PLC0415

        files_dir = Path(files_dir)

        if not files_dir.is_dir():
            self.__log.debug(
                "Creating decrypted submission files directory: %s...",
                files_dir,
            )
            files_dir.mkdir(mode=0o770, parents=False, exist_ok=False)

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
                    with (
                        open(encrypted_file_path, "rb") as src,
                        open(decrypted_file_path, "wb") as f,
                        tqdm(
                            total=os.stat(encrypted_file_path).st_size,
                            desc="DECRYPT ",
                            postfix={"file": Path(encrypted_file_path).name},
                            leave=False,
                            **TQDM_DEFAULTS,
                        ) as pbar,
                        TqdmObserver(src, pbar=pbar) as monitored,
                        Crypt4GHDecryptor(monitored, private_key) as decryptor,
                    ):
                        shutil.copyfileobj(decryptor, f)

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
