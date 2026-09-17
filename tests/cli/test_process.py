"""
Integration tests for grzctl process command.

This tests the full streaming pipeline with mocked S3 buckets (inbox + archive).
"""

import gzip
import hashlib
import json
import os
from collections import Counter
from dataclasses import dataclass, field
from io import BytesIO
from pathlib import Path

import boto3
import botocore.client
import botocore.exceptions
import click.testing
import crypt4gh.keys
import crypt4gh.lib
import grzctl.cli
import pytest
import yaml
from grz_db.models.submission import FailureReasonEnum, SubmissionDb, SubmissionStateEnum

# Path to test fixtures
MOCK_FILES_DIR = Path(__file__).parent.parent / "mock_files"
VALID_SUBMISSION_DIR = MOCK_FILES_DIR / "submissions" / "valid_submission"


@pytest.fixture
def process_config_content(
    crypt4gh_grz_private_key_file_path,
    crypt4gh_grz_public_key_file_path,
    migrated_db_config_content,
    tmp_path,
):
    db_config = migrated_db_config_content["db"]
    db_config["author"]["private_key_passphrase"] = "test"
    local_storage_dir = tmp_path / "local_storage"
    local_storage_dir.mkdir()

    return {
        "leistungserbringer": {
            "260914050": {
                "inbox_buckets": {
                    "inbox": {
                        "endpoint_url": "https://s3.amazonaws.com",
                        "access_key": "testing",
                        "secret": "testing",
                        "private_key_path": str(crypt4gh_grz_private_key_file_path),
                    }
                }
            }
        },
        "archives": {
            "consented": {
                "s3": {
                    "endpoint_url": "https://s3.amazonaws.com",
                    "bucket": "consented-archive",
                    "access_key": "testing",
                    "secret": "testing",
                },
                "public_key_path": str(crypt4gh_grz_public_key_file_path),
            },
            "non_consented": {
                "s3": {
                    "endpoint_url": "https://s3.amazonaws.com",
                    "bucket": "non-consented-archive",
                    "access_key": "testing",
                    "secret": "testing",
                },
                "public_key_path": str(crypt4gh_grz_public_key_file_path),
            },
            "interrogation": {
                "s3": {
                    "endpoint_url": "https://s3.amazonaws.com",
                    "bucket": "interrogation-archive",
                    "access_key": "testing",
                    "secret": "testing",
                },
                "keep_failed": False,
            },
        },
        "pruefbericht": {
            "authorization_url": "https://bfarm.localhost/token",
            "api_base_url": "https://bfarm.localhost/api/",
            "client_id": "pytest",
            "client_secret": "pysecret",
        },
        "db": db_config,
        "detailed_qc": {
            "salt": "salty",
            "target_percentage": "0.0",
            "local_storage": str(local_storage_dir),
            "auto_run": False,
        },
        "keys": {
            "grz_private_key_path": str(crypt4gh_grz_private_key_file_path),
            "grz_public_key_path": str(crypt4gh_grz_public_key_file_path),
        },
        "identifiers": {
            "grz": "GRZT00000",
        },
    }


@pytest.fixture
def temp_process_config_file_path(tmp_path, process_config_content) -> Path:
    """Write the process config to a YAML file."""
    config_file = tmp_path / "config.process.yaml"
    with open(config_file, "w") as fd:
        yaml.dump(process_config_content, fd)
    return config_file


@pytest.fixture
def s3_buckets(aws_credentials):
    """Create the four S3 buckets needed for processing."""
    conn = boto3.client("s3")

    # Create inbox, consented-archive, non-consented-archive and interrogation-archive buckets
    conn.create_bucket(Bucket="inbox")
    conn.create_bucket(Bucket="consented-archive")
    conn.create_bucket(Bucket="non-consented-archive")
    conn.create_bucket(Bucket="interrogation-archive")

    s3 = boto3.resource("s3")
    return {
        "inbox": s3.Bucket("inbox"),
        "consented": s3.Bucket("consented-archive"),
        "non_consented": s3.Bucket("non-consented-archive"),
        "interrogation": s3.Bucket("interrogation-archive"),
    }


def upload_submission_to_inbox(inbox_bucket, submission_id: str):
    """
    Upload the valid submission's encrypted files and metadata to the inbox bucket.

    The inbox structure is:
        {submission_id}/
            metadata/metadata.json
            files/*.c4gh
    """
    # Upload metadata
    metadata_path = VALID_SUBMISSION_DIR / "metadata" / "metadata.json"
    inbox_bucket.upload_file(
        Filename=str(metadata_path),
        Key=f"{submission_id}/metadata/metadata.json",
    )

    # Upload encrypted files
    encrypted_files_dir = VALID_SUBMISSION_DIR / "encrypted_files"
    for encrypted_file in encrypted_files_dir.glob("*.c4gh"):
        inbox_bucket.upload_file(
            Filename=str(encrypted_file),
            Key=f"{submission_id}/files/{encrypted_file.name}",
        )


def _run_process(config_file_path: Path, submission_id: str, output_dir: Path, *extra_args: str):
    args = [
        "--config",
        str(config_file_path),
        "process",
        "--submission-id",
        submission_id,
        "--output-dir",
        str(output_dir),
        "--no-submit-pruefbericht",
        *extra_args,
    ]
    return click.testing.CliRunner().invoke(grzctl.cli.build_cli(), args)


class TestGrzctlProcess:
    """Integration tests for grzctl process command."""

    def test_process_submission_basic(
        self,
        s3_buckets,
        temp_process_config_file_path,
        crypt4gh_grz_private_key_file_path,
        working_dir_path,
    ):
        """
        Test basic streaming pipeline processing.

        This test:
        1. Uploads a submission to the mock inbox bucket
        2. Runs grzctl process
        3. Verifies files are uploaded to the archive bucket
        """
        submission_id = "260914050_2024-07-15_c64603a7"

        # Upload submission to inbox
        upload_submission_to_inbox(s3_buckets["inbox"], submission_id)

        # Verify inbox has files
        inbox_keys = {o.key for o in s3_buckets["inbox"].objects.all()}
        assert f"{submission_id}/metadata/metadata.json" in inbox_keys
        assert any("files/" in key and ".c4gh" in key for key in inbox_keys)

        # Run grzctl process
        result = _run_process(temp_process_config_file_path, submission_id, working_dir_path, "--no-update-db")

        assert result.exit_code == 0, f"Process failed: {result.output}"

        # Check that files were uploaded to the consented archive
        # (the valid_submission has research consent)
        consented_keys = {o.key for o in s3_buckets["consented"].objects.all()}
        assert f"{submission_id}/metadata/metadata.json" in consented_keys
        _assert_archived(s3_buckets["consented"], submission_id, crypt4gh_grz_private_key_file_path)

    def test_process_archives_redacted_logs(
        self,
        s3_buckets,
        temp_process_config_file_path,
        working_dir_path,
    ):
        submission_id = "260914050_2024-07-15_c64603a7"
        upload_submission_to_inbox(s3_buckets["inbox"], submission_id)

        metadata_path = VALID_SUBMISSION_DIR / "metadata" / "metadata.json"
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        tan = metadata["submission"]["tanG"]
        local_case_id = metadata["submission"]["localCaseId"]

        logs_dir = working_dir_path / "logs"
        logs_dir.mkdir(parents=True, exist_ok=True)
        custom_log = logs_dir / "custom.log"
        custom_log.write_text(f"tan={tan}\nlocalCaseId={local_case_id}\n", encoding="utf-8")

        result = _run_process(temp_process_config_file_path, submission_id, working_dir_path, "--no-update-db")
        assert result.exit_code == 0, f"Process failed: {result.output}"

        archived_key = f"{submission_id}/logs/custom.log"
        archived_keys = {o.key for o in s3_buckets["consented"].objects.all()}
        assert archived_key in archived_keys

        archived_log = s3_buckets["consented"].Object(archived_key).get()["Body"].read().decode("utf-8")
        assert tan not in archived_log
        assert local_case_id not in archived_log
        assert "REDACTED_TAN_G" in archived_log
        assert "REDACTED_LOCAL_CASE_ID" in archived_log

    def test_process_creates_correct_directory_structure(
        self,
        s3_buckets,
        temp_process_config_file_path,
        working_dir_path,
    ):
        """
        Test that grzctl process creates the expected local directory structure.
        """
        submission_id = "260914050_2024-07-15_c64603a7"

        # Upload submission to inbox
        upload_submission_to_inbox(s3_buckets["inbox"], submission_id)

        # Run grzctl process
        result = _run_process(temp_process_config_file_path, submission_id, working_dir_path, "--no-update-db")

        assert result.exit_code == 0, f"Process failed: {result.output}"

        # Check directory structure
        assert (working_dir_path / "metadata").is_dir()
        assert (working_dir_path / "metadata" / "metadata.json").is_file()
        assert (working_dir_path / "logs").is_dir()

        # Check that progress logs were created
        log_files = list((working_dir_path / "logs").glob("*.cjson"))
        assert len(log_files) > 0, "No progress log files created"

    def test_process_submission_multi_inbox(
        self,
        aws_credentials,
        tmp_path,
        working_dir_path,
        process_config_content,
    ):
        """
        Test multi-inbox selection logic.
        """
        conn = boto3.client("s3")
        conn.create_bucket(Bucket="inbox-a")
        conn.create_bucket(Bucket="inbox-b")
        conn.create_bucket(Bucket="consented-archive")
        conn.create_bucket(Bucket="non-consented-archive")
        conn.create_bucket(Bucket="interrogation-archive")

        s3_resource = boto3.resource("s3")
        inbox_b = s3_resource.Bucket("inbox-b")

        submission_id = "260914050_2024-07-15_c64603a7"
        le_id = "260914050"

        # Upload submission to inbox-b
        upload_submission_to_inbox(inbox_b, submission_id)
        config_content = process_config_content.copy()

        base_inbox_config = config_content["leistungserbringer"][le_id]["inbox_buckets"]["inbox"]

        config_content["leistungserbringer"][le_id]["inbox_buckets"] = {
            "inbox-a": base_inbox_config.copy(),
            "inbox-b": base_inbox_config.copy(),
        }

        config_file = tmp_path / "config.multi.yaml"
        with open(config_file, "w") as fd:
            yaml.dump(config_content, fd)

        # Run grzctl process with --inbox-bucket inbox-b
        result = _run_process(
            config_file, submission_id, working_dir_path, "--inbox-bucket", "inbox-b", "--no-update-db"
        )

        assert result.exit_code == 0, f"Process failed: {result.output}"

        # Verify it worked (files should be in consented archive)
        consented_bucket = s3_resource.Bucket("consented-archive")
        consented_keys = {o.key for o in consented_bucket.objects.all()}
        assert any("metadata/metadata.json" in key for key in consented_keys)

        # Test failure if --inbox-bucket is missing when multiple are available
        result = _run_process(config_file, submission_id, working_dir_path / "fail", "--no-update-db")
        assert result.exit_code != 0
        assert "Multiple inboxes found" in result.output


class TestProcessValidationFailure:
    """Tests for validation failure handling during processing."""

    def test_validation_failure_aborts_upload(
        self,
        s3_buckets,
        temp_process_config_file_path,
        crypt4gh_grz_public_key_file_path,
        crypt4gh_submitter_private_key_file_path,
        working_dir_path,
    ):
        """A FASTQ file that fails format validation must fail the run and reach no archive."""
        submission_id = "260914050_2024-07-15_c64603a7"
        read1 = "aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_normal.read1.fastq.gz"
        upload_submission_to_inbox(s3_buckets["inbox"], submission_id)

        # replace read1 with a FASTQ whose only record lacks its quality line
        invalid_fastq = gzip.compress(b"@read1\nACGT\n+\n")
        sender_private_key = crypt4gh.keys.get_private_key(str(crypt4gh_submitter_private_key_file_path), lambda: None)
        grz_public_key = crypt4gh.keys.get_public_key(str(crypt4gh_grz_public_key_file_path))
        encrypted = BytesIO()
        crypt4gh.lib.encrypt([(0, sender_private_key, grz_public_key)], BytesIO(invalid_fastq), encrypted)
        s3_buckets["inbox"].put_object(Key=f"{submission_id}/files/{read1}.c4gh", Body=encrypted.getvalue())

        # update its checksum and size, so that it fails format validation rather than the checksum check
        metadata = json.loads((VALID_SUBMISSION_DIR / "metadata" / "metadata.json").read_text())
        for donor in metadata["donors"]:
            for lab_datum in donor["labData"]:
                for file in lab_datum.get("sequenceData", {}).get("files", []):
                    if file["filePath"] == read1:
                        file["fileChecksum"] = hashlib.sha256(invalid_fastq).hexdigest()
                        file["fileSizeInBytes"] = len(invalid_fastq)
        s3_buckets["inbox"].put_object(
            Key=f"{submission_id}/metadata/metadata.json", Body=json.dumps(metadata).encode()
        )

        result = _run_process(temp_process_config_file_path, submission_id, working_dir_path, "--no-update-db")

        assert result.exit_code != 0, f"Process should have failed but succeeded: {result.output}"
        logs = [path.read_text() for path in (working_dir_path / "logs").iterdir()]
        assert any("sequence length (4) does not match quality scores length (0)" in log for log in logs)

        for bucket in ("consented", "non_consented", "interrogation"):
            keys = {o.key for o in s3_buckets[bucket].objects.all()}
            assert not keys, f"The {bucket} bucket should be empty, has: {keys}"

    def test_checksum_mismatch_fails_processing(
        self,
        s3_buckets,
        temp_process_config_file_path,
        working_dir_path,
    ):
        """A file whose content does not match its metadata checksum must fail the run and reach no archive."""
        submission_id = "260914050_2024-07-15_c64603a7"
        upload_submission_to_inbox(s3_buckets["inbox"], submission_id)

        # replace the metadata with a copy in which one file has a wrong checksum; the file is
        # listed under several donors, so change every entry
        metadata = json.loads((VALID_SUBMISSION_DIR / "metadata" / "metadata.json").read_text())
        for donor in metadata["donors"]:
            for lab_datum in donor["labData"]:
                for file in lab_datum.get("sequenceData", {}).get("files", []):
                    if file["filePath"] == "target_regions.bed":
                        file["fileChecksum"] = "0" * 64
        s3_buckets["inbox"].put_object(
            Key=f"{submission_id}/metadata/metadata.json", Body=json.dumps(metadata).encode()
        )

        # with two threads, a second file gets staged while the first one fails, so the cleanup is tested
        result = _run_process(
            temp_process_config_file_path, submission_id, working_dir_path, "--no-update-db", "--threads", "2"
        )

        assert result.exit_code != 0, f"Process should have failed but succeeded: {result.output}"
        logs = [path.read_text() for path in (working_dir_path / "logs").iterdir()]
        assert any("Checksum mismatch" in log for log in logs)

        for bucket in ("consented", "non_consented", "interrogation"):
            keys = {o.key for o in s3_buckets[bucket].objects.all()}
            assert not keys, f"The {bucket} bucket should be empty, has: {keys}"


@pytest.fixture
def qc_process_config_file_path(tmp_path, process_config_content) -> Path:
    """Write a process config that selects submissions for detailed QC."""
    process_config_content["detailed_qc"]["target_percentage"] = "100.0"
    config_file = tmp_path / "config.process.qc.yaml"
    with open(config_file, "w") as fd:
        yaml.dump(process_config_content, fd)
    return config_file


def _upload_initial_submission_to_inbox(inbox_bucket, submission_id: str) -> None:
    """Upload the valid submission as an initial one, since only initial submissions are selected for QC."""
    upload_submission_to_inbox(inbox_bucket, submission_id)
    metadata = json.loads((VALID_SUBMISSION_DIR / "metadata" / "metadata.json").read_text())
    metadata["submission"]["submissionType"] = "initial"
    inbox_bucket.put_object(Key=f"{submission_id}/metadata/metadata.json", Body=json.dumps(metadata).encode())


def _qc_files(process_config_content: dict, submission_id: str) -> set[str]:
    files_dir = Path(process_config_content["detailed_qc"]["local_storage"]) / submission_id / "files"
    return {p.relative_to(files_dir).as_posix() for p in files_dir.rglob("*") if p.is_file()}


def _metadata_file_checksums() -> dict[str, str]:
    """Map each file path in the valid submission's metadata to its SHA256 checksum."""
    metadata = json.loads((VALID_SUBMISSION_DIR / "metadata" / "metadata.json").read_text())
    return {
        file["filePath"]: file["fileChecksum"]
        for donor in metadata["donors"]
        for lab_datum in donor["labData"]
        for file in lab_datum.get("sequenceData", {}).get("files", [])
    }


def _assert_archived(bucket, submission_id: str, private_key_path: Path) -> None:
    """Assert that ``bucket`` holds every file of the submission, each decrypting to its metadata checksum.

    Decrypts with the reference crypt4gh implementation rather than grzctl's own decryptor.
    """
    expected = _metadata_file_checksums()
    archived_keys = {o.key for o in bucket.objects.filter(Prefix=f"{submission_id}/files/")}
    assert archived_keys == {f"{submission_id}/files/{path}.c4gh" for path in expected}
    private_key = crypt4gh.keys.get_private_key(str(private_key_path), lambda: None)
    for path, checksum in expected.items():
        encrypted = bucket.Object(f"{submission_id}/files/{path}.c4gh").get()["Body"].read()
        decrypted = BytesIO()
        crypt4gh.lib.decrypt([(0, private_key, None)], BytesIO(encrypted), decrypted)
        assert hashlib.sha256(decrypted.getvalue()).hexdigest() == checksum, path


S3_WRITE_OPERATIONS = frozenset(
    {"PutObject", "CopyObject", "CreateMultipartUpload", "UploadPart", "UploadPartCopy", "CompleteMultipartUpload"}
)


@dataclass
class S3Requests:
    """S3 requests as botocore sends them, whichever code in grzctl makes them."""

    requests: list[tuple[str, str, str]] = field(default_factory=list)
    unavailable_bucket: str | None = None
    """Writes to this bucket fail, as if the bucket were unavailable."""

    def per_file(self, operations: set[str], bucket, submission_id: str) -> Counter[str]:
        """Count the requests with one of ``operations`` on each of the submission's files in ``bucket``."""
        prefix = f"{submission_id}/files/"
        return Counter(
            key.removeprefix(prefix).removesuffix(".c4gh")
            for operation, bucket_name, key in self.requests
            if operation in operations and bucket_name == bucket.name and key.startswith(prefix)
        )


@pytest.fixture
def s3_requests(monkeypatch) -> S3Requests:
    """Record every S3 request, and fail writes to :attr:`S3Requests.unavailable_bucket`."""
    recorder = S3Requests()
    make_api_call = botocore.client.BaseClient._make_api_call

    def record(client, operation_name, api_params):
        bucket = api_params.get("Bucket", "")
        recorder.requests.append((operation_name, bucket, api_params.get("Key", "")))
        if bucket == recorder.unavailable_bucket and operation_name in S3_WRITE_OPERATIONS:
            error = {"Error": {"Code": "ServiceUnavailable", "Message": "simulated outage"}}
            raise botocore.exceptions.ClientError(error, operation_name)
        return make_api_call(client, operation_name, api_params)

    monkeypatch.setattr(botocore.client.BaseClient, "_make_api_call", record)
    return recorder


class TestProcessDetailedQc:
    """Tests for the detailed QC prefetch: the main pass writes the QC copy when a selection is likely."""

    SUBMISSION_ID = "260914050_2024-07-15_c64603a7"

    def test_prefetch_replaces_the_qc_pass(
        self,
        s3_buckets,
        s3_requests,
        qc_process_config_file_path,
        process_config_content,
        working_dir_path,
    ):
        """When the guess is right, each file is downloaded once, and the QC pass finds it already in QC storage."""
        _upload_initial_submission_to_inbox(s3_buckets["inbox"], self.SUBMISSION_ID)

        result = _run_process(qc_process_config_file_path, self.SUBMISSION_ID, working_dir_path, "--update-db")

        assert result.exit_code == 0, f"Process failed: {result.output}"
        checksums = _metadata_file_checksums()
        downloads = s3_requests.per_file({"GetObject"}, s3_buckets["inbox"], self.SUBMISSION_ID)
        assert downloads == Counter(dict.fromkeys(checksums, 1))
        assert _qc_files(process_config_content, self.SUBMISSION_ID) == set(checksums)
        metadata = json.loads((VALID_SUBMISSION_DIR / "metadata" / "metadata.json").read_text())
        qc_dir = Path(process_config_content["detailed_qc"]["local_storage"]) / self.SUBMISSION_ID
        uploaded_metadata = {**metadata, "submission": {**metadata["submission"], "submissionType": "initial"}}
        assert json.loads((qc_dir / "metadata" / "metadata.json").read_text()) == uploaded_metadata

    def test_prefetched_files_are_deleted_when_not_selected(
        self,
        s3_buckets,
        qc_process_config_file_path,
        process_config_content,
        working_dir_path,
        monkeypatch,
    ):
        """A wrong guess leaves no decrypted files behind."""
        upload_submission_to_inbox(s3_buckets["inbox"], self.SUBMISSION_ID)

        # guess "selected", then decide "not selected"
        monkeypatch.setattr(SubmissionDb, "should_qc", lambda self, *args, predict=False, **kwargs: predict)

        result = _run_process(qc_process_config_file_path, self.SUBMISSION_ID, working_dir_path, "--update-db")

        assert result.exit_code == 0, f"Process failed: {result.output}"
        assert _qc_files(process_config_content, self.SUBMISSION_ID) == set()

    def test_no_update_db_skips_qc_selection(
        self,
        s3_buckets,
        qc_process_config_file_path,
        process_config_content,
        working_dir_path,
    ):
        """The selection stores its decision in the DB, so ``--no-update-db`` runs without detailed QC."""
        _upload_initial_submission_to_inbox(s3_buckets["inbox"], self.SUBMISSION_ID)

        result = _run_process(qc_process_config_file_path, self.SUBMISSION_ID, working_dir_path, "--no-update-db")

        assert result.exit_code == 0, f"Process failed: {result.output}"
        assert _qc_files(process_config_content, self.SUBMISSION_ID) == set()


class TestProcessDuplicateInitial:
    """An initial submission fails basic QC when its case already has a QC-passed initial submission."""

    FIRST_ID = "260914050_2024-07-15_c64603a7"
    DUPLICATE_TAN_G = "bbbbbbbb00000000bbbbbbbb00000000bbbbbbbb00000000bbbbbbbb00000000"

    def _process_first_and_upload_duplicate(self, inbox_bucket, config_file_path: Path, working_dir_path: Path) -> str:
        """Process an initial submission with ``--update-db``, then upload a second initial submission of its case."""
        _upload_initial_submission_to_inbox(inbox_bucket, self.FIRST_ID)
        result = _run_process(config_file_path, self.FIRST_ID, working_dir_path / "first", "--update-db")
        assert result.exit_code == 0, f"Process failed: {result.output}"

        # same submitter and local case ID, but its own tanG and therefore its own submission ID
        duplicate_id = f"260914050_2024-07-15_{hashlib.sha256(self.DUPLICATE_TAN_G.encode()).hexdigest()[:8]}"
        upload_submission_to_inbox(inbox_bucket, duplicate_id)
        metadata = json.loads((VALID_SUBMISSION_DIR / "metadata" / "metadata.json").read_text())
        metadata["submission"]["submissionType"] = "initial"
        metadata["submission"]["tanG"] = self.DUPLICATE_TAN_G
        inbox_bucket.put_object(Key=f"{duplicate_id}/metadata/metadata.json", Body=json.dumps(metadata).encode())
        return duplicate_id

    @staticmethod
    def _assert_failed_basic_qc(process_config_content: dict, submission_id: str) -> None:
        db = SubmissionDb(db_url=process_config_content["db"]["database_url"], author=None)
        submission = db.get_submission(submission_id)
        assert submission is not None
        assert submission.basic_qc_passed is False
        assert submission.states[-1].state == SubmissionStateEnum.ERROR
        assert submission.states[-1].failure_reason == FailureReasonEnum.DUPLICATE_INITIAL

    def test_update_db_fails_basic_qc_before_processing(
        self,
        s3_buckets,
        s3_requests,
        temp_process_config_file_path,
        process_config_content,
        working_dir_path,
    ):
        """The pre-check rejects the duplicate before any of its files is downloaded."""
        duplicate_id = self._process_first_and_upload_duplicate(
            s3_buckets["inbox"], temp_process_config_file_path, working_dir_path
        )

        result = _run_process(
            temp_process_config_file_path, duplicate_id, working_dir_path / "duplicate", "--update-db"
        )

        assert result.exit_code != 0
        assert s3_requests.per_file({"GetObject"}, s3_buckets["inbox"], duplicate_id) == Counter()
        self._assert_failed_basic_qc(process_config_content, duplicate_id)

    def test_update_db_records_duplicate_detected_after_processing(
        self,
        s3_buckets,
        temp_process_config_file_path,
        process_config_content,
        working_dir_path,
        monkeypatch,
    ):
        """The ``basic_qc_passed`` write catches a competing initial submission that the pre-check missed."""
        duplicate_id = self._process_first_and_upload_duplicate(
            s3_buckets["inbox"], temp_process_config_file_path, working_dir_path
        )
        monkeypatch.setattr(SubmissionDb, "assert_no_duplicate_initial", lambda self, *args, **kwargs: None)

        result = _run_process(
            temp_process_config_file_path, duplicate_id, working_dir_path / "duplicate", "--update-db"
        )

        assert result.exit_code != 0
        self._assert_failed_basic_qc(process_config_content, duplicate_id)


@pytest.fixture
def rerun_config_file_path(tmp_path, process_config_content) -> Path:
    """Write a process config that selects submissions for detailed QC and keeps the staged files of a failed run."""
    process_config_content["detailed_qc"]["target_percentage"] = "100.0"
    process_config_content["archives"]["interrogation"]["keep_failed"] = True
    config_file = tmp_path / "config.process.rerun.yaml"
    config_file.write_text(yaml.dump(process_config_content))
    return config_file


class TestProcessRerun:
    """A rerun into the same output directory redoes exactly the work whose output is missing.

    These tests observe only what an operator could: the exit code, the S3 requests, the buckets and
    local storage. They do not depend on how grzctl tracks progress.
    """

    SUBMISSION_ID = "260914050_2024-07-15_c64603a7"
    VCF = "aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_normal.vcf"
    READ1 = "aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_normal.read1.fastq.gz"
    READ2 = "aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_normal.read2.fastq.gz"

    @pytest.mark.parametrize("staged", [True, False], ids=["staged", "not-staged"])
    @pytest.mark.parametrize("on_local_storage", [True, False], ids=["on-local-storage", "not-on-local-storage"])
    def test_rerun_redoes_only_missing_outputs(
        self,
        s3_buckets,
        s3_requests,
        rerun_config_file_path,
        process_config_content,
        crypt4gh_grz_private_key_file_path,
        working_dir_path,
        staged,
        on_local_storage,
    ):
        """A rerun downloads, stages and writes a file to local storage only as far as its outputs are missing.

        The first run fails while copying to the archive. By then every file is staged in the
        interrogation bucket and, since the submission is selected for detailed QC, on local storage.
        The test then removes the staged copy, the local copy, or both, of one file.
        """
        sid = self.SUBMISSION_ID
        checksums = _metadata_file_checksums()
        local_files = Path(process_config_content["detailed_qc"]["local_storage"]) / sid / "files"
        _upload_initial_submission_to_inbox(s3_buckets["inbox"], sid)

        s3_requests.unavailable_bucket = s3_buckets["consented"].name
        result = _run_process(rerun_config_file_path, sid, working_dir_path, "--update-db")
        s3_requests.unavailable_bucket = None
        assert result.exit_code != 0, "the first run should fail while copying to the archive"
        staged_keys = {o.key for o in s3_buckets["interrogation"].objects.filter(Prefix=f"{sid}/files/")}
        assert staged_keys == {f"{sid}/files/{path}.c4gh" for path in checksums}
        assert _qc_files(process_config_content, sid) == set(checksums)

        if not staged:
            s3_buckets["interrogation"].Object(f"{sid}/files/{self.VCF}.c4gh").delete()
        if not on_local_storage:
            (local_files / self.VCF).unlink()
        for path in local_files.iterdir():
            os.utime(path, ns=(0, 0))  # a local copy the rerun writes gets a newer modification time
        s3_requests.requests.clear()

        result = _run_process(rerun_config_file_path, sid, working_dir_path, "--update-db")

        assert result.exit_code == 0, f"Rerun failed: {result.output}"
        outputs_missing = not (staged and on_local_storage)
        downloads = s3_requests.per_file({"GetObject"}, s3_buckets["inbox"], sid)
        assert downloads == Counter({self.VCF: 1} if outputs_missing else {})
        uploads = s3_requests.per_file({"PutObject", "CompleteMultipartUpload"}, s3_buckets["interrogation"], sid)
        assert uploads == Counter({} if staged else {self.VCF: 1})
        written = {path.name for path in local_files.iterdir() if path.stat().st_mtime_ns != 0}
        assert written == (set() if on_local_storage else {self.VCF})
        for path, checksum in checksums.items():
            assert hashlib.sha256((local_files / path).read_bytes()).hexdigest() == checksum, path
        _assert_archived(s3_buckets["consented"], sid, crypt4gh_grz_private_key_file_path)

    @pytest.mark.parametrize("keep_failed", [True, False], ids=["staged-files-kept", "staged-files-deleted"])
    def test_rerun_completes_a_pair_whose_read2_failed(
        self,
        s3_buckets,
        s3_requests,
        tmp_path,
        process_config_content,
        crypt4gh_grz_private_key_file_path,
        working_dir_path,
        keep_failed,
    ):
        """The first run processes read1 and fails on the missing read2; the rerun completes the pair.

        With ``keep_failed`` the staged read1 survives, so the rerun must not download it again, and
        the read-pair check of read2 must still pass. Without it, the rerun downloads read1 again.
        """
        sid = self.SUBMISSION_ID
        process_config_content["archives"]["interrogation"]["keep_failed"] = keep_failed
        config_file_path = tmp_path / "config.process.yaml"
        config_file_path.write_text(yaml.dump(process_config_content))
        upload_submission_to_inbox(s3_buckets["inbox"], sid)
        read2_key = f"{sid}/files/{self.READ2}.c4gh"
        s3_buckets["inbox"].Object(read2_key).delete()

        # with one thread, read1 comes before read2; the second assertion checks that it did
        result = _run_process(config_file_path, sid, working_dir_path, "--no-update-db", "--threads", "1")
        assert result.exit_code != 0, "the first run should fail on the missing read2"
        downloads = s3_requests.per_file({"GetObject"}, s3_buckets["inbox"], sid)
        assert downloads[self.READ1] == 1, "the first run should download read1"

        s3_buckets["inbox"].upload_file(
            Filename=str(VALID_SUBMISSION_DIR / "encrypted_files" / f"{self.READ2}.c4gh"), Key=read2_key
        )
        s3_requests.requests.clear()
        result = _run_process(config_file_path, sid, working_dir_path, "--no-update-db", "--threads", "1")

        assert result.exit_code == 0, f"Rerun failed: {result.output}"
        downloads = s3_requests.per_file({"GetObject"}, s3_buckets["inbox"], sid)
        assert downloads[self.READ1] == (0 if keep_failed else 1)
        assert downloads[self.READ2] == 1
        _assert_archived(s3_buckets["consented"], sid, crypt4gh_grz_private_key_file_path)
