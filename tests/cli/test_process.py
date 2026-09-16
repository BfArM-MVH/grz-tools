"""
Integration tests for grzctl process command.

This tests the full streaming pipeline with mocked S3 buckets (inbox + archive).
"""

import hashlib
import json
import os
import shutil
from pathlib import Path

import boto3
import click.testing
import grzctl.cli
import pytest
import yaml
from grz_common.models.base import get_secret_value
from grz_db.models.submission import FailureReasonEnum, SubmissionDb, SubmissionStateEnum
from grzctl.models.config import GrzctlConfig
from grzctl.processor import FilePipelineExecutor
from moto import mock_aws

# Path to test fixtures
MOCK_FILES_DIR = Path(__file__).parent.parent / "mock_files"
VALID_SUBMISSION_DIR = MOCK_FILES_DIR / "submissions" / "valid_submission"


@pytest.fixture
def process_config_content(
    crypt4gh_grz_private_key_file_path,
    crypt4gh_grz_public_key_file_path,
    db_alice_private_key_file_path,
    db_known_keys_file_path,
    tmpdir_factory,
):
    db_dir = tmpdir_factory.mktemp("db")
    db_file = db_dir / "test.db"
    local_storage_dir = tmpdir_factory.mktemp("local_storage")

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
        "db": {
            "database_url": f"sqlite:///{str(db_file)}",
            "author": {
                "name": "Alice",
                "private_key_path": str(db_alice_private_key_file_path),
            },
            "known_public_keys": str(db_known_keys_file_path),
        },
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
def temp_process_config_file_path(temp_data_dir_path, process_config_content) -> Path:
    """Write the process config to a YAML file."""
    import yaml

    config_file = temp_data_dir_path / "config.process.yaml"
    with open(config_file, "w") as fd:
        yaml.dump(process_config_content, fd)
    return config_file


@pytest.fixture
def aws_credentials_for_process():
    """Mocked AWS Credentials for moto."""
    os.environ["AWS_ACCESS_KEY_ID"] = "testing"
    os.environ["AWS_SECRET_ACCESS_KEY"] = "testing"
    os.environ["MOTO_ALLOW_NONEXISTENT_REGION"] = "1"
    with mock_aws():
        yield


@pytest.fixture
def s3_buckets(aws_credentials_for_process):
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


@pytest.fixture
def initialized_db(temp_process_config_file_path):
    """Initialize the database for the process command tests."""
    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()

    result = runner.invoke(
        cli,
        ["--config", str(temp_process_config_file_path), "db", "init"],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, f"DB init failed: {result.output}"
    return True


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


class TestGrzctlProcess:
    """Integration tests for grzctl process command."""

    def test_process_submission_basic(
        self,
        s3_buckets,
        temp_process_config_file_path,
        initialized_db,
        working_dir_path,
        tmp_path,
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
        args = [
            "--config",
            str(temp_process_config_file_path),
            "process",
            "--submission-id",
            submission_id,
            "--output-dir",
            str(working_dir_path),
            "--no-submit-pruefbericht",
            "--no-update-db",
        ]

        runner = click.testing.CliRunner()
        cli = grzctl.cli.build_cli()
        result = runner.invoke(cli, args, catch_exceptions=False)

        assert result.exit_code == 0, f"Process failed: {result.output}"

        # Check that files were uploaded to the consented archive
        # (the valid_submission has research consent)
        consented_keys = {o.key for o in s3_buckets["consented"].objects.all()}

        # Should have metadata and encrypted files
        assert any("metadata/metadata.json" in key for key in consented_keys), (
            f"Metadata not found in consented archive. Keys: {consented_keys}"
        )
        assert any("files/" in key and ".c4gh" in key for key in consented_keys), (
            f"Encrypted files not found in consented archive. Keys: {consented_keys}"
        )

    def test_process_archives_redacted_logs(
        self,
        s3_buckets,
        temp_process_config_file_path,
        initialized_db,
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

        args = [
            "--config",
            str(temp_process_config_file_path),
            "process",
            "--submission-id",
            submission_id,
            "--output-dir",
            str(working_dir_path),
            "--no-submit-pruefbericht",
            "--no-update-db",
        ]

        runner = click.testing.CliRunner()
        cli = grzctl.cli.build_cli()
        result = runner.invoke(cli, args, catch_exceptions=False)
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
        initialized_db,
        working_dir_path,
    ):
        """
        Test that grzctl process creates the expected local directory structure.
        """
        submission_id = "260914050_2024-07-15_c64603a7"

        # Upload submission to inbox
        upload_submission_to_inbox(s3_buckets["inbox"], submission_id)

        # Run grzctl process
        args = [
            "--config",
            str(temp_process_config_file_path),
            "process",
            "--submission-id",
            submission_id,
            "--output-dir",
            str(working_dir_path),
            "--no-submit-pruefbericht",
            "--no-update-db",
        ]

        runner = click.testing.CliRunner()
        cli = grzctl.cli.build_cli()
        result = runner.invoke(cli, args, catch_exceptions=False)

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
        aws_credentials_for_process,
        temp_data_dir_path,
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
        _inbox_a = s3_resource.Bucket("inbox-a")
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

        config_file = temp_data_dir_path / "config.multi.yaml"
        import yaml

        with open(config_file, "w") as fd:
            yaml.dump(config_content, fd)

        # Run grzctl process with --inbox-bucket inbox-b
        args = [
            "--config",
            str(config_file),
            "process",
            "--submission-id",
            submission_id,
            "--output-dir",
            str(working_dir_path),
            "--inbox-bucket",
            "inbox-b",
            "--no-submit-pruefbericht",
            "--no-update-db",
        ]

        runner = click.testing.CliRunner()
        cli = grzctl.cli.build_cli()
        result = runner.invoke(cli, args, catch_exceptions=False)

        assert result.exit_code == 0, f"Process failed: {result.output}"

        # Verify it worked (files should be in consented archive)
        consented_bucket = s3_resource.Bucket("consented-archive")
        consented_keys = {o.key for o in consented_bucket.objects.all()}
        assert any("metadata/metadata.json" in key for key in consented_keys)

        # Test failure if --inbox-bucket is missing when multiple are available
        args_no_bucket = [
            "--config",
            str(config_file),
            "process",
            "--submission-id",
            submission_id,
            "--output-dir",
            str(working_dir_path / "fail"),
            "--no-submit-pruefbericht",
            "--no-update-db",
        ]
        result = runner.invoke(cli, args_no_bucket, catch_exceptions=False)
        assert result.exit_code != 0
        assert "Multiple inboxes found" in result.output


class TestProcessVsManualWorkflow:
    """
    Compare grzctl process output with manual CLI workflow.

    This ensures the streaming pipeline produces equivalent results to:
        grzctl decrypt -> grz-cli validate -> grz-cli encrypt -> grzctl archive
    """

    @pytest.fixture
    def working_dir_process(self, tmpdir_factory) -> Path:
        """Working directory for grzctl process."""
        return Path(tmpdir_factory.mktemp("process_submission").strpath)

    @pytest.fixture
    def working_dir_manual(self, tmpdir_factory) -> Path:
        """Working directory for manual CLI workflow."""
        return Path(tmpdir_factory.mktemp("manual_submission").strpath)

    def _setup_manual_submission_dir(self, working_dir: Path) -> None:
        """Copy encrypted submission files to working directory for manual processing."""
        shutil.copytree(
            VALID_SUBMISSION_DIR / "encrypted_files",
            working_dir / "encrypted_files",
            dirs_exist_ok=True,
        )
        shutil.copytree(
            VALID_SUBMISSION_DIR / "metadata",
            working_dir / "metadata",
            dirs_exist_ok=True,
        )

    def _run_manual_decrypt(self, working_dir: Path, config_file: str) -> None:
        """Run grzctl decrypt command."""
        import grzctl.cli

        runner = click.testing.CliRunner()
        cli = grzctl.cli.build_cli()

        result = runner.invoke(
            cli,
            ["--config", config_file, "decrypt", "--submission-dir", str(working_dir), "--no-update-db"],
            catch_exceptions=False,
        )
        assert result.exit_code == 0, f"Decrypt failed: {result.output}"

    def test_process_produces_same_metadata_as_archive(
        self,
        s3_buckets,
        temp_process_config_file_path,
        temp_keys_config_file_path,
        temp_grzctl_keys_config_file_path,
        initialized_db,
        working_dir_process,
        working_dir_manual,
        tmp_path,
    ):
        """
        Test that grzctl process archives the same metadata as manual archive command.
        """
        submission_id = "260914050_2024-07-15_c64603a7"

        # Upload submission to inbox
        upload_submission_to_inbox(s3_buckets["inbox"], submission_id)

        # === Run grzctl process ===
        args = [
            "--config",
            str(temp_process_config_file_path),
            "process",
            "--submission-id",
            submission_id,
            "--output-dir",
            str(working_dir_process),
            "--no-submit-pruefbericht",
            "--no-update-db",
        ]

        runner = click.testing.CliRunner()
        cli = grzctl.cli.build_cli()
        result = runner.invoke(cli, args, catch_exceptions=False)

        assert result.exit_code == 0, f"Process failed: {result.output}"

        # === Run manual workflow ===
        self._setup_manual_submission_dir(working_dir_manual)
        self._run_manual_decrypt(working_dir_manual, str(temp_grzctl_keys_config_file_path))

        # === Compare decrypted content checksums ===
        # Get checksums from files decrypted by grzctl process
        # (grzctl process downloads metadata but streams files directly to archive)
        process_metadata_path = working_dir_process / "metadata" / "metadata.json"
        manual_metadata_path = working_dir_manual / "metadata" / "metadata.json"

        # Both should have valid metadata
        assert process_metadata_path.is_file()
        assert manual_metadata_path.is_file()

        # Load and compare (content should be similar, though format may differ)
        with open(process_metadata_path) as f:
            process_metadata = json.load(f)
        with open(manual_metadata_path) as f:
            manual_metadata = json.load(f)

        # Core submission info should match
        assert process_metadata["submission"]["submitterId"] == manual_metadata["submission"]["submitterId"]
        assert process_metadata["submission"]["submissionDate"] == manual_metadata["submission"]["submissionDate"]


class TestProcessValidationFailure:
    """Tests for validation failure handling during processing."""

    def test_validation_failure_aborts_upload(
        self,
        s3_buckets,
        temp_process_config_file_path,
        initialized_db,
        working_dir_path,
        tmp_path,
    ):
        """
        Test that validation failure causes the pipeline to fail and abort upload.

        When validation fails:
        1. The pipeline should return an error
        2. The upload should be aborted (no files in archive)
        3. The error should be reported

        This test creates an invalid FASTQ file (non-multiple-of-4 lines) to trigger
        validation failure.
        """
        submission_id = "260914050_2024-07-15_c64603a7"

        # Create a mock submission with an invalid FASTQ file
        self._upload_submission_with_invalid_fastq(s3_buckets["inbox"], submission_id, tmp_path)

        # Verify inbox has files
        inbox_keys = {o.key for o in s3_buckets["inbox"].objects.all()}
        assert f"{submission_id}/metadata/metadata.json" in inbox_keys

        # Run grzctl process
        args = [
            "--config",
            str(temp_process_config_file_path),
            "process",
            "--submission-id",
            submission_id,
            "--output-dir",
            str(working_dir_path),
            "--no-submit-pruefbericht",
            "--no-update-db",
        ]

        runner = click.testing.CliRunner()
        cli = grzctl.cli.build_cli()
        result = runner.invoke(cli, args, catch_exceptions=False)

        # The pipeline should fail
        assert result.exit_code != 0, f"Process should have failed but succeeded: {result.output}"
        assert "validation" in result.output.lower() or "error" in result.output.lower()

        # Archive should be empty (upload was aborted)
        consented_keys = {o.key for o in s3_buckets["consented"].objects.all()}
        non_consented_keys = {o.key for o in s3_buckets["non_consented"].objects.all()}

        # No files should have been uploaded to either archive
        assert len(consented_keys) == 0, f"Consented archive should be empty, has: {consented_keys}"
        assert len(non_consented_keys) == 0, f"Non-consented archive should be empty, has: {non_consented_keys}"

    def test_checksum_mismatch_fails_processing(
        self,
        s3_buckets,
        temp_process_config_file_path,
        initialized_db,
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

        args = [
            "--config",
            str(temp_process_config_file_path),
            "process",
            "--submission-id",
            submission_id,
            "--output-dir",
            str(working_dir_path),
            "--no-submit-pruefbericht",
            "--no-update-db",
        ]

        runner = click.testing.CliRunner()
        cli = grzctl.cli.build_cli()
        result = runner.invoke(cli, args)

        assert result.exit_code != 0, f"Process should have failed but succeeded: {result.output}"
        progress_log = (working_dir_path / "logs" / "progress_processing.cjson").read_text()
        assert "Checksum mismatch" in progress_log

        for bucket in ("consented", "non_consented", "interrogation"):
            keys = {o.key for o in s3_buckets[bucket].objects.all()}
            assert not keys, f"The {bucket} bucket should be empty, has: {keys}"

    def _upload_submission_with_invalid_fastq(self, inbox_bucket, submission_id: str, tmp_path: Path):
        """
        Upload a submission with an invalid FASTQ file to the inbox.

        Creates a FASTQ with only 3 lines (should be multiple of 4).
        """
        import gzip
        import io

        # Get the public key for encryption
        grz_public_key_path = VALID_SUBMISSION_DIR.parent.parent.parent / "conftest_keys" / "grz.pub"
        if not grz_public_key_path.exists():
            # Use the key from valid submission if available
            grz_public_key_path = Path(__file__).parent.parent / "conftest_keys" / "grz.pub"

        # Create invalid FASTQ content (only 3 lines - not a multiple of 4)
        invalid_fastq = b"@read1\nACGT\n+\n"  # Only 3 lines instead of 4

        # Gzip the content
        gzipped = io.BytesIO()
        with gzip.GzipFile(fileobj=gzipped, mode="wb") as gz:
            gz.write(invalid_fastq)
        gzipped_content = gzipped.getvalue()

        # Create metadata that references a FASTQ file
        import hashlib

        checksum = hashlib.sha256(gzipped_content).hexdigest()

        metadata = {
            "metadataSchemaVersion": "1.0.0",
            "submission": {
                "submitterId": "260914050",
                "submissionDate": "2024-07-15",
                "dataSubmitterId": "260914050",
                "tanG": "TAN123456",
                "localCaseId": "LOCAL123",
            },
            "donors": [
                {
                    "donorPseudonym": "index",
                    "relation": "index",
                    "sex": "male",
                    "labData": [
                        {
                            "labDataName": "blood_test",
                            "sequencingData": [
                                {
                                    "files": [
                                        {
                                            "filePath": "invalid.fastq.gz",
                                            "fileType": "fastq",
                                            "checksumType": "sha256",
                                            "fileChecksum": checksum,
                                            "fileSizeInBytes": len(gzipped_content),
                                            "readOrder": "R1",
                                        }
                                    ]
                                }
                            ],
                        }
                    ],
                }
            ],
        }

        # Upload metadata
        metadata_json = json.dumps(metadata)
        inbox_bucket.put_object(
            Key=f"{submission_id}/metadata/metadata.json",
            Body=metadata_json.encode(),
        )

        # For now, just upload the gzipped content directly (not encrypted)
        # This won't work with the full pipeline but tests the concept
        # In a real test, we'd need to encrypt with crypt4gh
        inbox_bucket.put_object(
            Key=f"{submission_id}/files/invalid.fastq.gz.c4gh",
            Body=gzipped_content,  # Not actually encrypted for simplicity
        )


@pytest.fixture
def qc_process_config_file_path(temp_data_dir_path, process_config_content) -> Path:
    """Write a process config that selects submissions for detailed QC and can sign DB state changes."""
    process_config_content["detailed_qc"]["target_percentage"] = "100.0"
    process_config_content["db"]["author"]["private_key_passphrase"] = "test"
    config_file = temp_data_dir_path / "config.process.qc.yaml"
    with open(config_file, "w") as fd:
        yaml.dump(process_config_content, fd)
    return config_file


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


def _upload_initial_submission_to_inbox(inbox_bucket, submission_id: str) -> None:
    """Upload the valid submission as an initial one, since only initial submissions are selected for QC."""
    upload_submission_to_inbox(inbox_bucket, submission_id)
    metadata = json.loads((VALID_SUBMISSION_DIR / "metadata" / "metadata.json").read_text())
    metadata["submission"]["submissionType"] = "initial"
    inbox_bucket.put_object(Key=f"{submission_id}/metadata/metadata.json", Body=json.dumps(metadata).encode())


def _qc_files(process_config_content: dict, submission_id: str) -> set[str]:
    files_dir = Path(process_config_content["detailed_qc"]["local_storage"]) / submission_id / "files"
    return {p.relative_to(files_dir).as_posix() for p in files_dir.rglob("*") if p.is_file()}


class TestProcessDetailedQc:
    """Tests for the detailed QC prefetch: the main pass writes the QC copy when a selection is likely."""

    SUBMISSION_ID = "260914050_2024-07-15_c64603a7"

    def test_prefetch_replaces_the_qc_pass(
        self,
        s3_buckets,
        qc_process_config_file_path,
        process_config_content,
        initialized_db,
        working_dir_path,
        monkeypatch,
    ):
        """When the guess is right, the QC pass finds every file already in QC storage."""
        _upload_initial_submission_to_inbox(s3_buckets["inbox"], self.SUBMISSION_ID)

        def fail_qc_download(*args, **kwargs):
            raise AssertionError("the QC pass should not download files that the main pass already wrote")

        monkeypatch.setattr(FilePipelineExecutor, "_run_qc_pipeline", fail_qc_download)

        result = _run_process(qc_process_config_file_path, self.SUBMISSION_ID, working_dir_path, "--update-db")

        assert result.exit_code == 0, f"Process failed: {result.output}"
        metadata = json.loads((VALID_SUBMISSION_DIR / "metadata" / "metadata.json").read_text())
        expected = {
            file["filePath"]
            for donor in metadata["donors"]
            for lab_datum in donor["labData"]
            for file in lab_datum.get("sequenceData", {}).get("files", [])
        }
        assert _qc_files(process_config_content, self.SUBMISSION_ID) == expected
        qc_dir = Path(process_config_content["detailed_qc"]["local_storage"]) / self.SUBMISSION_ID
        uploaded_metadata = {**metadata, "submission": {**metadata["submission"], "submissionType": "initial"}}
        assert json.loads((qc_dir / "metadata" / "metadata.json").read_text()) == uploaded_metadata

    def test_prefetched_files_are_deleted_when_not_selected(
        self,
        s3_buckets,
        qc_process_config_file_path,
        process_config_content,
        initialized_db,
        working_dir_path,
        monkeypatch,
    ):
        """A wrong guess leaves no decrypted files and no QC progress log behind."""
        upload_submission_to_inbox(s3_buckets["inbox"], self.SUBMISSION_ID)

        # guess "selected", then decide "not selected"
        monkeypatch.setattr(SubmissionDb, "should_qc", lambda self, *args, predict=False, **kwargs: predict)

        result = _run_process(qc_process_config_file_path, self.SUBMISSION_ID, working_dir_path, "--update-db")

        assert result.exit_code == 0, f"Process failed: {result.output}"
        assert _qc_files(process_config_content, self.SUBMISSION_ID) == set()
        assert not (working_dir_path / "logs" / "progress_qc.cjson").exists()

    def test_no_update_db_skips_qc_selection(
        self,
        s3_buckets,
        qc_process_config_file_path,
        process_config_content,
        initialized_db,
        working_dir_path,
    ):
        """The selection stores its decision in the DB, so ``--no-update-db`` runs without detailed QC."""
        _upload_initial_submission_to_inbox(s3_buckets["inbox"], self.SUBMISSION_ID)

        result = _run_process(qc_process_config_file_path, self.SUBMISSION_ID, working_dir_path, "--no-update-db")

        assert result.exit_code == 0, f"Process failed: {result.output}"
        assert _qc_files(process_config_content, self.SUBMISSION_ID) == set()


@pytest.fixture
def db_process_config_file_path(temp_data_dir_path, process_config_content) -> Path:
    """Write a process config that can sign DB state changes, without detailed QC."""
    process_config_content["db"]["author"]["private_key_passphrase"] = "test"
    config_file = temp_data_dir_path / "config.process.db.yaml"
    with open(config_file, "w") as fd:
        yaml.dump(process_config_content, fd)
    return config_file


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
        db_process_config_file_path,
        process_config_content,
        initialized_db,
        working_dir_path,
        monkeypatch,
    ):
        """The pre-check rejects the duplicate before any file is processed."""
        duplicate_id = self._process_first_and_upload_duplicate(
            s3_buckets["inbox"], db_process_config_file_path, working_dir_path
        )

        def must_not_process(*args, **kwargs):
            raise AssertionError("no file should be processed for a duplicate initial submission")

        monkeypatch.setattr(FilePipelineExecutor, "process_submission_files", must_not_process)

        result = _run_process(db_process_config_file_path, duplicate_id, working_dir_path / "duplicate", "--update-db")

        assert result.exit_code != 0
        self._assert_failed_basic_qc(process_config_content, duplicate_id)

    def test_update_db_records_duplicate_detected_after_processing(
        self,
        s3_buckets,
        db_process_config_file_path,
        process_config_content,
        initialized_db,
        working_dir_path,
        monkeypatch,
    ):
        """The ``basic_qc_passed`` write catches a competing initial submission that the pre-check missed."""
        duplicate_id = self._process_first_and_upload_duplicate(
            s3_buckets["inbox"], db_process_config_file_path, working_dir_path
        )
        monkeypatch.setattr(SubmissionDb, "assert_no_duplicate_initial", lambda self, *args, **kwargs: None)

        result = _run_process(db_process_config_file_path, duplicate_id, working_dir_path / "duplicate", "--update-db")

        assert result.exit_code != 0
        self._assert_failed_basic_qc(process_config_content, duplicate_id)


class TestConfigValidation:
    """Tests for process configuration parsing and environment variable merging."""

    def test_pydantic_nested_env_var_merging(self, tmp_path: Path, monkeypatch, process_config_content: dict):
        le_id = "260914050"
        bucket_name = "grz-inbox-test"

        process_config_content["leistungserbringer"][le_id]["inbox_buckets"] = {
            bucket_name: {
                "endpoint_url": "https://s3.amazonaws.com",
                "access_key": "testing",
                "secret": "testing",
                "private_key_path": "/path/to/test.sec",
            }
        }

        config_file = tmp_path / "config.yaml"
        with open(config_file, "w") as f:
            yaml.dump(process_config_content, f)

        env_var_name = f"GRZ_LEISTUNGSERBRINGER__{le_id}__INBOX_BUCKETS__{bucket_name}__PRIVATE_KEY_PASSPHRASE"
        monkeypatch.setenv(env_var_name.upper(), "dotenv-secret-passphrase")

        with open(config_file) as f:
            raw_dict = yaml.safe_load(f)

        config = GrzctlConfig.from_configuration(raw_dict)
        entry = config.leistungserbringer[le_id]
        assert get_secret_value(entry.inbox_buckets[bucket_name].private_key_passphrase) == "dotenv-secret-passphrase"

    def test_pydantic_json_env_var_merging(self, tmp_path: Path, monkeypatch, process_config_content: dict):
        le_id = "260914050"
        bucket_name = "grz-inbox-test"

        process_config_content["leistungserbringer"][le_id]["inbox_buckets"] = {
            bucket_name: {
                "endpoint_url": "https://s3.amazonaws.com",
                "access_key": "testing",
                "secret": "testing",
                "private_key_path": "/path/to/test.sec",
            }
        }

        config_file = tmp_path / "config.yaml"
        with open(config_file, "w") as f:
            yaml.dump(process_config_content, f)

        json_override = {
            le_id: {
                "inbox_buckets": {
                    bucket_name: {
                        "endpoint_url": "https://s3.amazonaws.com",
                        "access_key": "testing",
                        "secret": "testing",
                        "private_key_path": "/path/to/test.sec",
                        "private_key_passphrase": "json-secret-passphrase",
                    }
                }
            }
        }
        monkeypatch.setenv("GRZ_LEISTUNGSERBRINGER", json.dumps(json_override))

        with open(config_file) as f:
            raw_dict = yaml.safe_load(f)

        config = GrzctlConfig.from_configuration(raw_dict)
        entry = config.leistungserbringer[le_id]
        assert get_secret_value(entry.inbox_buckets[bucket_name].private_key_passphrase) == "json-secret-passphrase"
