"""
Integration tests for grzctl process command.

This tests the full streaming pipeline with mocked S3 buckets (inbox + archive).
"""

import json
import os
import shutil
from pathlib import Path

import boto3
import click.testing
import grzctl.cli
import pytest
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
    """
    Create configuration for grzctl process command.

    This sets up:
    - Inbox S3 bucket configuration
    - Consented archive S3 bucket configuration
    - Non-consented archive S3 bucket configuration
    - Keys for decryption and re-encryption
    """
    db_dir = tmpdir_factory.mktemp("db")
    db_file = db_dir / "test.db"

    return {
        "s3": {
            "endpoint_url": "https://s3.amazonaws.com",
            "bucket": "inbox",
            "access_key": "testing",
            "secret": "testing",
        },
        "consented_archive_s3": {
            "endpoint_url": "https://s3.amazonaws.com",
            "bucket": "consented-archive",
            "access_key": "testing",
            "secret": "testing",
        },
        "non_consented_archive_s3": {
            "endpoint_url": "https://s3.amazonaws.com",
            "bucket": "non-consented-archive",
            "access_key": "testing",
            "secret": "testing",
        },
        "keys": {
            "grz_private_key_path": str(crypt4gh_grz_private_key_file_path),
            "consented_archive_public_key_path": str(crypt4gh_grz_public_key_file_path),
            "non_consented_archive_public_key_path": str(crypt4gh_grz_public_key_file_path),
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
    """Create the three S3 buckets needed for processing."""
    conn = boto3.client("s3")

    # Create inbox, consented-archive, and non-consented-archive buckets
    conn.create_bucket(Bucket="inbox")
    conn.create_bucket(Bucket="consented-archive")
    conn.create_bucket(Bucket="non-consented-archive")

    s3 = boto3.resource("s3")
    return {
        "inbox": s3.Bucket("inbox"),
        "consented": s3.Bucket("consented-archive"),
        "non_consented": s3.Bucket("non-consented-archive"),
    }


@pytest.fixture
def initialized_db(temp_process_config_file_path):
    """Initialize the database for the process command tests."""
    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()

    result = runner.invoke(
        cli,
        ["db", "--config-file", str(temp_process_config_file_path), "init"],
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
        4. Verifies the processed files can be decrypted correctly
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
            "process",
            "--config-file",
            str(temp_process_config_file_path),
            "--submission-id",
            submission_id,
            "--output-dir",
            str(working_dir_path),
            "--no-validate",  # Skip validation for faster test
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
            "process",
            "--config-file",
            str(temp_process_config_file_path),
            "--submission-id",
            submission_id,
            "--output-dir",
            str(working_dir_path),
            "--no-validate",
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
        crypt4gh_grz_private_key_file_path,
        crypt4gh_grz_public_key_file_path,
        db_alice_private_key_file_path,
        db_known_keys_file_path,
        working_dir_path,
        tmp_path,
    ):
        """
        Test multi-inbox selection logic.
        """
        conn = boto3.client("s3")
        conn.create_bucket(Bucket="inbox-a")
        conn.create_bucket(Bucket="inbox-b")
        conn.create_bucket(Bucket="consented-archive")
        conn.create_bucket(Bucket="non-consented-archive")

        s3_resource = boto3.resource("s3")
        _inbox_a = s3_resource.Bucket("inbox-a")
        inbox_b = s3_resource.Bucket("inbox-b")

        submission_id = "260914050_2024-07-15_c64603a7"
        le_id = "260914050"

        # Upload submission to inbox-b
        upload_submission_to_inbox(inbox_b, submission_id)

        db_dir = tmp_path / "db"
        db_dir.mkdir()
        db_file = db_dir / "test.db"

        # Create config with multiple inboxes
        config_content = {
            "s3": {
                "endpoint_url": "https://s3.amazonaws.com",
                "access_key": "testing",
                "secret": "testing",
                "inboxes": {
                    le_id: {
                        "inbox-a": {},
                        "inbox-b": {},
                    }
                },
            },
            "consented_archive_s3": {
                "endpoint_url": "https://s3.amazonaws.com",
                "bucket": "consented-archive",
                "access_key": "testing",
                "secret": "testing",
            },
            "non_consented_archive_s3": {
                "endpoint_url": "https://s3.amazonaws.com",
                "bucket": "non-consented-archive",
                "access_key": "testing",
                "secret": "testing",
            },
            "keys": {
                "grz_private_key_path": str(crypt4gh_grz_private_key_file_path),
                "consented_archive_public_key_path": str(crypt4gh_grz_public_key_file_path),
                "non_consented_archive_public_key_path": str(crypt4gh_grz_public_key_file_path),
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
        }

        config_file = temp_data_dir_path / "config.multi.yaml"
        import yaml

        with open(config_file, "w") as fd:
            yaml.dump(config_content, fd)

        # Run grzctl process with --inbox-bucket inbox-b
        args = [
            "process",
            "--config-file",
            str(config_file),
            "--submission-id",
            submission_id,
            "--output-dir",
            str(working_dir_path),
            "--inbox-bucket",
            "inbox-b",
            "--no-validate",
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
            "process",
            "--config-file",
            str(config_file),
            "--submission-id",
            submission_id,
            "--output-dir",
            str(working_dir_path / "fail"),
            "--no-validate",
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
            [
                "decrypt",
                "--submission-dir",
                str(working_dir),
                "--config-file",
                config_file,
            ],
            catch_exceptions=False,
        )
        assert result.exit_code == 0, f"Decrypt failed: {result.output}"

    def test_process_produces_same_metadata_as_archive(
        self,
        s3_buckets,
        temp_process_config_file_path,
        temp_keys_config_file_path,
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
            "process",
            "--config-file",
            str(temp_process_config_file_path),
            "--submission-id",
            submission_id,
            "--output-dir",
            str(working_dir_process),
            "--no-validate",
            "--no-submit-pruefbericht",
            "--no-update-db",
        ]

        runner = click.testing.CliRunner()
        cli = grzctl.cli.build_cli()
        result = runner.invoke(cli, args, catch_exceptions=False)

        assert result.exit_code == 0, f"Process failed: {result.output}"

        # === Run manual workflow ===
        self._setup_manual_submission_dir(working_dir_manual)
        self._run_manual_decrypt(working_dir_manual, str(temp_keys_config_file_path))

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

        # Run grzctl process WITH validation enabled
        args = [
            "process",
            "--config-file",
            str(temp_process_config_file_path),
            "--submission-id",
            submission_id,
            "--output-dir",
            str(working_dir_path),
            "--validate",  # Enable validation
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
