"""Fixtures for the tests."""

import json
import os
from datetime import datetime
from importlib.metadata import version
from os import PathLike
from pathlib import Path
from shutil import copyfile

import boto3
import grz_cli.models.config
import grz_common.models.s3
import grzctl.models.config
import pytest
import yaml
from grz_common.utils.crypt import Crypt4GH
from grz_common.workers.submission import EncryptedSubmission, SubmissionMetadata
from grz_db import testing as grz_db_testing
from grz_db.models.submission import SubmissionDb
from moto import mock_aws

# These names are bound by assignment rather than imported directly, because fixtures below take
# ``db_test_connection`` as a parameter name, which ruff would otherwise read as redefining an
# imported name.
db_backend = grz_db_testing.db_backend
db_test_connection = grz_db_testing.db_test_connection
migrated_db_connection = grz_db_testing.migrated_db_connection
postgresql_proc = grz_db_testing.postgresql_proc
_migrated_postgresql_template = grz_db_testing._migrated_postgresql_template
_migrated_sqlite_template = grz_db_testing._migrated_sqlite_template
_pg_test_dbnames = grz_db_testing._pg_test_dbnames


@pytest.fixture(scope="session", autouse=True)
def mock_home(tmp_path_factory):
    # Create a temporary directory for the session
    temp_home = tmp_path_factory.mktemp("fake_home")

    # Set the environment variable
    os.environ["HOME"] = str(temp_home)
    os.environ["USERPROFILE"] = str(temp_home)  # For Windows compatibility

    return temp_home


config_path = "tests/mock_files/mock_config.yaml"
small_file_input_path = "tests/mock_files/mock_small_input_file.bed"
metadata_path = "tests/mock_files/submissions/valid_submission/metadata/metadata.json"

crypt4gh_grz_private_key_file = "tests/mock_files/grz_mock_private_key.sec"
crypt4gh_grz_public_key_file = "tests/mock_files/grz_mock_public_key.pub"
crypt4gh_submitter_private_key_file = "tests/mock_files/submitter_mock_private_key.sec"
crypt4gh_submitter_public_key_file = "tests/mock_files/submitter_mock_public_key.pub"
db_alice_private_key_file = "tests/mock_files/db/alice_mock_private_key.sec"
db_known_keys_file = "tests/mock_files/db/known_keys"


@pytest.fixture()
def crypt4gh_grz_private_key_file_path():
    return Path(crypt4gh_grz_private_key_file)


@pytest.fixture()
def crypt4gh_grz_public_key_file_path():
    return Path(crypt4gh_grz_public_key_file)


@pytest.fixture()
def crypt4gh_submitter_private_key_file_path():
    return Path(crypt4gh_submitter_private_key_file)


@pytest.fixture()
def crypt4gh_submitter_public_key_file_path():
    return Path(crypt4gh_submitter_public_key_file)


@pytest.fixture()
def db_alice_private_key_file_path():
    return Path(db_alice_private_key_file)


@pytest.fixture()
def db_known_keys_file_path():
    return Path(db_known_keys_file)


@pytest.fixture()
def initiated_db_test_connection(db_test_connection):
    """The database behind :func:`db_config_content`, brought to the latest schema.

    This is not the shared ``migrated_db_connection`` fixture: that one hands back a separate
    database, whereas the tests here need the specific database their config file already points at.
    """
    submission_db = SubmissionDb(db_test_connection, author=None)
    submission_db.initialize_schema()
    return db_test_connection


@pytest.fixture(scope="session")
def temp_data_dir(tmpdir_factory: pytest.TempdirFactory):
    """Create temporary folder for the session"""
    datadir = tmpdir_factory.mktemp("data")
    return datadir


@pytest.fixture
def temp_data_dir_path(temp_data_dir) -> Path:
    return Path(temp_data_dir.strpath)


def copy_file_to_tempdir(input_path: str | PathLike, datadir: str | PathLike) -> Path:
    filename = os.path.basename(input_path)
    target = Path(datadir) / filename

    copyfile(input_path, target)

    return target


@pytest.fixture
def temp_small_file_path(temp_data_dir_path) -> Path:
    return copy_file_to_tempdir(small_file_input_path, temp_data_dir_path)


@pytest.fixture()
def temp_small_file_md5sum():
    return "710781ec9efd25b87bfbf8d6cf4030e9"


@pytest.fixture()
def temp_small_file_sha256sum():
    return "78858035d88f0c66d27984789ddd8fa8a8fc633cf7689ac2b4b1e2e7b37ee3be"


@pytest.fixture
def remote_bucket_with_version(remote_bucket):
    """Mock S3 bucket with version.json file at root."""
    current_version = version("grz-cli")

    version_content = {
        "schema_version": 1,
        "grzcli_version": [
            {
                "minimal_version": current_version,
                "recommended_version": current_version,
                "max_version": current_version,
                "enforced_from": datetime.now().isoformat(),
            }
        ],
        "metadata_version": [
            {
                "minimal_version": "1.3.0",
                "recommended_version": "1.3.0",
                "max_version": None,
                "enforced_from": datetime.now().isoformat(),
            }
        ],
    }

    remote_bucket.put_object(
        Key="version.json",
        Body=json.dumps(version_content),
    )

    return remote_bucket


def generate_fastq(file_path: str | PathLike, target_size: int) -> int:
    """Create a FASTQ file of at least *target_size* bytes.

    The reads repeat rather than being sampled per base: the tests using this only round-trip the
    file through S3 and compare checksums, so the content never has to look like real sequencing
    data. The size does matter: it is what pushes the transfer over the multipart threshold.

    :param file_path: Path to the FASTQ file.
    :param target_size: Target size in bytes.
    :return: Actual bytes written
    """
    bases_per_read = 100
    record = f"@SEQ_ID\n{'ACGT' * (bases_per_read // 4)}\n+\n{'I' * bases_per_read}\n"
    repeats = -(-target_size // len(record))  # ceiling division

    with open(file_path, "w") as fastq_file:
        return fastq_file.write(record * repeats)


@pytest.fixture
def temp_fastq_file_path(temp_data_dir_path) -> Path:
    file_name = "5M.fastq"
    temp_fastq_gz_path = temp_data_dir_path / file_name
    target_size = 1024 * 1024 * 6  # create 5MB file, multiupload limit is 5MB

    generate_fastq(temp_fastq_gz_path, target_size)

    return temp_fastq_gz_path


@pytest.fixture
def temp_fastq_file_md5sum(temp_fastq_file_path):
    import hashlib

    with open(temp_fastq_file_path, "rb") as f:
        file_hash = hashlib.md5()
        while chunk := f.read(8192):
            file_hash.update(chunk)

    return file_hash.hexdigest()


@pytest.fixture
def temp_fastq_file_sha256sum(temp_fastq_file_path):
    import hashlib

    with open(temp_fastq_file_path, "rb") as f:
        file_hash = hashlib.sha256()
        while chunk := f.read(8192):
            file_hash.update(chunk)

    return file_hash.hexdigest()


# @pytest.fixture
# def crypt4gh_private_key_file(temp_data_dir_path):
#     return copy_file_to_tempdir(private_key_path, temp_data_dir_path)
#
#
# @pytest.fixture
# def crypt4gh_public_key_file(temp_data_dir_path):
#     return copy_file_to_tempdir(public_key_path, temp_data_dir_path)


@pytest.fixture
def temp_metadata_file_path(temp_data_dir_path) -> Path:
    """Metadata naming a file that is never created: validation reads paths; it does not open them."""
    with open(metadata_path) as fd:
        metadata = json.load(fd)

    metadata["donors"][0]["labData"][0]["sequenceData"]["files"][0]["filePath"] = "temp_large_input_file.bed"

    metadata_file_path = temp_data_dir_path / "metadata.json"
    with open(metadata_file_path, "w") as fd:
        json.dump(metadata, fd)

    return metadata_file_path


@pytest.fixture
def keys_config_content(
    crypt4gh_grz_public_key_file_path,
    crypt4gh_grz_private_key_file_path,
    crypt4gh_submitter_private_key_file_path,
):
    return {
        "keys": {
            "grz_public_key_path": str(crypt4gh_grz_public_key_file_path),
            "grz_private_key_path": str(crypt4gh_grz_private_key_file_path),
            "submitter_private_key_path": str(crypt4gh_submitter_private_key_file_path),
        }
    }


@pytest.fixture
def s3_config_content():
    return {
        "s3": {
            "endpoint_url": "https://s3.amazonaws.com",
            "bucket": "testing",
            "access_key": "testing",
            "secret": "testing",
        }
    }


def _db_config_content(private_key_path, known_keys_path, database_url):
    return {
        "db": {
            "database_url": database_url,
            "author": {"name": "Alice", "private_key_path": str(private_key_path)},
            "known_public_keys": str(known_keys_path),
        }
    }


@pytest.fixture
def db_config_content(
    db_alice_private_key_file_path,
    db_known_keys_file_path,
    db_test_connection,
):
    """Config for a database with no schema, for tests that run ``db init`` themselves."""
    return _db_config_content(db_alice_private_key_file_path, db_known_keys_file_path, db_test_connection)


@pytest.fixture
def migrated_db_config_content(
    db_alice_private_key_file_path,
    db_known_keys_file_path,
    migrated_db_connection,
):
    """Config for a database already on the latest schema."""
    return _db_config_content(db_alice_private_key_file_path, db_known_keys_file_path, migrated_db_connection)


@pytest.fixture
def pruefbericht_config_content():
    return {
        "pruefbericht": {
            "authorization_url": "https://bfarm.localhost/token",
            "api_base_url": "https://bfarm.localhost/api/",
            "client_id": "pytest",
            "client_secret": "pysecret",
        }
    }


@pytest.fixture
def identifiers_config_content():
    return {
        "identifiers": {
            "grz": "GRZK00007",
            "le": "260914050",
        }
    }


@pytest.fixture
def s3_config_model(s3_config_content):
    return grz_common.models.s3.S3ConfigModel(**s3_config_content)


@pytest.fixture
def encrypt_config_model(keys_config_content):
    return grz_cli.models.config.EncryptConfig(**keys_config_content)


@pytest.fixture
def db_config_model(db_config_content):
    return grzctl.models.config.DbConfig(**db_config_content)


@pytest.fixture
def migrated_db_config_model(migrated_db_config_content):
    return grzctl.models.config.DbConfig(**migrated_db_config_content)


@pytest.fixture
def identifiers_config_model(identifiers_config_content):
    return grz_cli.models.config.ValidateConfig(**identifiers_config_content)


@pytest.fixture
def pruefbericht_config_model(pruefbericht_config_content):
    return grzctl.models.config.PruefberichtConfig(**pruefbericht_config_content)


@pytest.fixture
def temp_s3_db_config_file_path(temp_data_dir_path, s3_config_model, db_config_model) -> Path:
    config_file = temp_data_dir_path / "config.db_s3.yaml"

    combined = {
        **s3_config_model.model_dump(mode="json", exclude_none=True, exclude_unset=True, exclude_defaults=True),
        **db_config_model.model_dump(mode="json", exclude_none=True, exclude_unset=True, exclude_defaults=True),
    }

    with open(config_file, "w") as fd:
        yaml.safe_dump(combined, fd, sort_keys=False)

    return config_file


@pytest.fixture
def temp_s3_config_file_path(temp_data_dir_path, s3_config_model) -> Path:
    config_file = temp_data_dir_path / "config.s3.yaml"
    with open(config_file, "w") as fd:
        s3_config_model.to_yaml(fd)
    return config_file


@pytest.fixture
def temp_db_config_file_path(temp_data_dir_path, db_config_model) -> Path:
    """Config file for a database with no schema, for tests that run ``db init`` themselves."""
    config_file = temp_data_dir_path / "config.db.yaml"
    with open(config_file, "w") as fd:
        db_config_model.to_yaml(fd)
    return config_file


@pytest.fixture
def temp_migrated_db_config_file_path(temp_data_dir_path, migrated_db_config_model) -> Path:
    """Config file for a database already on the latest schema."""
    config_file = temp_data_dir_path / "config.migrated_db.yaml"
    with open(config_file, "w") as fd:
        migrated_db_config_model.to_yaml(fd)
    return config_file


@pytest.fixture
def temp_keys_config_file_path(temp_data_dir_path, encrypt_config_model) -> Path:
    config_file = temp_data_dir_path / "config.keys.yaml"
    with open(config_file, "w") as fd:
        encrypt_config_model.to_yaml(fd)
    return config_file


@pytest.fixture
def temp_identifiers_config_file_path(temp_data_dir_path, identifiers_config_model) -> Path:
    config_file = temp_data_dir_path / "config.ids.yaml"
    with open(config_file, "w") as fd:
        identifiers_config_model.to_yaml(fd)
    return config_file


@pytest.fixture
def temp_pruefbericht_config_file_path(temp_data_dir_path, pruefbericht_config_model) -> Path:
    config_file = temp_data_dir_path / "config.pruefbericht.yaml"
    with open(config_file, "w") as fd:
        pruefbericht_config_model.to_yaml(fd)
    return config_file


@pytest.fixture
def crypt4gh_grz_public_keys(crypt4gh_grz_public_key_file_path, crypt4gh_submitter_private_key_file_path):
    keys = Crypt4GH.prepare_c4gh_keys(
        recipient_key_file_path=crypt4gh_grz_public_key_file_path,
        sender_private_key=crypt4gh_submitter_private_key_file_path,
    )
    return keys


@pytest.fixture
def aws_credentials(s3_config_model):
    """Mocked AWS Credentials for moto."""
    os.environ["AWS_ACCESS_KEY_ID"] = s3_config_model.s3.access_key
    os.environ["AWS_SECRET_ACCESS_KEY"] = s3_config_model.s3.secret
    os.environ["MOTO_ALLOW_NONEXISTENT_REGION"] = "1"
    with mock_aws():
        yield


@pytest.fixture
def boto_s3_client(aws_credentials):
    conn = boto3.client("s3")
    yield conn


@pytest.fixture
def remote_bucket(boto_s3_client, s3_config_model):
    # create bucket
    boto_s3_client.create_bucket(Bucket=s3_config_model.s3.bucket)

    return boto3.resource("s3").Bucket(s3_config_model.s3.bucket)


@pytest.fixture
def submission_metadata_dir() -> Path:
    return Path("tests/mock_files/submissions/valid_submission/metadata")


@pytest.fixture
def submission_metadata(submission_metadata_dir) -> SubmissionMetadata:
    return SubmissionMetadata(submission_metadata_dir / "metadata.json")


@pytest.fixture
def encrypted_files_dir() -> Path:
    return Path("tests/mock_files/submissions/valid_submission/encrypted_files")


@pytest.fixture
def encrypted_submission(submission_metadata_dir, encrypted_files_dir) -> EncryptedSubmission:
    submission = EncryptedSubmission(submission_metadata_dir, encrypted_files_dir)
    return submission


@pytest.fixture
def working_dir(tmpdir_factory: pytest.TempdirFactory):
    """Create temporary folder for the session"""
    datadir = tmpdir_factory.mktemp("submission")
    return datadir


@pytest.fixture
def working_dir_path(working_dir) -> Path:
    return Path(working_dir.strpath)
