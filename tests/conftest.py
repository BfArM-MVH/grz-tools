"""Fixtures for the tests."""

import json
import os
from os import PathLike
from pathlib import Path
from shutil import copyfile, which

import boto3
import grz_cli.models.config
import grz_common.models.s3
import grzctl.models.config
import numpy as np
import psycopg
import pytest
from grz_common.utils.crypt import Crypt4GH
from grz_common.workers.submission import EncryptedSubmission, SubmissionMetadata
from moto import mock_aws

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


@pytest.fixture(
    params=[
        "sqlite",
        pytest.param(
            "postgresql",
            marks=pytest.mark.skipif(condition=which("pg_config") is None, reason="postgresql not detected"),
        ),
    ]
)
def db_test_connection(request: pytest.FixtureRequest):
    if request.param == "sqlite":
        tmpdir_factory: pytest.TempdirFactory = request.getfixturevalue("tmpdir_factory")
        db_dir = tmpdir_factory.mktemp("db")
        db_file = db_dir / "test.db"
        yield f"sqlite:///{str(db_file)}"
    elif request.param == "postgresql":
        postgresql: psycopg.Connection = request.getfixturevalue("postgresql")
        yield f"postgresql+psycopg://{postgresql.info.user}:@{postgresql.info.host}:{postgresql.info.port}/{postgresql.info.dbname}"


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


def create_large_file(content: str | bytes, output_file: str | PathLike, target_size: int) -> int:
    """
    Write some content repeatedly to a file until some target size is reached.

    :param content: Content of the file. Will be repeatedly written to output_file until `target_size` is reached.
    :param output_file: Path to the output file.
    :param target_size: target size in bytes.
    :return: Actual bytes written
    """
    # Initialize the size of the new file and open it for writing
    current_size = 0
    with open(output_file, "w") as outfile:
        while current_size < target_size:
            bytes_written = outfile.write(content)
            current_size += bytes_written
    return current_size


@pytest.fixture
def temp_large_file_path(temp_data_dir_path) -> Path:
    temp_large_file_path = temp_data_dir_path / "temp_large_input_file.bed"
    target_size = 1024 * 1024 * 6  # create 5MB file, multiupload limit is 5MB

    with open(small_file_input_path) as fd:
        content = fd.read()

    create_large_file(content, temp_large_file_path, target_size)

    return temp_large_file_path


def generate_random_fastq(file_path: str | PathLike, target_size: int) -> int:
    """
    Create a random FASTQ file.

    :param file_path: Path to the FASTQ file.
    :param target_size: Target size in bytes.
    :return: Actual bytes written
    """
    nucleotides = np.array(["A", "T", "C", "G"])
    quality_scores = np.array(list("!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI"))
    bases_per_read = 100  # Length of each read

    with open(file_path, "w") as fastq_file:
        total_bytes_written = 0

        while total_bytes_written < target_size:
            # Generate a random sequence of nucleotides using numpy
            seq = "".join(np.random.choice(nucleotides, bases_per_read))
            qual = "".join(np.random.choice(quality_scores, bases_per_read))

            # FASTQ entry format
            entry = f"@SEQ_ID_{np.random.randint(1, 10**6)}\n{seq}\n+\n{qual}\n"

            actual_bytes_written = fastq_file.write(entry)
            total_bytes_written += actual_bytes_written

    return total_bytes_written


@pytest.fixture
def temp_fastq_file_path(temp_data_dir_path) -> Path:
    file_name = "5M.fastq"
    temp_fastq_gz_path = temp_data_dir_path / file_name
    target_size = 1024 * 1024 * 6  # create 5MB file, multiupload limit is 5MB

    generate_random_fastq(temp_fastq_gz_path, target_size)

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
def temp_metadata_file_path(temp_data_dir_path, temp_large_file_path) -> Path:
    with open(metadata_path) as fd:
        metadata = json.load(fd)

    # insert large file
    metadata["donors"][0]["labData"][0]["sequenceData"]["files"][0]["filePath"] = temp_large_file_path.name

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


@pytest.fixture
def db_config_content(
    db_alice_private_key_file_path,
    db_known_keys_file_path,
    db_test_connection,
):
    return {
        "db": {
            "database_url": db_test_connection,
            "author": {"name": "Alice", "private_key_path": str(db_alice_private_key_file_path)},
            "known_public_keys": str(db_known_keys_file_path),
        }
    }


@pytest.fixture
def pruefbericht_config_content():
    return {
        "pruefbericht": {
            "authorization_url": "https://bfarm.localhost/token",
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
def identifiers_config_model(identifiers_config_content):
    return grz_cli.models.config.ValidateConfig(**identifiers_config_content)


@pytest.fixture
def pruefbericht_config_model(pruefbericht_config_content):
    return grzctl.models.config.PruefberichtConfig(**pruefbericht_config_content)


@pytest.fixture
def temp_s3_config_file_path(temp_data_dir_path, s3_config_model) -> Path:
    config_file = temp_data_dir_path / "config.s3.yaml"
    with open(config_file, "w") as fd:
        s3_config_model.to_yaml(fd)
    return config_file


@pytest.fixture
def temp_db_config_file_path(temp_data_dir_path, db_config_model) -> Path:
    config_file = temp_data_dir_path / "config.db.yaml"
    with open(config_file, "w") as fd:
        db_config_model.to_yaml(fd)
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
