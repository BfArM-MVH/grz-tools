"""S3MultipartUploader stores every non-empty stream as a multipart upload, and an empty one with a PUT."""

import io

import boto3
import pytest
from grz_common.pipeline.components import ReadStream
from grz_common.pipeline.components.s3 import S3MultipartUploader
from moto import mock_aws

BUCKET = "archive"
KEY = "subA/files/f.c4gh"


@pytest.fixture
def s3():
    with mock_aws():
        # "us-east-1" is only a placeholder for the moto mock; it is S3's default region, so
        # create_bucket() works without a LocationConstraint.
        client = boto3.client("s3", region_name="us-east-1")
        client.create_bucket(Bucket=BUCKET)
        yield client


def test_empty_stream_creates_an_empty_object(s3):
    """A multipart upload needs at least one part, so an empty stream must still create an object."""
    ReadStream(io.BytesIO(b"")) >> S3MultipartUploader(s3, BUCKET, KEY)

    stored = s3.get_object(Bucket=BUCKET, Key=KEY)
    assert stored["Body"].read() == b""


def test_explicit_empty_write_creates_an_empty_object(s3):
    """An empty write() must not start a multipart upload that close() cannot complete."""
    uploader = S3MultipartUploader(s3, BUCKET, KEY)
    uploader.write(b"")
    uploader.close()

    assert s3.get_object(Bucket=BUCKET, Key=KEY)["Body"].read() == b""


@pytest.mark.parametrize("size", [1, 1024, 9 * 1024 * 1024])
def test_non_empty_stream_uses_multipart(s3, size):
    """Any non-empty stream, however small, is stored as a multipart upload."""
    data = b"x" * size

    ReadStream(io.BytesIO(data)) >> S3MultipartUploader(s3, BUCKET, KEY)

    stored = s3.get_object(Bucket=BUCKET, Key=KEY)
    assert stored["Body"].read() == data
    assert stored["ETag"].strip('"').endswith("-1"), "a one-part multipart upload has an ETag ending in -1"
