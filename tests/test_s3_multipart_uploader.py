"""S3MultipartUploader stores every non-empty stream as a multipart upload, and an empty one with a PUT."""

import io

import boto3
import botocore.client
import botocore.exceptions
import pytest
from grz_common.constants import MULTIPART_MIN_PART_SIZE
from grz_common.exceptions import UploadError
from grz_common.pipeline.components import ReadStream, UploadIntegrityError
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


@pytest.fixture
def wrong_etag_from(monkeypatch):
    """Let an S3 operation answer with an ETag of the wrong content, as a corrupted transfer would."""

    def install(operation: str) -> None:
        original_call = botocore.client.BaseClient._make_api_call

        def answer_with_a_wrong_etag(self, operation_name, kwargs):
            response = original_call(self, operation_name, kwargs)
            if operation_name == operation:
                response["ETag"] = '"{}"'.format("0" * 32)
            return response

        monkeypatch.setattr(botocore.client.BaseClient, "_make_api_call", answer_with_a_wrong_etag)

    return install


@pytest.fixture
def failing_upload_part(monkeypatch):
    """Let every UploadPart fail, as an S3 outage would."""
    original_call = botocore.client.BaseClient._make_api_call

    def fail_the_part(self, operation_name, kwargs):
        if operation_name == "UploadPart":
            error = {"Error": {"Code": "ServiceUnavailable", "Message": "simulated outage"}}
            raise botocore.exceptions.ClientError(error, operation_name)
        return original_call(self, operation_name, kwargs)

    monkeypatch.setattr(botocore.client.BaseClient, "_make_api_call", fail_the_part)


def _pending_uploads(s3) -> list:
    """The multipart uploads the bucket still holds parts for."""
    return s3.list_multipart_uploads(Bucket=BUCKET).get("Uploads", [])


def _assert_nothing_was_stored(s3) -> None:
    assert _pending_uploads(s3) == [], "a failed upload must not leave its parts staged"
    assert "Contents" not in s3.list_objects_v2(Bucket=BUCKET, Prefix=KEY), "a failed upload must store no object"


def test_a_part_stored_differently_fails_the_upload(s3, wrong_etag_from):
    """A part whose ETag does not match what was sent fails as a data integrity error."""
    wrong_etag_from("UploadPart")

    with pytest.raises(UploadIntegrityError):
        ReadStream(io.BytesIO(b"x" * 1024)) >> S3MultipartUploader(s3, BUCKET, KEY)

    _assert_nothing_was_stored(s3)


def test_an_object_assembled_differently_fails_the_upload(s3, wrong_etag_from):
    """A completed object whose ETag does not match the parts fails as a data integrity error."""
    wrong_etag_from("CompleteMultipartUpload")

    with pytest.raises(UploadIntegrityError):
        ReadStream(io.BytesIO(b"x" * 1024)) >> S3MultipartUploader(s3, BUCKET, KEY)

    _assert_nothing_was_stored(s3)


def test_a_part_that_cannot_be_uploaded_fails_the_upload(s3, failing_upload_part):
    """An S3 error while uploading a part fails as an upload error, not as an integrity error."""
    with pytest.raises(UploadError) as exception:
        ReadStream(io.BytesIO(b"x" * 1024)) >> S3MultipartUploader(s3, BUCKET, KEY)

    assert not isinstance(exception.value, UploadIntegrityError)
    _assert_nothing_was_stored(s3)


def test_a_stream_of_whole_parts_uploads_no_empty_last_part(s3):
    """A stream that is an exact multiple of the part size is stored in exactly that many parts."""
    data = b"x" * (2 * MULTIPART_MIN_PART_SIZE)

    ReadStream(io.BytesIO(data)) >> S3MultipartUploader(s3, BUCKET, KEY, part_size=MULTIPART_MIN_PART_SIZE)

    stored = s3.get_object(Bucket=BUCKET, Key=KEY)
    assert stored["Body"].read() == data
    assert stored["ETag"].strip('"').endswith("-2"), "two whole parts, and no empty third one"


def test_an_empty_object_stored_differently_fails_the_upload(s3, wrong_etag_from):
    """The PUT that stores an empty stream is checked the same way, and leaves nothing behind either."""
    wrong_etag_from("PutObject")

    with pytest.raises(UploadIntegrityError):
        ReadStream(io.BytesIO(b"")) >> S3MultipartUploader(s3, BUCKET, KEY)

    _assert_nothing_was_stored(s3)
