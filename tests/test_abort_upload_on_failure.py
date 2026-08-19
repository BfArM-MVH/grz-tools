"""
Regression test for failure-path upload handling in `Pipeable.__rshift__`.

When the pipeline fails (here: a validator that only fails at ``close()``, mirroring the
deferred FASTQ/BAM/checksum validators), the S3 multipart upload must be **aborted**:

- no partial object may be committed to the archive, and
- no incomplete multipart upload may be left dangling on the server.

This is expected to FAIL until `__rshift__` aborts the sink on the failure path
(currently `self.close()` raises in the `finally` and `other.close()` is skipped, so the
upload is neither completed nor aborted).
"""

import io
import logging

import boto3
import pytest
from botocore.exceptions import ClientError
from grz_common.pipeline.components import Observer, ReadStream, Tee
from grz_common.pipeline.components.s3 import S3MultipartUploader
from moto import mock_aws

BUCKET = "archive"
KEY = "subA/files/f.c4gh"


class _FailOnCloseObserver(Observer):
    """Stand-in for a validator that accepts all data but only fails at close()."""

    def observe(self, chunk: bytes) -> None:
        pass

    def close(self) -> None:
        if not self.closed:
            try:
                raise RuntimeError("simulated validation failure at close()")
            finally:
                super().close()


@mock_aws
def test_failed_pipeline_aborts_upload():
    # "us-east-1" is just a placeholder for the moto mock -- no real AWS or region is involved.
    # We pick it specifically because it's S3's default region, so create_bucket() below works
    # without an extra CreateBucketConfiguration/LocationConstraint argument.
    s3 = boto3.client("s3", region_name="us-east-1")
    s3.create_bucket(Bucket=BUCKET)

    # 9 MiB is larger than the default 8 MiB part size, so at least one part is actually
    # uploaded before the failure -- not just the create_multipart_upload call.
    source = ReadStream(io.BytesIO(b"x" * (9 * 1024 * 1024)))
    chain = source | Tee(_FailOnCloseObserver())
    uploader = S3MultipartUploader(s3, BUCKET, KEY)

    with pytest.raises(RuntimeError):
        chain >> uploader

    # A failure must not commit a (partial) object. S3 omits "Contents" entirely for an empty
    # bucket, so assert the key is absent rather than hiding it behind a default empty list.
    objects = s3.list_objects_v2(Bucket=BUCKET)
    assert "Contents" not in objects, f"failure must not commit an object, found: {objects.get('Contents')}"

    # A failure must not leave an incomplete multipart upload dangling on the server.
    uploads = s3.list_multipart_uploads(Bucket=BUCKET)
    assert "Uploads" not in uploads, f"failure must abort the upload, dangling: {uploads.get('Uploads')}"


@mock_aws
def test_failed_pipeline_surfaces_original_error_when_abort_is_denied(caplog):
    """On a bucket where we may write parts but lack s3:AbortMultipartUpload, a failed abort
    must not hide the real error: the original validation failure still has to surface, and the
    denied abort is only logged as a warning. The dangling upload is then left for a bucket
    lifecycle rule or an admin to reap.
    """
    s3 = boto3.client("s3", region_name="us-east-1")  # moto placeholder, see note above
    s3.create_bucket(Bucket=BUCKET)

    source = ReadStream(io.BytesIO(b"x" * (9 * 1024 * 1024)))
    chain = source | Tee(_FailOnCloseObserver())
    uploader = S3MultipartUploader(s3, BUCKET, KEY)

    def _deny_abort(*args, **kwargs):
        raise ClientError(
            {"Error": {"Code": "AccessDenied", "Message": "abort not permitted"}},
            "AbortMultipartUpload",
        )

    uploader.s3.abort_multipart_upload = _deny_abort

    # The real validation failure must still surface -- a denied abort must not mask it.
    with caplog.at_level(logging.WARNING):
        with pytest.raises(RuntimeError, match="simulated validation failure"):
            chain >> uploader

    # The denied abort is reported, not raised.
    assert any("Could not abort multipart upload" in r.message for r in caplog.records)

    # Even when abort is denied, no (partial) object may be committed.
    objects = s3.list_objects_v2(Bucket=BUCKET)
    assert "Contents" not in objects, f"failure must not commit an object, found: {objects.get('Contents')}"
