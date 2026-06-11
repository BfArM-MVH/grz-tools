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

import boto3
import pytest
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
    s3 = boto3.client("s3", region_name="us-east-1")
    s3.create_bucket(Bucket=BUCKET)

    source = ReadStream(io.BytesIO(b"x" * (1024 * 1024)))  # 1 MiB -> starts a multipart upload
    chain = source | Tee(_FailOnCloseObserver())
    uploader = S3MultipartUploader(s3, BUCKET, KEY)

    with pytest.raises(RuntimeError):
        chain >> uploader

    # A failure must not commit a (partial) object.
    objects = [o["Key"] for o in s3.list_objects_v2(Bucket=BUCKET).get("Contents", [])]
    assert objects == [], f"failure must not commit an object, found: {objects}"

    # A failure must not leave an incomplete multipart upload dangling on the server.
    uploads = [u["Key"] for u in s3.list_multipart_uploads(Bucket=BUCKET).get("Uploads", [])]
    assert uploads == [], f"failure must abort the upload, dangling multipart upload(s): {uploads}"
