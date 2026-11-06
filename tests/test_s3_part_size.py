"""calculate_s3_part_size keeps a multipart upload within the part limit, and S3Downloader reports the object size."""

import math

import boto3
import pytest
from grz_common.constants import MULTIPART_DEFAULT_PART_SIZE, MULTIPART_MAX_PARTS, MULTIPART_MIN_PART_SIZE
from grz_common.pipeline.components.s3 import S3Downloader, calculate_s3_part_size
from moto import mock_aws

MiB = 1024 * 1024
PREFERRED = 64 * MiB


@pytest.mark.parametrize("file_size", [None, 0, 10 * MiB, PREFERRED * MULTIPART_MAX_PARTS])
def test_uses_preferred_part_size_when_the_file_fits(file_size):
    assert calculate_s3_part_size(file_size, PREFERRED) == PREFERRED


def test_preferred_part_size_defaults_to_the_default_part_size():
    assert calculate_s3_part_size(10 * MiB) == MULTIPART_DEFAULT_PART_SIZE


def test_part_size_is_never_below_the_minimum():
    assert calculate_s3_part_size(10 * MiB, 1) == MULTIPART_MIN_PART_SIZE


@pytest.mark.parametrize("file_size", [PREFERRED * MULTIPART_MAX_PARTS + 1, 300 * 1024 * MiB])
def test_part_size_grows_to_stay_within_the_part_limit(file_size):
    part_size = calculate_s3_part_size(file_size, PREFERRED)

    assert part_size > PREFERRED
    assert math.ceil(file_size / part_size) <= MULTIPART_MAX_PARTS


@mock_aws
def test_downloader_length_is_the_object_size():
    # "us-east-1" is only a placeholder for the moto mock; it is S3's default region, so
    # create_bucket() works without a LocationConstraint.
    s3 = boto3.client("s3", region_name="us-east-1")
    s3.create_bucket(Bucket="inbox")
    s3.put_object(Bucket="inbox", Key="f.c4gh", Body=b"x" * 1234)

    with S3Downloader(s3, "inbox", "f.c4gh") as downloader:
        assert downloader.length == 1234
