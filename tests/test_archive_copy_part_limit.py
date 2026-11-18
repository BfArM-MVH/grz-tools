"""The copy into the final archive stays within ``MULTIPART_MAX_PARTS`` parts."""

import boto3
import botocore.client
import grz_common.pipeline.components.s3
from grzctl.processor import _copy_object
from moto import mock_aws

SOURCE = "interrogation"
TARGET = "archive"
KEY = "subA/files/f.c4gh"


@mock_aws
def test_copy_stays_within_part_limit(monkeypatch):
    """An object that boto3's 8 MiB default would copy in too many parts is copied in fewer, larger parts."""
    # with a limit of 2 parts, the 8 MiB default would copy 20 MiB in 3 parts
    monkeypatch.setattr(grz_common.pipeline.components.s3, "MULTIPART_MAX_PARTS", 2)
    s3 = boto3.client("s3", region_name="us-east-1")
    s3.create_bucket(Bucket=SOURCE)
    s3.create_bucket(Bucket=TARGET)
    body = bytes(20 * 1024 * 1024)
    s3.put_object(Bucket=SOURCE, Key=KEY, Body=body)

    part_copies = []
    make_api_call = botocore.client.BaseClient._make_api_call

    def record(client, operation_name, api_params):
        if operation_name == "UploadPartCopy":
            part_copies.append(api_params["CopySourceRange"])
        return make_api_call(client, operation_name, api_params)

    monkeypatch.setattr(botocore.client.BaseClient, "_make_api_call", record)

    _copy_object(s3, SOURCE, TARGET, KEY, preferred_part_size=5 * 1024 * 1024)

    assert len(part_copies) == 2
    assert s3.get_object(Bucket=TARGET, Key=KEY)["Body"].read() == body
