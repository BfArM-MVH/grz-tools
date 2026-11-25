from collections.abc import Iterator
from typing import Any

import boto3
import pytest
from grz_common.models.s3 import S3Options
from grzctl.commands.delete import storage
from grzctl.commands.delete.storage import (
    Holdings,
    StorageDeletionError,
    archive_holdings,
    delete_from_archive,
    delete_from_inbox,
    has_object_lock,
    inbox_holdings,
)
from moto import mock_aws

REGION = "us-east-1"
SUBMISSION_ID = "000000000_2025-01-01_00000000"
OTHER_SUBMISSION_ID = "000000000_2025-02-02_11111111"

ARCHIVED = {
    f"{SUBMISSION_ID}/files/sample_R1.fastq.gz.c4gh": b"x" * 100,
    f"{SUBMISSION_ID}/files/sample_R2.fastq.gz.c4gh": b"y" * 50,
    f"{SUBMISSION_ID}/metadata/metadata.json": b'{"tanG": "..."}',
    f"{SUBMISSION_ID}/logs/progress_encrypt.cjson": b"log",
    f"{SUBMISSION_ID}/logs/progress_upload.cjson": b"log",
}
UPLOADED = {
    f"{SUBMISSION_ID}/files/sample_R1.fastq.gz.c4gh": b"x" * 100,
    f"{SUBMISSION_ID}/metadata/metadata.json": b'{"tanG": "..."}',
    f"{SUBMISSION_ID}/version": b"3.0.3",
}


@pytest.fixture
def s3() -> Iterator[Any]:
    with mock_aws():
        yield boto3.client("s3", region_name=REGION)


def _bucket(s3: Any, name: str, objects: dict[str, bytes], versioned: bool = False) -> S3Options:
    s3.create_bucket(Bucket=name)
    if versioned:
        s3.put_bucket_versioning(Bucket=name, VersioningConfiguration={"Status": "Enabled"})
    for key, body in objects.items():
        s3.put_object(Bucket=name, Key=key, Body=body)
    return S3Options(bucket=name, region_name=REGION, access_key="testing", secret="testing")


def _keys(s3: Any, bucket: str) -> set[str]:
    return {o["Key"] for o in s3.list_objects_v2(Bucket=bucket).get("Contents", [])}


def _versions(s3: Any, bucket: str) -> int:
    listing = s3.list_object_versions(Bucket=bucket)
    return len(listing.get("Versions", [])) + len(listing.get("DeleteMarkers", []))


class TestDeleteFromArchive:
    def test_deletes_all_but_the_logs_and_leaves_a_tombstone(self, s3: Any):
        other = {f"{OTHER_SUBMISSION_ID}/files/a.c4gh": b"z"}
        archive = _bucket(s3, "consented", {**ARCHIVED, **other})

        deletion = delete_from_archive(archive, "consented archive", SUBMISSION_ID)

        assert _keys(s3, "consented") == {
            f"{SUBMISSION_ID}/logs/progress_encrypt.cjson",
            f"{SUBMISSION_ID}/logs/progress_upload.cjson",
            f"{SUBMISSION_ID}/deleted",
            f"{OTHER_SUBMISSION_ID}/files/a.c4gh",
        }
        assert (deletion.objects_deleted, deletion.bytes_deleted, deletion.logs_kept) == (3, 165, 2)
        assert archive_holdings(archive, SUBMISSION_ID).objects == 0

    def test_none_for_an_archive_that_never_held_the_submission(self, s3: Any):
        archive = _bucket(s3, "non-consented", {f"{OTHER_SUBMISSION_ID}/files/a.c4gh": b"z"})

        assert delete_from_archive(archive, "non-consented archive", SUBMISSION_ID) is None
        assert _keys(s3, "non-consented") == {f"{OTHER_SUBMISSION_ID}/files/a.c4gh"}

    def test_deletes_old_versions_and_delete_markers(self, s3: Any):
        archive = _bucket(s3, "consented", ARCHIVED, versioned=True)
        s3.put_object(Bucket="consented", Key=f"{SUBMISSION_ID}/metadata/metadata.json", Body=b"{}")
        s3.delete_object(Bucket="consented", Key=f"{SUBMISSION_ID}/files/sample_R2.fastq.gz.c4gh")

        deletion = delete_from_archive(archive, "consented archive", SUBMISSION_ID)

        assert (deletion.objects_deleted, deletion.delete_markers_removed) == (4, 1)
        assert archive_holdings(archive, SUBMISSION_ID).objects == 0
        assert _versions(s3, "consented") == 3  # two logs and the tombstone

    def test_continues_an_interrupted_deletion(self, s3: Any):
        archive = _bucket(s3, "consented", {**ARCHIVED, f"{SUBMISSION_ID}/deleted": b""})

        deletion = delete_from_archive(archive, "consented archive", SUBMISSION_ID)

        assert deletion.objects_deleted == 3
        assert archive_holdings(archive, SUBMISSION_ID).objects == 0

    def test_running_it_again_deletes_nothing_more(self, s3: Any):
        archive = _bucket(s3, "consented", ARCHIVED)
        delete_from_archive(archive, "consented archive", SUBMISSION_ID)

        deletion = delete_from_archive(archive, "consented archive", SUBMISSION_ID)

        assert (deletion.objects_deleted, deletion.logs_kept) == (0, 2)

    def test_holdings_count_what_deleting_removes(self, s3: Any):
        archive = _bucket(s3, "consented", ARCHIVED)

        assert archive_holdings(archive, SUBMISSION_ID) == Holdings(objects=3, bytes=165)

    def test_refused_deletion_raises(self, s3: Any, monkeypatch: pytest.MonkeyPatch):
        archive = _bucket(s3, "consented", ARCHIVED)
        real_init = storage.init_s3_client

        def refusing_client(s3_options: S3Options) -> Any:
            client = real_init(s3_options)
            client.delete_objects = lambda **_: {"Errors": [{"Key": "k", "Code": "AccessDenied"}]}
            return client

        monkeypatch.setattr(storage, "init_s3_client", refusing_client)

        with pytest.raises(StorageDeletionError, match="AccessDenied"):
            delete_from_archive(archive, "consented archive", SUBMISSION_ID)


class TestObjectLock:
    def test_detects_object_lock(self, s3: Any):
        s3.create_bucket(Bucket="locked", ObjectLockEnabledForBucket=True)
        locked = S3Options(bucket="locked", region_name=REGION, access_key="testing", secret="testing")

        assert has_object_lock(locked) is True

    def test_a_plain_bucket_has_none(self, s3: Any):
        assert has_object_lock(_bucket(s3, "plain", {})) is False


class TestDeleteFromInbox:
    def test_cleans_the_inbox(self, s3: Any):
        inbox = _bucket(s3, "inbox", UPLOADED)

        deletion = delete_from_inbox(inbox, "inbox 'inbox'", SUBMISSION_ID)

        assert _keys(s3, "inbox") == {f"{SUBMISSION_ID}/metadata/metadata.json", f"{SUBMISSION_ID}/cleaned"}
        metadata = s3.get_object(Bucket="inbox", Key=f"{SUBMISSION_ID}/metadata/metadata.json")
        assert metadata["Body"].read() == b""
        assert deletion.objects_deleted == 3
        assert inbox_holdings(inbox, SUBMISSION_ID).objects == 0

    def test_a_cleaned_inbox_is_left_as_it_is(self, s3: Any):
        inbox = _bucket(s3, "inbox", {f"{SUBMISSION_ID}/metadata/metadata.json": b"", f"{SUBMISSION_ID}/cleaned": b""})

        deletion = delete_from_inbox(inbox, "inbox 'inbox'", SUBMISSION_ID)

        assert deletion.objects_deleted == 0
        assert inbox_holdings(inbox, SUBMISSION_ID).objects == 0

    def test_none_for_an_inbox_without_the_submission(self, s3: Any):
        inbox = _bucket(s3, "inbox", {})

        assert delete_from_inbox(inbox, "inbox 'inbox'", SUBMISSION_ID) is None
        assert _keys(s3, "inbox") == set()

    def test_deletes_the_versions_cleaning_leaves_behind(self, s3: Any):
        inbox = _bucket(s3, "inbox", UPLOADED, versioned=True)

        delete_from_inbox(inbox, "inbox 'inbox'", SUBMISSION_ID)

        assert inbox_holdings(inbox, SUBMISSION_ID).objects == 0
        assert _versions(s3, "inbox") == 2  # the empty metadata.json and the cleaned marker
        metadata = s3.get_object(Bucket="inbox", Key=f"{SUBMISSION_ID}/metadata/metadata.json")
        assert metadata["Body"].read() == b""

    def test_a_cleaned_versioned_inbox_loses_the_old_versions(self, s3: Any):
        """An inbox cleaned before versioning was looked at still keeps the uploaded data in old versions."""
        inbox = _bucket(s3, "inbox", UPLOADED, versioned=True)
        s3.put_object(Bucket="inbox", Key=f"{SUBMISSION_ID}/metadata/metadata.json", Body=b"")
        s3.delete_object(Bucket="inbox", Key=f"{SUBMISSION_ID}/files/sample_R1.fastq.gz.c4gh")
        s3.delete_object(Bucket="inbox", Key=f"{SUBMISSION_ID}/version")
        s3.put_object(Bucket="inbox", Key=f"{SUBMISSION_ID}/cleaned", Body=b"")
        assert inbox_holdings(inbox, SUBMISSION_ID).objects > 0

        delete_from_inbox(inbox, "inbox 'inbox'", SUBMISSION_ID)

        assert inbox_holdings(inbox, SUBMISSION_ID).objects == 0
