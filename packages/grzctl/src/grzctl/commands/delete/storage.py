"""Delete what the archives and the inbox hold of a submission.

A bucket that keeps old object versions would keep the data in them, so those are deleted too.
Not every S3-compatible store answers the versioning and object lock queries; one that does not
is taken to have neither. The check after deleting, :func:`archive_holdings` and
:func:`inbox_holdings`, catches anything that is left all the same.
"""

from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from typing import Any

from botocore.exceptions import ClientError
from grz_common.models.s3 import S3Options
from grz_common.transfer import init_s3_client, s3_errors

from ..clean import clean_submission_from_inbox

TOMBSTONE = "deleted"
"""Name of the empty object ``grzctl delete run`` leaves in the archive that held a submission."""

_DELETE_BATCH_SIZE = 1000
"""Most objects one DeleteObjects request may name."""


@dataclass(frozen=True)
class BucketDeletion:
    """What was deleted of a submission in one bucket.

    :param location: The bucket's role, e.g. ``consented archive``.
    :param objects_deleted: How many object versions were deleted, the current ones included.
    :param bytes_deleted: Their total size.
    :param delete_markers_removed: How many delete markers of a versioned bucket were removed.
    :param logs_kept: How many log files were kept.
    """

    location: str
    objects_deleted: int
    bytes_deleted: int
    delete_markers_removed: int
    logs_kept: int


class StorageDeletionError(Exception):
    """S3 refused to delete some objects of a submission."""


@dataclass(frozen=True)
class _Entry:
    key: str
    version_id: str | None
    size: int
    is_delete_marker: bool
    is_latest: bool


def _is_versioned(client: Any, bucket: str) -> bool:
    """Whether a bucket keeps, or has kept, old object versions."""
    try:
        status = client.get_bucket_versioning(Bucket=bucket).get("Status")
    except ClientError:
        return False
    return status in ("Enabled", "Suspended")


def has_object_lock(s3_options: S3Options) -> bool:
    """Whether a bucket has object lock enabled, which keeps objects from being deleted.

    :param s3_options: S3 connection to the bucket.
    """
    client = init_s3_client(s3_options)
    try:
        configuration = client.get_object_lock_configuration(Bucket=s3_options.bucket)
    except ClientError:
        return False
    return configuration.get("ObjectLockConfiguration", {}).get("ObjectLockEnabled") == "Enabled"


def _entries(client: Any, bucket: str, prefix: str, versioned: bool) -> list[_Entry]:
    """Every object under *prefix*, with each old version and delete marker if *versioned*."""
    entries = []
    if versioned:
        for page in client.get_paginator("list_object_versions").paginate(Bucket=bucket, Prefix=prefix):
            entries += [
                _Entry(v["Key"], v["VersionId"], v.get("Size", 0), False, v["IsLatest"])
                for v in page.get("Versions", [])
            ]
            entries += [_Entry(m["Key"], m["VersionId"], 0, True, m["IsLatest"]) for m in page.get("DeleteMarkers", [])]
    else:
        for page in client.get_paginator("list_objects_v2").paginate(Bucket=bucket, Prefix=prefix):
            entries += [_Entry(o["Key"], None, o.get("Size", 0), False, True) for o in page.get("Contents", [])]
    return entries


def _batches(entries: Sequence[_Entry]) -> Iterator[Sequence[_Entry]]:
    for start in range(0, len(entries), _DELETE_BATCH_SIZE):
        yield entries[start : start + _DELETE_BATCH_SIZE]


def _delete(client: Any, bucket: str, entries: Sequence[_Entry]) -> None:
    """Delete each entry, each version permanently.

    :raises StorageDeletionError: if S3 refuses to delete any of them.
    """
    for batch in _batches(entries):
        objects: list[dict[str, str]] = [
            {"Key": e.key, "VersionId": e.version_id} if e.version_id is not None else {"Key": e.key} for e in batch
        ]
        response = client.delete_objects(Bucket=bucket, Delete={"Objects": objects, "Quiet": True})
        if errors := response.get("Errors"):
            codes = sorted({error.get("Code", "unknown") for error in errors})
            raise StorageDeletionError(f"S3 refused to delete {len(errors)} objects in {bucket}: {', '.join(codes)}")


def _deletion(location: str, deleted: Sequence[_Entry], logs_kept: int) -> BucketDeletion:
    return BucketDeletion(
        location=location,
        objects_deleted=sum(not e.is_delete_marker for e in deleted),
        bytes_deleted=sum(e.size for e in deleted),
        delete_markers_removed=sum(e.is_delete_marker for e in deleted),
        logs_kept=logs_kept,
    )


def delete_from_archive(s3_options: S3Options, location: str, submission_id: str) -> BucketDeletion | None:
    """Delete what an archive holds of a submission, except its log files.

    Before deleting, leaves an empty ``<submission ID>/deleted`` object, so that the archive shows
    which submission it held even when the deletion is interrupted. Running it again continues an
    interrupted deletion.

    :param s3_options: S3 connection to the archive.
    :param location: The archive's role, for messages and the returned record.
    :param submission_id: ID of the submission.
    :returns: What was deleted, or ``None`` if the archive never held the submission.
    :raises StorageDeletionError: if S3 refuses to delete some objects.
    """
    client = init_s3_client(s3_options)
    bucket = s3_options.bucket
    prefix = f"{submission_id}/"
    tombstone = f"{prefix}{TOMBSTONE}"
    with s3_errors(f"Deleting {submission_id} from the {location}"):
        entries = _entries(client, bucket, prefix, _is_versioned(client, bucket))
        logs = [e for e in entries if e.key.startswith(f"{prefix}logs/")]
        to_delete = [e for e in entries if not e.key.startswith(f"{prefix}logs/") and e.key != tombstone]
        has_tombstone = any(e.key == tombstone and e.is_latest and not e.is_delete_marker for e in entries)
        if not to_delete and not has_tombstone:
            return None
        if not has_tombstone:
            client.put_object(Bucket=bucket, Key=tombstone, Body=b"")
        _delete(client, bucket, to_delete)
    return _deletion(location, to_delete, logs_kept=sum(e.is_latest and not e.is_delete_marker for e in logs))


def delete_from_inbox(s3_options: S3Options, location: str, submission_id: str) -> BucketDeletion | None:
    """Clean a submission from an inbox, as ``grzctl clean`` does, and delete what versions keep of it.

    Cleaning leaves an empty ``metadata.json``, which keeps the submission ID from being uploaded
    again, and a ``cleaned`` marker. An inbox that is cleaned already is not cleaned again.

    :param s3_options: S3 connection to the inbox.
    :param location: The inbox's name, for messages and the returned record.
    :param submission_id: ID of the submission.
    :returns: What was deleted, or ``None`` if the inbox holds nothing of the submission.
    :raises StorageDeletionError: if S3 refuses to delete some objects.
    """
    client = init_s3_client(s3_options)
    bucket = s3_options.bucket
    prefix = f"{submission_id}/"
    with s3_errors(f"Deleting {submission_id} from the {location}"):
        versioned = _is_versioned(client, bucket)
        before = _entries(client, bucket, prefix, versioned)
        if not before:
            return None
        current = [e for e in before if e.is_latest and not e.is_delete_marker]
        deleted: list[_Entry] = []
        if not any(e.key == f"{prefix}cleaned" for e in current):
            # what cleaning deletes or, for metadata.json, overwrites
            deleted = [e for e in current if e.key not in (f"{prefix}cleaning", f"{prefix}cleaned")]
            clean_submission_from_inbox(bucket, s3_options, submission_id, location)
        if versioned:
            # in a versioned bucket, cleaning only hides the objects behind delete markers and
            # keeps the original metadata.json as an old version: deleting the versions removes them
            kept = {f"{prefix}metadata/metadata.json", f"{prefix}cleaned"}
            deleted = [
                e
                for e in _entries(client, bucket, prefix, versioned=True)
                if not (e.is_latest and not e.is_delete_marker and e.key in kept)
            ]
            _delete(client, bucket, deleted)
    return _deletion(location, deleted, logs_kept=0)


@dataclass(frozen=True)
class Holdings:
    """What a bucket holds of a submission that deleting it removes.

    :param objects: How many object versions and delete markers.
    :param bytes: Their total size.
    """

    objects: int
    bytes: int


def _holdings(entries: Sequence[_Entry]) -> Holdings:
    return Holdings(objects=len(entries), bytes=sum(e.size for e in entries))


def archive_holdings(s3_options: S3Options, submission_id: str) -> Holdings:
    """What an archive holds of a submission besides its log files and the tombstone.

    Before :func:`delete_from_archive`, this is what it will delete; after it, what it left.

    :param s3_options: S3 connection to the archive.
    :param submission_id: ID of the submission.
    """
    client = init_s3_client(s3_options)
    prefix = f"{submission_id}/"
    entries = _entries(client, s3_options.bucket, prefix, _is_versioned(client, s3_options.bucket))
    return _holdings(
        [
            e
            for e in entries
            if not e.key.startswith(f"{prefix}logs/") and not (e.key == f"{prefix}{TOMBSTONE}" and e.is_latest)
        ]
    )


def inbox_holdings(s3_options: S3Options, submission_id: str) -> Holdings:
    """What an inbox holds of a submission besides the empty ``metadata.json`` and the ``cleaned`` marker.

    Before :func:`delete_from_inbox`, this is what it will delete; after it, what it left.

    :param s3_options: S3 connection to the inbox.
    :param submission_id: ID of the submission.
    """
    client = init_s3_client(s3_options)
    prefix = f"{submission_id}/"
    entries = _entries(client, s3_options.bucket, prefix, _is_versioned(client, s3_options.bucket))
    return _holdings(
        [
            e
            for e in entries
            if not (
                e.is_latest
                and not e.is_delete_marker
                and (e.key == f"{prefix}cleaned" or (e.key == f"{prefix}metadata/metadata.json" and e.size == 0))
            )
        ]
    )
