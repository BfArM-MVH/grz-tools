"""Delete what the archives and the inbox hold of a submission.

A bucket that keeps old object versions would keep the data in them, so those are deleted too.
Each listing lists versions, which needs the ``s3:ListBucketVersions`` permission, and deletes
each version by its ID, which needs ``s3:DeleteObjectVersion``; without them, deleting stops with
an error rather than leaving old versions behind. Not every S3-compatible store answers the
object lock query; one that does not is taken to have none, and the check after deleting,
:func:`archive_holdings` and :func:`inbox_holdings`, catches anything that is left.
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
class StoredObject:
    """One version of an object in a bucket, or a delete marker.

    :param key: The object's full path in the bucket.
    :param version_id: The version's ID, ``null`` in a bucket that never kept versions.
    :param size: Its size in bytes, 0 for a delete marker.
    :param is_delete_marker: Whether it is a delete marker, which hides older versions of the key.
    :param is_latest: Whether it is the version that S3 returns for the key.
    """

    key: str
    version_id: str | None
    size: int
    is_delete_marker: bool
    is_latest: bool


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


def _entries(client: Any, bucket: str, prefix: str) -> list[StoredObject]:
    """Every version and delete marker of every object under *prefix*.

    Lists versions whether the bucket keeps them or not, as a bucket that never kept any lists
    each object once, with the version ID ``null``. So no old version is missed when the bucket's
    versioning cannot be queried.

    Keeps only what is under *prefix* even when the store returns more: what this lists is
    deleted, and an S3-compatible store that ignored the prefix would otherwise hand over other
    submissions.

    :raises ValueError: if *prefix* is not a submission ID followed by ``/``, which could match
        other submissions or the whole bucket.
    """
    if not prefix.endswith("/") or not prefix.strip("/"):
        raise ValueError(f"Refusing to list '{prefix}': it is no submission ID followed by '/'.")
    entries = []
    for page in client.get_paginator("list_object_versions").paginate(Bucket=bucket, Prefix=prefix):
        entries += [
            StoredObject(v["Key"], v["VersionId"], v.get("Size", 0), False, v["IsLatest"])
            for v in page.get("Versions", [])
        ]
        entries += [
            StoredObject(m["Key"], m["VersionId"], 0, True, m["IsLatest"]) for m in page.get("DeleteMarkers", [])
        ]
    return [entry for entry in entries if entry.key.startswith(prefix)]


def _is_current(entry: StoredObject, key: str) -> bool:
    """Whether *entry* is the object that S3 currently returns for *key*."""
    return entry.key == key and entry.is_latest and not entry.is_delete_marker


def _archive_payload(entry: StoredObject, prefix: str) -> bool:
    """Whether deleting a submission removes *entry* from an archive: all but its log files and its current tombstone."""
    return not entry.key.startswith(f"{prefix}logs/") and not _is_current(entry, f"{prefix}{TOMBSTONE}")


def _inbox_payload(entry: StoredObject, prefix: str) -> bool:
    """Whether deleting a submission removes *entry* from an inbox: all but an empty ``metadata.json`` and the ``cleaned`` marker."""
    kept = _is_current(entry, f"{prefix}cleaned") or (
        _is_current(entry, f"{prefix}metadata/metadata.json") and entry.size == 0
    )
    return not kept


def _batches(entries: Sequence[StoredObject]) -> Iterator[Sequence[StoredObject]]:
    for start in range(0, len(entries), _DELETE_BATCH_SIZE):
        yield entries[start : start + _DELETE_BATCH_SIZE]


def _delete(client: Any, bucket: str, entries: Sequence[StoredObject]) -> None:
    """Delete each entry, each version permanently.

    :raises StorageDeletionError: if S3 refuses to delete any of them.
    """
    for batch in _batches(entries):
        objects = [{"Key": e.key, "VersionId": e.version_id} for e in batch if e.version_id is not None]
        response = client.delete_objects(Bucket=bucket, Delete={"Objects": objects, "Quiet": True})
        if errors := response.get("Errors"):
            codes = sorted({error.get("Code", "unknown") for error in errors})
            raise StorageDeletionError(f"S3 refused to delete {len(errors)} objects in {bucket}: {', '.join(codes)}")


def _deletion(location: str, deleted: Sequence[StoredObject], logs_kept: int) -> BucketDeletion:
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
        entries = _entries(client, bucket, prefix)
        to_delete = [e for e in entries if _archive_payload(e, prefix)]
        has_tombstone = any(_is_current(e, tombstone) for e in entries)
        if not to_delete and not has_tombstone:
            return None
        if not has_tombstone:
            client.put_object(Bucket=bucket, Key=tombstone, Body=b"")
        _delete(client, bucket, to_delete)
    logs_kept = sum(e.key.startswith(f"{prefix}logs/") and e.is_latest and not e.is_delete_marker for e in entries)
    return _deletion(location, to_delete, logs_kept=logs_kept)


def delete_from_inbox(s3_options: S3Options, location: str, submission_id: str) -> BucketDeletion | None:
    """Clean a submission from an inbox, as ``grzctl clean`` does, and delete what versions keep of it.

    Cleaning leaves an empty ``metadata.json``, which keeps the submission ID from being uploaded
    again, and a ``cleaned`` marker. It cleans whenever the inbox holds more than these, also an
    inbox cleaned before into which a restore from a backup brought the submission back.

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
        before = _entries(client, bucket, prefix)
        if not before:
            return None
        # what cleaning deletes or, for metadata.json, empties
        cleaned = [e for e in before if e.is_latest and not e.is_delete_marker and _inbox_payload(e, prefix)]
        if cleaned:
            clean_submission_from_inbox(bucket, s3_options, submission_id, location)
        # in a bucket that keeps versions, cleaning only hides the objects behind delete markers and
        # keeps the original metadata.json as an old version: deleting the versions removes them
        old_versions = [e for e in _entries(client, bucket, prefix) if _inbox_payload(e, prefix)]
        _delete(client, bucket, old_versions)
    # an object that cleaning deleted can show up again as an old version: count each version once
    deleted = list({(e.key, e.version_id): e for e in cleaned + old_versions if not e.is_delete_marker}.values())
    deleted += [e for e in old_versions if e.is_delete_marker]
    return _deletion(location, deleted, logs_kept=0)


@dataclass(frozen=True)
class Holdings:
    """What a bucket holds of a submission that deleting it removes.

    :param objects: How many object versions and delete markers.
    :param bytes: Their total size.
    """

    objects: int
    bytes: int


def _holdings(entries: Sequence[StoredObject]) -> Holdings:
    return Holdings(objects=len(entries), bytes=sum(e.size for e in entries))


@dataclass(frozen=True)
class Listing:
    """What a bucket holds of a submission, split into what deleting it removes and what it keeps.

    :param to_delete: Every object version and delete marker that deleting the submission removes.
    :param kept: The objects it keeps, such as the log files in an archive.
    """

    to_delete: list[StoredObject]
    kept: list[StoredObject]


def archive_listing(s3_options: S3Options, submission_id: str) -> Listing:
    """What an archive holds of a submission: everything but its log files and the tombstone goes.

    :param s3_options: S3 connection to the archive.
    :param submission_id: ID of the submission.
    """
    prefix = f"{submission_id}/"
    entries = _entries(init_s3_client(s3_options), s3_options.bucket, prefix)
    return Listing(
        to_delete=[e for e in entries if _archive_payload(e, prefix)],
        kept=[e for e in entries if not _archive_payload(e, prefix)],
    )


def inbox_listing(s3_options: S3Options, submission_id: str) -> Listing:
    """What an inbox holds of a submission: everything but an empty ``metadata.json`` and the ``cleaned`` marker goes.

    :param s3_options: S3 connection to the inbox.
    :param submission_id: ID of the submission.
    """
    prefix = f"{submission_id}/"
    entries = _entries(init_s3_client(s3_options), s3_options.bucket, prefix)
    return Listing(
        to_delete=[e for e in entries if _inbox_payload(e, prefix)],
        kept=[e for e in entries if not _inbox_payload(e, prefix)],
    )


def archive_holdings(s3_options: S3Options, submission_id: str) -> Holdings:
    """What an archive holds of a submission besides its log files and the tombstone.

    Before :func:`delete_from_archive`, this is what it will delete; after it, what it left.

    :param s3_options: S3 connection to the archive.
    :param submission_id: ID of the submission.
    """
    return _holdings(archive_listing(s3_options, submission_id).to_delete)


def inbox_holdings(s3_options: S3Options, submission_id: str) -> Holdings:
    """What an inbox holds of a submission besides the empty ``metadata.json`` and the ``cleaned`` marker.

    Before :func:`delete_from_inbox`, this is what it will delete; after it, what it left.

    :param s3_options: S3 connection to the inbox.
    :param submission_id: ID of the submission.
    """
    return _holdings(inbox_listing(s3_options, submission_id).to_delete)


def tombstoned_submissions(s3_options: S3Options) -> set[str]:
    """The submissions an archive holds a tombstone of, which ``grzctl delete run`` left when deleting them.

    :param s3_options: S3 connection to the archive.
    :returns: Their IDs.
    """
    client = init_s3_client(s3_options)
    submission_ids = set()
    for page in client.get_paginator("list_objects_v2").paginate(Bucket=s3_options.bucket):
        for obj in page.get("Contents", []):
            submission_id, _, name = obj["Key"].partition("/")
            if name == TOMBSTONE:
                submission_ids.add(submission_id)
    return submission_ids
