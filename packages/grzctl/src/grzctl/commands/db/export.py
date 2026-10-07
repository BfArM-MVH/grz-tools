"""
Logic for exporting archived metadata.json files with their redacted fields restored.

The archive buckets hold each metadata.json as the submitter sent it, with tanG and localCaseId
redacted. The database holds those two values, so the export reads the document from the archive
and the values from the database.
"""

import copy
import datetime
import hashlib
import io
import json
import re
import tempfile
import zipfile
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import crypt4gh.lib
from cryptography.hazmat.primitives.asymmetric.x25519 import X25519PublicKey
from grz_common.utils.crypt import Crypt4GH
from grz_db.models.submission import Submission
from grz_pydantic_models.submission.metadata import (
    REDACTED_TAN,
    SCHEMA_URL_PATTERN,
    is_redacted_local_case_id,
)

from ... import get_versions

MANIFEST_NAME = "manifest.json"


class ArchivedMetadataError(Exception):
    """A submission's archived metadata.json cannot be exported, although both archives could be read.

    Not a :class:`~grz_common.exceptions.GrzError`: those record a failure reason in a submission's
    state, and the export records nothing.
    """


@dataclass(frozen=True)
class ExportEntry:
    """One metadata.json to export, with what the manifest says about it."""

    submission_id: str
    archive: str
    content: dict[str, Any]
    unrestored: frozenset[str]
    metadata_version: str | None
    submission_uploaded_date: datetime.date | None
    latest_state: str | None


@dataclass(frozen=True)
class SkippedSubmission:
    """A submission left out of the export, and why."""

    submission_id: str
    reason: str


def restore_metadata_dict(
    redacted: dict[str, Any], *, tan_g: str | None, local_case_id: str | None
) -> tuple[dict[str, Any], frozenset[str]]:
    """Put the submitter's tanG and localCaseId back into a redacted metadata JSON document, in a copy.

    The dict-level counterpart to
    :meth:`~grz_pydantic_models.submission.metadata.GrzSubmissionMetadata.restore_redacted_fields`,
    so the exported document keeps the submitter's structure instead of a re-serialised model.
    Only placeholders are replaced, so a field that already holds a real value stays as it is.

    The index donor's pseudonym, which redaction also replaces, cannot be restored:
    nothing stores the original value.

    :param redacted: Redacted metadata, as read from an archive bucket. Not modified.
    :param tan_g: The submitter's tanG from the database, or ``None`` if unknown.
    :param local_case_id: The submitter's localCaseId from the database, or ``None`` if unknown.
    :returns: The restored copy, and the names of the fields still redacted because no value was known.
    """
    restored = copy.deepcopy(redacted)
    submission = restored["submission"]
    unrestored = set()

    if submission["tanG"] == REDACTED_TAN:
        if tan_g:
            submission["tanG"] = tan_g
        else:
            unrestored.add("tan_g")

    if is_redacted_local_case_id(submission["localCaseId"]):
        if not is_redacted_local_case_id(local_case_id):
            submission["localCaseId"] = local_case_id
        else:
            unrestored.add("local_case_id")

    return restored, frozenset(unrestored)


def _schema_version(content: dict[str, Any]) -> str | None:
    """Read the metadata schema version, such as ``1.3.0``, from the document's ``$schema`` URL.

    Formats it like :meth:`~grz_pydantic_models.submission.metadata.GrzSubmissionMetadata.get_schema_version`,
    which needs a validated document; the export reads the version without validating it.

    :param content: Metadata document.
    :returns: The version, or ``None`` if the URL is missing, not a string, or not a known schema URL.
    """
    schema_url = content.get("$schema")
    if not isinstance(schema_url, str):
        return None
    match = re.fullmatch(SCHEMA_URL_PATTERN, schema_url)
    return ".".join(match.groups(default="0")) if match else None


def build_export_entry(submission: Submission, archive: str, raw_json: str) -> ExportEntry:
    """Restore a submission's archived metadata.json and describe it for the manifest.

    :param submission: The submission, with its ``states`` loaded, for :meth:`Submission.get_latest_state`.
    :param archive: Name of the archive the document was read from, such as ``consented``.
    :param raw_json: The redacted document, as archiving wrote it.
    :returns: The entry to export.
    :raises ArchivedMetadataError: If the document is not JSON, or holds no submission with tanG and localCaseId.
    """
    try:
        content, unrestored = restore_metadata_dict(
            json.loads(raw_json), tan_g=submission.tan_g, local_case_id=submission.local_case_id
        )
    except (ValueError, KeyError, TypeError) as e:
        raise ArchivedMetadataError(f"metadata.json in the {archive} archive cannot be read") from e

    latest_state = submission.get_latest_state()
    return ExportEntry(
        submission_id=submission.id,
        archive=archive,
        content=content,
        unrestored=unrestored,
        metadata_version=_schema_version(content),
        submission_uploaded_date=submission.submission_uploaded_date,
        latest_state=latest_state.state.value if latest_state else None,
    )


def build_manifest(
    entries: Sequence[ExportEntry],
    skipped: Sequence[SkippedSubmission],
    checksums: dict[str, str],
    *,
    created_at: datetime.datetime,
    include_test_submissions: bool,
) -> dict[str, Any]:
    """Describe an export: when and how it was made, each exported file, and what was left out.

    Holds no tanG, localCaseId or other content of the documents, so it can be inspected
    without exposing what the export protects.

    :param entries: The exported submissions.
    :param skipped: The submissions left out.
    :param checksums: SHA-256 hex digest of each exported file, by submission ID.
    :param created_at: When the export was made.
    :param include_test_submissions: Whether test submissions were exported.
    :returns: The manifest, ready to be written as JSON.
    """
    return {
        "export": {
            "created_at": created_at.isoformat(),
            "grzctl_versions": get_versions(),
            "exported_count": len(entries),
            "skipped_count": len(skipped),
            "test_submissions_included": include_test_submissions,
        },
        "submissions": [
            {
                "submission_id": entry.submission_id,
                "path": f"{entry.submission_id}/metadata.json",
                "archive": entry.archive,
                "sha256": checksums[entry.submission_id],
                "metadata_version": entry.metadata_version,
                "submission_uploaded_date": (
                    entry.submission_uploaded_date.isoformat() if entry.submission_uploaded_date else None
                ),
                "latest_state": entry.latest_state,
                "unrestored_fields": sorted(entry.unrestored),
            }
            for entry in entries
        ],
        "skipped": [{"submission_id": skip.submission_id, "reason": skip.reason} for skip in skipped],
    }


def write_metadata_zip(  # noqa: PLR0913
    entries: Sequence[ExportEntry],
    skipped: Sequence[SkippedSubmission],
    output: Path,
    *,
    include_test_submissions: bool,
    recipient_public_key: X25519PublicKey | None = None,
    created_at: datetime.datetime | None = None,
) -> dict[str, Any]:
    """Write each exported metadata.json and the manifest into a new zip file, encrypted if a key is given.

    The zip is built in memory, so with a key its unencrypted content never reaches the disk.
    It is written under a temporary name next to *output* and only renamed once complete,
    so a failed export never leaves a file that looks finished.

    :param entries: The submissions to export.
    :param skipped: The submissions left out, for the manifest.
    :param output: Path of the file to create. Must not exist yet.
    :param include_test_submissions: Whether test submissions were exported, for the manifest.
    :param recipient_public_key: Crypt4GH public key to encrypt the zip for, or ``None`` to write it unencrypted.
    :param created_at: When the export was made; now, in UTC, when not given.
    :returns: The manifest written into the zip.
    :raises FileExistsError: If *output* already exists. Nothing is written then.
    """
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite existing export '{output}'.")

    files = {
        entry.submission_id: json.dumps(entry.content, indent=2, ensure_ascii=False).encode("utf-8")
        for entry in entries
    }
    checksums = {submission_id: hashlib.sha256(data).hexdigest() for submission_id, data in files.items()}
    manifest = build_manifest(
        entries,
        skipped,
        checksums,
        created_at=created_at or datetime.datetime.now(datetime.UTC),
        include_test_submissions=include_test_submissions,
    )

    with tempfile.NamedTemporaryFile(dir=output.parent, prefix=f".{output.name}.", suffix=".tmp", delete=False) as tmp:
        tmp_path = Path(tmp.name)
    try:
        zip_buffer = io.BytesIO()
        with zipfile.ZipFile(zip_buffer, mode="w", compression=zipfile.ZIP_DEFLATED) as archive:
            for submission_id, data in files.items():
                archive.writestr(f"{submission_id}/metadata.json", data)
            archive.writestr(MANIFEST_NAME, json.dumps(manifest, indent=2, ensure_ascii=False))
        zip_buffer.seek(0)

        with open(tmp_path, "wb") as out_fd:
            if recipient_public_key is None:
                out_fd.write(zip_buffer.getbuffer())
            else:
                keys = Crypt4GH.prepare_c4gh_keys(recipient_public_key)
                crypt4gh.lib.encrypt(keys=keys, infile=zip_buffer, outfile=out_fd)
        tmp_path.replace(output)
    except BaseException:
        tmp_path.unlink(missing_ok=True)
        raise

    return manifest
