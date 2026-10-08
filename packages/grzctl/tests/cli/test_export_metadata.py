"""Tests for `grzctl db export-metadata`.

The unit tests feed the export an archive lookup directly. The command tests at the end run it
against a moto-mocked S3 with both archives, as the backfill tests do.
"""

import copy
import datetime
import hashlib
import importlib.resources
import io
import json
import zipfile
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import boto3
import botocore.exceptions
import crypt4gh.keys
import crypt4gh.keys.c4gh
import crypt4gh.lib
import grzctl.cli
import grzctl.commands.db.cli
import pytest
from click.testing import CliRunner
from grz_db.models.submission import Submission
from grz_pydantic_models.submission.metadata import (
    REDACTED_LOCAL_CASE_ID,
    REDACTED_TAN,
    GrzSubmissionMetadata,
    SubmissionType,
    redact_metadata_dict,
)
from grz_pydantic_models_testing.example_metadata import grzctl as grzctl_metadata
from grzctl.commands.db.export import (
    ArchivedMetadataError,
    ExportEntry,
    SkippedSubmission,
    build_export_entry,
    restore_metadata_dict,
    write_metadata_zip,
)
from moto import mock_aws

TEST_METADATA_PATH = importlib.resources.files(grzctl_metadata).joinpath("metadata.json")


def _original() -> dict:
    return json.loads(TEST_METADATA_PATH.read_text())


def _archived_json(submitted: dict) -> str:
    """The copy archiving uploads: the submitted document, redacted, written as JSON again."""
    return json.dumps(redact_metadata_dict(submitted), indent=2)


def _submission(suffix: str, **fields) -> Submission:
    """A Submission row, not stored anywhere, holding the example's tanG and localCaseId."""
    original = _original()
    defaults = {
        "id": f"260914050_2025-09-15_c64603{suffix}",
        "tan_g": original["submission"]["tanG"],
        "local_case_id": original["submission"]["localCaseId"],
        "submission_type": SubmissionType.initial,
    }
    return Submission(**(defaults | fields))


ARCHIVED_JSON = _archived_json(_original())
"""The example metadata as an archive holds it."""


def _entry(submission: Submission, raw_json: str = ARCHIVED_JSON, archive: str = "consented") -> ExportEntry:
    return build_export_entry(submission, archive, raw_json)


def test_restore_round_trip():
    """Restoring a redacted document with the database values gives back what the submitter sent."""
    original = _original()
    redacted = redact_metadata_dict(original)

    restored, unrestored = restore_metadata_dict(
        redacted,
        tan_g=original["submission"]["tanG"],
        local_case_id=original["submission"]["localCaseId"],
    )

    assert restored == original
    assert unrestored == frozenset()


def test_restore_keeps_real_values():
    """A field that already holds a real value is not overwritten by the database value."""
    original = _original()

    restored, unrestored = restore_metadata_dict(original, tan_g="b" * 64, local_case_id="other")

    assert restored == original
    assert unrestored == frozenset()


def test_restore_empty_local_case_id_is_placeholder():
    """An empty localCaseId, which the archive also holds as a placeholder, is restored."""
    original = _original()
    redacted = redact_metadata_dict(original)
    redacted["submission"]["localCaseId"] = ""

    restored, unrestored = restore_metadata_dict(
        redacted,
        tan_g=original["submission"]["tanG"],
        local_case_id=original["submission"]["localCaseId"],
    )

    assert restored["submission"]["localCaseId"] == "valid_submission"
    assert unrestored == frozenset()


def test_restore_reports_missing_values():
    """Without database values the placeholders stay, and the fields are reported as unrestored."""
    redacted = redact_metadata_dict(_original())

    restored, unrestored = restore_metadata_dict(redacted, tan_g=None, local_case_id=None)

    assert restored["submission"]["tanG"] == REDACTED_TAN
    assert restored["submission"]["localCaseId"] == REDACTED_LOCAL_CASE_ID
    assert unrestored == frozenset({"tan_g", "local_case_id"})


def test_restore_leaves_input_untouched():
    """The redacted document passed in is not modified."""
    original = _original()
    redacted = redact_metadata_dict(original)
    before = copy.deepcopy(redacted)

    restore_metadata_dict(
        redacted,
        tan_g=original["submission"]["tanG"],
        local_case_id=original["submission"]["localCaseId"],
    )

    assert redacted == before


def test_build_entry_exports_what_was_submitted():
    """The archived document, restored, is what the submitter sent, numbers written as they were."""
    entry = _entry(_submission("a1"))

    assert entry.submission_id == "260914050_2025-09-15_c64603a1"
    assert entry.archive == "consented"
    assert entry.content == _original()
    # == treats 30 and 30.0 as equal; the JSON text does not
    assert json.dumps(entry.content) == json.dumps(_original())
    assert entry.unrestored == frozenset()
    assert entry.metadata_version == "1.3.0"
    assert entry.latest_state is None


@pytest.mark.parametrize(
    "raw_json",
    ["{not json", "[]", '{"donors": []}'],
    ids=["not-json", "not-an-object", "no-submission"],
)
def test_build_entry_rejects_an_unreadable_archived_copy(raw_json: str):
    """An archived copy that is not a metadata document cannot be exported, naming its archive."""
    with pytest.raises(ArchivedMetadataError, match=r"metadata\.json in the non_consented archive cannot be read"):
        _entry(_submission("a3"), raw_json, "non_consented")


def test_build_entry_reports_unrestored_fields():
    """A submission whose tanG the database does not hold is exported with the field reported."""
    entry = _entry(_submission("a6", tan_g=None))

    assert entry.content["submission"]["tanG"] == REDACTED_TAN
    assert entry.unrestored == frozenset({"tan_g"})


SCHEMA_URL = "https://raw.githubusercontent.com/BfArM-MVH/MVGenomseq/refs/tags/{}/GRZ/grz-schema.json"
NO_SCHEMA = object()
"""Stands for a document without a $schema key."""


@pytest.mark.parametrize(
    ("schema", "version"),
    [
        (SCHEMA_URL.format("v1.2.1"), "1.2.1"),
        (SCHEMA_URL.format("v1.1.9"), "1.1.9"),
        (SCHEMA_URL.format("v1.3"), "1.3.0"),
        (NO_SCHEMA, None),
        ("", None),
        ("https://example.org/schema.json", None),
        (None, None),
        (13, None),
    ],
    ids=["three-part", "older", "two-part", "missing", "empty", "unknown-url", "null", "not-a-string"],
)
def test_build_entry_reads_the_schema_version(schema: object, version: str | None):
    """The version comes from the document's own $schema; without a known URL it is exported without one."""
    document = redact_metadata_dict(_original())
    if schema is NO_SCHEMA:
        del document["$schema"]
    else:
        document["$schema"] = schema

    entry = _entry(_submission("a7"), json.dumps(document))

    assert entry.metadata_version == version


CREATED_AT = datetime.datetime(2026, 10, 6, 9, 30, tzinfo=datetime.UTC)


def _export(tmp_path: Path, entries: list[ExportEntry], skipped=()) -> tuple[Path, dict]:
    """Write an export into *tmp_path*; returns the zip's path and the manifest."""
    output = tmp_path / "export.zip"
    manifest = write_metadata_zip(entries, list(skipped), output, include_test_submissions=False, created_at=CREATED_AT)
    return output, manifest


def test_zip_holds_each_metadata_and_the_manifest(tmp_path: Path):
    """The zip holds one metadata.json per exported submission, restored, plus the manifest."""
    entries = [_entry(_submission("b1")), _entry(_submission("b2"))]

    output, _ = _export(tmp_path, entries)

    with zipfile.ZipFile(output) as archive:
        assert sorted(archive.namelist()) == [
            "260914050_2025-09-15_c64603b1/metadata.json",
            "260914050_2025-09-15_c64603b2/metadata.json",
            "manifest.json",
        ]
        assert json.loads(archive.read("260914050_2025-09-15_c64603b1/metadata.json")) == _original()


def test_manifest_describes_the_export(tmp_path: Path):
    """The manifest in the zip lists every exported and skipped submission, with checksums that match."""
    entries = [_entry(_submission("c1", tan_g=None), archive="non_consented")]
    skipped = [SkippedSubmission("260914050_2025-09-15_c64603c2", "metadata.json found in neither archive")]

    output, manifest = _export(tmp_path, entries, skipped)

    with zipfile.ZipFile(output) as archive:
        assert json.loads(archive.read("manifest.json")) == manifest
        exported_bytes = archive.read("260914050_2025-09-15_c64603c1/metadata.json")

    assert manifest["export"]["created_at"] == "2026-10-06T09:30:00+00:00"
    assert manifest["export"]["exported_count"] == 1
    assert manifest["export"]["skipped_count"] == 1
    assert manifest["export"]["test_submissions_included"] is False
    assert manifest["export"]["grzctl_versions"]["grzctl"]
    [described] = manifest["submissions"]
    assert described == {
        "submission_id": "260914050_2025-09-15_c64603c1",
        "path": "260914050_2025-09-15_c64603c1/metadata.json",
        "archive": "non_consented",
        "sha256": hashlib.sha256(exported_bytes).hexdigest(),
        "metadata_version": "1.3.0",
        "submission_uploaded_date": None,
        "latest_state": None,
        "unrestored_fields": ["tan_g"],
    }
    assert manifest["skipped"] == [
        {"submission_id": "260914050_2025-09-15_c64603c2", "reason": "metadata.json found in neither archive"}
    ]


def test_manifest_holds_no_tan_g_or_local_case_id(tmp_path: Path):
    """The manifest can be inspected without exposing what the export protects."""
    entries = [_entry(_submission("d1"))]

    output, _ = _export(tmp_path, entries)

    with zipfile.ZipFile(output) as archive:
        manifest_text = archive.read("manifest.json").decode("utf-8")
    assert _original()["submission"]["tanG"] not in manifest_text
    assert _original()["submission"]["localCaseId"] not in manifest_text


def test_refuses_to_overwrite_an_existing_export(tmp_path: Path):
    """An existing file at the output path is left as it is."""
    entries = [_entry(_submission("e1"))]
    output = tmp_path / "export.zip"
    output.write_text("an earlier export")

    with pytest.raises(FileExistsError):
        write_metadata_zip(entries, [], output, include_test_submissions=False)

    assert output.read_text() == "an earlier export"


def test_failed_export_leaves_nothing_behind(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """A failure halfway through writing leaves neither the zip nor its temporary file."""
    entries = [_entry(_submission("f1"))]

    def disk_full(*args, **kwargs):
        raise OSError("No space left on device")

    monkeypatch.setattr(zipfile.ZipFile, "writestr", disk_full)

    with pytest.raises(OSError, match="No space left"):
        _export(tmp_path, entries)

    assert list(tmp_path.iterdir()) == []


UPLOAD_DATE = "2025-09-16"


@pytest.fixture
def archives() -> Iterator[Any]:
    """A moto-mocked S3 with the consented and non_consented archive buckets the test config names."""
    with mock_aws():
        s3_client = boto3.client("s3", region_name="us-east-1")
        for bucket in ("consented", "non_consented"):
            s3_client.create_bucket(Bucket=bucket)
        yield s3_client


def _put_archived(s3_client: Any, bucket: str, submission_id: str, submitted: dict) -> None:
    """Put a submission's metadata.json into an archive bucket, redacted, as archiving does."""
    s3_client.put_object(
        Bucket=bucket, Key=f"{submission_id}/metadata/metadata.json", Body=_archived_json(submitted).encode()
    )


@pytest.fixture
def recipient_key_paths(tmp_path: Path) -> tuple[Path, Path]:
    """A fresh crypt4gh key pair of the export's recipient, as (public key path, private key path)."""
    public_key_path, private_key_path = tmp_path / "recipient.pub", tmp_path / "recipient.sec"
    crypt4gh.keys.c4gh.generate(private_key_path, public_key_path, None, comment=None)
    return public_key_path, private_key_path


def _decrypt(encrypted_path: Path, private_key_path: Path) -> zipfile.ZipFile:
    """Decrypt an export in memory, as its recipient does, and open the zip inside."""
    private_key = crypt4gh.keys.get_private_key(private_key_path, lambda: None)
    decrypted = io.BytesIO()
    with open(encrypted_path, "rb") as encrypted:
        crypt4gh.lib.decrypt([(0, private_key, None)], encrypted, decrypted)
    decrypted.seek(0)
    return zipfile.ZipFile(decrypted)


def _invoke(config_path: Path, *args: str):
    """Run grzctl with *config_path* and fail the test right away if the command fails."""
    result = CliRunner().invoke(grzctl.cli.build_cli(), ["--config", str(config_path), *args])
    assert result.exit_code == 0, result.output
    return result


def _populate(config_path: Path, tmp_path: Path, **submission_fields) -> tuple[str, dict]:
    """Add and populate one submission from the example metadata, as the processing pipeline does.

    :param submission_fields: Fields of the metadata's ``submission`` object to override.
    :returns: The submission ID and the metadata exactly as it was submitted.
    """
    submitted = _original()
    submitted["submission"].update(submission_fields)
    submission_id = GrzSubmissionMetadata.model_validate(submitted).submission_id

    metadata_path = tmp_path / f"{submission_id}.metadata.json"
    metadata_path.write_text(json.dumps(submitted))
    _invoke(config_path, "db", "submission", "add", submission_id)
    _invoke(
        config_path,
        *("db", "submission", "populate", "--no-confirm", "--submission-date", UPLOAD_DATE),
        *(submission_id, str(metadata_path)),
    )
    return submission_id, submitted


def _unwrapped(text: str) -> str:
    """*text* on one line: errors are wrapped at the terminal width, which a long temporary path can cross."""
    return " ".join(text.split())


def _export_args(output: Path, public_key_path: Path, *extra: str) -> list[str]:
    return ["db", "export-metadata", "--output", str(output), "--public-key", str(public_key_path), *extra]


def test_command_exports_what_was_submitted(
    migrated_database_config_path: Path, tmp_path: Path, recipient_key_paths: tuple[Path, Path], archives: Any
):
    """The decrypted metadata.json is what the submitter sent, tanG and localCaseId included, from either archive."""
    public_key_path, private_key_path = recipient_key_paths
    first_id, first_submitted = _populate(
        migrated_database_config_path, tmp_path, submissionType="initial", tanG="b" * 64
    )
    second_id, second_submitted = _populate(
        migrated_database_config_path, tmp_path, submissionType="initial", tanG="c" * 64, localCaseId="case-2"
    )
    _put_archived(archives, "consented", first_id, first_submitted)
    _put_archived(archives, "non_consented", second_id, second_submitted)
    _invoke(migrated_database_config_path, "db", "submission", "update", first_id, "Finished")
    output = tmp_path / "export.zip.c4gh"

    result = _invoke(migrated_database_config_path, *_export_args(output, public_key_path))

    assert "Exported 2 metadata.json file(s)" in result.stdout
    assert hashlib.sha256(output.read_bytes()).hexdigest() in result.stdout
    with _decrypt(output, private_key_path) as archive:
        first_exported = json.loads(archive.read(f"{first_id}/metadata.json"))
        second_exported = json.loads(archive.read(f"{second_id}/metadata.json"))
        manifest = json.loads(archive.read("manifest.json"))
    # == treats 30 and 30.0 as equal; the JSON text does not
    assert json.dumps(first_exported) == json.dumps(first_submitted)
    assert json.dumps(second_exported) == json.dumps(second_submitted)
    described = {submission["submission_id"]: submission for submission in manifest["submissions"]}
    assert described[first_id]["archive"] == "consented"
    assert described[second_id]["archive"] == "non_consented"
    assert described[first_id]["latest_state"] == "Finished"
    assert described[second_id]["latest_state"] is None
    assert described[first_id]["submission_uploaded_date"] == UPLOAD_DATE
    assert described[first_id]["unrestored_fields"] == []
    assert manifest["skipped"] == []


def test_command_writes_nothing_readable_without_the_private_key(
    migrated_database_config_path: Path, tmp_path: Path, recipient_key_paths: tuple[Path, Path], archives: Any
):
    """The export is no readable zip, and holds no tanG in plain text."""
    public_key_path, _ = recipient_key_paths
    submission_id, submitted = _populate(
        migrated_database_config_path, tmp_path, submissionType="initial", tanG="d" * 64
    )
    _put_archived(archives, "consented", submission_id, submitted)
    output = tmp_path / "export.zip.c4gh"

    _invoke(migrated_database_config_path, *_export_args(output, public_key_path))

    assert not zipfile.is_zipfile(output)
    assert submitted["submission"]["tanG"].encode() not in output.read_bytes()
    assert list(tmp_path.glob(".export.zip.c4gh.*")) == []


def test_command_skips_what_it_cannot_export(
    migrated_database_config_path: Path, tmp_path: Path, recipient_key_paths: tuple[Path, Path], archives: Any
):
    """Test submissions, and submissions in neither or both archives, are skipped on stderr and in the manifest."""
    public_key_path, private_key_path = recipient_key_paths
    test_id, test_submitted = _populate(migrated_database_config_path, tmp_path)  # the example is a test submission
    _put_archived(archives, "consented", test_id, test_submitted)
    unarchived_id = "260914050_2025-09-15_0000abcd"
    _invoke(migrated_database_config_path, "db", "submission", "add", unarchived_id)
    doubled_id, doubled_submitted = _populate(
        migrated_database_config_path, tmp_path, submissionType="initial", tanG="e" * 64
    )
    _put_archived(archives, "consented", doubled_id, doubled_submitted)
    _put_archived(archives, "non_consented", doubled_id, doubled_submitted)
    output = tmp_path / "export.zip.c4gh"

    result = _invoke(migrated_database_config_path, *_export_args(output, public_key_path))

    reasons = {
        test_id: "test submission",
        unarchived_id: "metadata.json found in neither archive",
        doubled_id: "metadata.json found in both consented and non_consented archives",
    }
    for submission_id, reason in reasons.items():
        assert f"Skipped {submission_id}: {reason}." in result.stderr
    with _decrypt(output, private_key_path) as archive:
        assert archive.namelist() == ["manifest.json"]
        manifest = json.loads(archive.read("manifest.json"))
    assert {skip["submission_id"]: skip["reason"] for skip in manifest["skipped"]} == reasons


def test_command_skips_when_an_archive_cannot_be_read(
    migrated_database_config_path: Path,
    tmp_path: Path,
    recipient_key_paths: tuple[Path, Path],
    archives: Any,
    monkeypatch: pytest.MonkeyPatch,
):
    """An archive that cannot be read might hold a second copy, so the copy in the other one is not exported."""
    public_key_path, _ = recipient_key_paths
    submission_id, submitted = _populate(
        migrated_database_config_path, tmp_path, submissionType="initial", tanG="f" * 64
    )
    _put_archived(archives, "non_consented", submission_id, submitted)
    fetch = grzctl.commands.db.cli._fetch_metadata_json

    def consented_archive_fails(s3_client: Any, bucket: str, submission_id: str) -> str | None:
        if bucket == "consented":
            raise botocore.exceptions.ClientError({"Error": {"Code": "InternalError"}}, "GetObject")
        return fetch(s3_client, bucket, submission_id)

    monkeypatch.setattr(grzctl.commands.db.cli, "_fetch_metadata_json", consented_archive_fails)
    output = tmp_path / "export.zip.c4gh"

    result = _invoke(migrated_database_config_path, *_export_args(output, public_key_path))

    assert f"Skipped {submission_id}: Reading metadata.json from the consented archive failed" in result.stderr
    assert "Exported 0 metadata.json file(s)" in result.stdout


def test_command_stops_when_an_archive_bucket_is_missing(
    migrated_database_config_path: Path, tmp_path: Path, recipient_key_paths: tuple[Path, Path], archives: Any
):
    """A missing bucket is a faulty setup that fails every submission alike, so the export stops."""
    public_key_path, _ = recipient_key_paths
    submission_id, submitted = _populate(
        migrated_database_config_path, tmp_path, submissionType="initial", tanG="f" * 64
    )
    _put_archived(archives, "non_consented", submission_id, submitted)
    archives.delete_bucket(Bucket="consented")
    output = tmp_path / "export.zip.c4gh"

    result = CliRunner().invoke(
        grzctl.cli.build_cli(),
        ["--config", str(migrated_database_config_path), *_export_args(output, public_key_path)],
    )

    assert result.exit_code != 0
    assert "Reading metadata.json from the consented archive failed" in _unwrapped(result.stderr)
    assert not output.exists()


def test_command_skips_test_submissions_without_reading_the_archives(
    migrated_database_config_path: Path, tmp_path: Path, recipient_key_paths: tuple[Path, Path], archives: Any
):
    """A test submission is skipped before its archives are read: with both buckets gone, the export still runs."""
    public_key_path, _ = recipient_key_paths
    test_id, _ = _populate(migrated_database_config_path, tmp_path)
    for bucket in ("consented", "non_consented"):
        archives.delete_bucket(Bucket=bucket)

    result = _invoke(migrated_database_config_path, *_export_args(tmp_path / "export.zip.c4gh", public_key_path))

    assert f"Skipped {test_id}: test submission." in result.stderr


def test_command_includes_test_submissions_on_request(
    migrated_database_config_path: Path, tmp_path: Path, recipient_key_paths: tuple[Path, Path], archives: Any
):
    """--include-test-submissions exports test submissions too."""
    public_key_path, private_key_path = recipient_key_paths
    test_id, submitted = _populate(migrated_database_config_path, tmp_path)
    _put_archived(archives, "consented", test_id, submitted)
    output = tmp_path / "export.zip.c4gh"

    _invoke(migrated_database_config_path, *_export_args(output, public_key_path, "--include-test-submissions"))

    with _decrypt(output, private_key_path) as archive:
        assert json.loads(archive.read(f"{test_id}/metadata.json")) == submitted
        assert json.loads(archive.read("manifest.json"))["export"]["test_submissions_included"] is True


def test_command_refuses_to_overwrite_before_reading_the_archives(
    migrated_database_config_path: Path, tmp_path: Path, recipient_key_paths: tuple[Path, Path]
):
    """An existing output file aborts the command before any archive is read, and is left as it is."""
    public_key_path, _ = recipient_key_paths
    output = tmp_path / "export.zip.c4gh"
    output.write_text("an earlier export")

    # no mocked S3: reading an archive would fail
    result = CliRunner().invoke(
        grzctl.cli.build_cli(),
        ["--config", str(migrated_database_config_path), *_export_args(output, public_key_path)],
    )

    assert result.exit_code != 0
    assert "Refusing to overwrite" in _unwrapped(result.stderr)
    assert output.read_text() == "an earlier export"


def test_command_rejects_an_unreadable_public_key(migrated_database_config_path: Path, tmp_path: Path):
    """A file that is no crypt4gh public key aborts the command before anything is written."""
    not_a_key = tmp_path / "not-a-key.pub"
    not_a_key.write_text("this is no key")
    output = tmp_path / "export.zip.c4gh"

    result = CliRunner().invoke(
        grzctl.cli.build_cli(),
        ["--config", str(migrated_database_config_path), *_export_args(output, not_a_key)],
    )

    assert result.exit_code != 0
    assert "cannot be read" in _unwrapped(result.stderr)
    assert not output.exists()
