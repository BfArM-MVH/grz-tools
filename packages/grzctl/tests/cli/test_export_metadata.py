"""Tests for `grzctl db export-metadata`."""

import copy
import datetime
import hashlib
import importlib.resources
import io
import json
import zipfile
from pathlib import Path

import crypt4gh.keys
import crypt4gh.keys.c4gh
import crypt4gh.lib
import grzctl.cli
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
    ExportEntry,
    _schema_version,
    collect_export_entries,
    restore_metadata_dict,
    write_metadata_zip,
)

TEST_METADATA_PATH = importlib.resources.files(grzctl_metadata).joinpath("metadata.json")


def _original() -> dict:
    return json.loads(TEST_METADATA_PATH.read_text())


def _submission(suffix: str, **fields) -> Submission:
    """A Submission row, not stored anywhere, holding the redacted example metadata."""
    original = _original()
    defaults = {
        "id": f"260914050_2025-09-15_c64603{suffix}",
        "tan_g": original["submission"]["tanG"],
        "local_case_id": original["submission"]["localCaseId"],
        "submission_type": SubmissionType.initial,
        "submission_metadata": redact_metadata_dict(original),
    }
    return Submission(**(defaults | fields))


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


def test_collect_exports_restored_metadata():
    """A submission with stored metadata is exported, restored, with what the manifest needs."""
    entries, skipped = collect_export_entries([_submission("a1")], include_test_submissions=False)

    assert skipped == []
    [entry] = entries
    assert entry.submission_id == "260914050_2025-09-15_c64603a1"
    assert entry.content == _original()
    assert entry.unrestored == frozenset()
    assert entry.metadata_version == "1.3"
    assert entry.latest_state is None


def test_collect_skips_missing_metadata():
    """A submission without stored metadata is skipped, with the reason."""
    entries, skipped = collect_export_entries(
        [_submission("a2", submission_metadata=None)], include_test_submissions=False
    )

    assert entries == []
    [skip] = skipped
    assert skip.submission_id == "260914050_2025-09-15_c64603a2"
    assert skip.reason == "no metadata stored in the database"


def test_collect_skips_test_submissions_by_default():
    """Test submissions are left out unless asked for, and listed as skipped."""
    entries, skipped = collect_export_entries(
        [_submission("a3", submission_type=SubmissionType.test)], include_test_submissions=False
    )

    assert entries == []
    [skip] = skipped
    assert skip.reason == "test submission"


def test_collect_includes_test_submissions_on_request():
    """Test submissions are exported when asked for."""
    entries, skipped = collect_export_entries(
        [_submission("a3", submission_type=SubmissionType.test)], include_test_submissions=True
    )

    assert len(entries) == 1
    assert skipped == []


def test_collect_reports_unrestored_fields():
    """A submission whose tanG the database does not hold is exported with the field reported."""
    entries, _ = collect_export_entries([_submission("a4", tan_g=None)], include_test_submissions=False)

    [entry] = entries
    assert entry.content["submission"]["tanG"] == REDACTED_TAN
    assert entry.unrestored == frozenset({"tan_g"})


def test_schema_version_reads_any_version():
    """The version comes from each document's own URL, in both two- and three-part form."""
    url = "https://raw.githubusercontent.com/BfArM-MVH/MVGenomseq/refs/tags/{}/GRZ/grz-schema.json"
    assert _schema_version({"$schema": url.format("v1.2.1")}) == "1.2.1"
    assert _schema_version({"$schema": url.format("v1.1.9")}) == "1.1.9"


@pytest.mark.parametrize(
    "content",
    [
        {},
        {"$schema": ""},
        {"$schema": "https://example.org/schema.json"},
        {"$schema": None},
        {"$schema": 13},
    ],
    ids=["missing", "empty", "unknown-url", "null", "not-a-string"],
)
def test_schema_version_without_known_url(content: dict):
    """A missing, empty, unknown or non-string $schema gives no version instead of an error."""
    assert _schema_version(content) is None


def test_collect_exports_metadata_without_schema_url():
    """A document without a $schema URL is still exported, only without a version."""
    metadata = redact_metadata_dict(_original())
    del metadata["$schema"]

    entries, skipped = collect_export_entries(
        [_submission("a5", submission_metadata=metadata)], include_test_submissions=False
    )

    assert skipped == []
    [entry] = entries
    assert entry.metadata_version is None


CREATED_AT = datetime.datetime(2026, 10, 6, 9, 30, tzinfo=datetime.UTC)


def _export(tmp_path: Path, entries: list[ExportEntry], skipped=()) -> tuple[Path, dict]:
    """Write an export into *tmp_path*; returns the zip's path and the manifest."""
    output = tmp_path / "export.zip"
    manifest = write_metadata_zip(entries, list(skipped), output, include_test_submissions=False, created_at=CREATED_AT)
    return output, manifest


def test_zip_holds_each_metadata_and_the_manifest(tmp_path: Path):
    """The zip holds one metadata.json per exported submission, restored, plus the manifest."""
    entries, _ = collect_export_entries([_submission("b1"), _submission("b2")], include_test_submissions=False)

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
    entries, skipped = collect_export_entries(
        [_submission("c1", tan_g=None), _submission("c2", submission_metadata=None)],
        include_test_submissions=False,
    )

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
        "sha256": hashlib.sha256(exported_bytes).hexdigest(),
        "metadata_version": "1.3",
        "submission_uploaded_date": None,
        "latest_state": None,
        "unrestored_fields": ["tan_g"],
    }
    assert manifest["skipped"] == [
        {"submission_id": "260914050_2025-09-15_c64603c2", "reason": "no metadata stored in the database"}
    ]


def test_manifest_holds_no_tan_g_or_local_case_id(tmp_path: Path):
    """The manifest can be inspected without exposing what the export protects."""
    entries, _ = collect_export_entries([_submission("d1")], include_test_submissions=False)

    output, _ = _export(tmp_path, entries)

    with zipfile.ZipFile(output) as archive:
        manifest_text = archive.read("manifest.json").decode("utf-8")
    assert _original()["submission"]["tanG"] not in manifest_text
    assert _original()["submission"]["localCaseId"] not in manifest_text


def test_refuses_to_overwrite_an_existing_export(tmp_path: Path):
    """An existing file at the output path is left as it is."""
    entries, _ = collect_export_entries([_submission("e1")], include_test_submissions=False)
    output = tmp_path / "export.zip"
    output.write_text("an earlier export")

    with pytest.raises(FileExistsError):
        write_metadata_zip(entries, [], output, include_test_submissions=False)

    assert output.read_text() == "an earlier export"


def test_failed_export_leaves_nothing_behind(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """A failure halfway through writing leaves neither the zip nor its temporary file."""
    entries, _ = collect_export_entries([_submission("f1")], include_test_submissions=False)

    def disk_full(*args, **kwargs):
        raise OSError("No space left on device")

    monkeypatch.setattr(zipfile.ZipFile, "writestr", disk_full)

    with pytest.raises(OSError, match="No space left"):
        _export(tmp_path, entries)

    assert list(tmp_path.iterdir()) == []


UPLOAD_DATE = "2025-09-16"


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


def _export_args(output: Path, public_key_path: Path, *extra: str) -> list[str]:
    return ["db", "export-metadata", "--output", str(output), "--public-key", str(public_key_path), *extra]


def test_command_exports_what_was_submitted(
    migrated_database_config_path: Path, tmp_path: Path, recipient_key_paths: tuple[Path, Path]
):
    """The decrypted metadata.json equals what the submitter sent, tanG and localCaseId included."""
    public_key_path, private_key_path = recipient_key_paths
    first_id, first_submitted = _populate(
        migrated_database_config_path, tmp_path, submissionType="initial", tanG="b" * 64
    )
    second_id, second_submitted = _populate(
        migrated_database_config_path, tmp_path, submissionType="initial", tanG="c" * 64, localCaseId="case-2"
    )
    _invoke(migrated_database_config_path, "db", "submission", "update", first_id, "Finished")
    output = tmp_path / "export.zip.c4gh"

    result = _invoke(migrated_database_config_path, *_export_args(output, public_key_path))

    assert "Exported 2 metadata.json file(s)" in result.stdout
    assert hashlib.sha256(output.read_bytes()).hexdigest() in result.stdout
    with _decrypt(output, private_key_path) as archive:
        assert json.loads(archive.read(f"{first_id}/metadata.json")) == first_submitted
        assert json.loads(archive.read(f"{second_id}/metadata.json")) == second_submitted
        manifest = json.loads(archive.read("manifest.json"))
    described = {submission["submission_id"]: submission for submission in manifest["submissions"]}
    assert described[first_id]["latest_state"] == "Finished"
    assert described[second_id]["latest_state"] is None
    assert described[first_id]["submission_uploaded_date"] == UPLOAD_DATE
    assert described[first_id]["unrestored_fields"] == []
    assert manifest["skipped"] == []


def test_command_writes_nothing_readable_without_the_private_key(
    migrated_database_config_path: Path, tmp_path: Path, recipient_key_paths: tuple[Path, Path]
):
    """The export is no readable zip, and holds no tanG in plain text."""
    public_key_path, _ = recipient_key_paths
    _, submitted = _populate(migrated_database_config_path, tmp_path, submissionType="initial", tanG="d" * 64)
    output = tmp_path / "export.zip.c4gh"

    _invoke(migrated_database_config_path, *_export_args(output, public_key_path))

    assert not zipfile.is_zipfile(output)
    assert submitted["submission"]["tanG"].encode() not in output.read_bytes()
    assert list(tmp_path.glob(".export.zip.c4gh.*")) == []


def test_command_skips_test_and_unpopulated_submissions(
    migrated_database_config_path: Path, tmp_path: Path, recipient_key_paths: tuple[Path, Path]
):
    """Test submissions and submissions without metadata are listed as skipped, on stderr and in the manifest."""
    public_key_path, private_key_path = recipient_key_paths
    test_id, _ = _populate(migrated_database_config_path, tmp_path)  # the example metadata is a test submission
    unpopulated_id = "260914050_2025-09-15_0000abcd"
    _invoke(migrated_database_config_path, "db", "submission", "add", unpopulated_id)
    output = tmp_path / "export.zip.c4gh"

    result = _invoke(migrated_database_config_path, *_export_args(output, public_key_path))

    assert f"Skipped {test_id}: test submission." in result.stderr
    assert f"Skipped {unpopulated_id}: no metadata stored in the database." in result.stderr
    with _decrypt(output, private_key_path) as archive:
        assert archive.namelist() == ["manifest.json"]
        manifest = json.loads(archive.read("manifest.json"))
    assert {skip["submission_id"]: skip["reason"] for skip in manifest["skipped"]} == {
        test_id: "test submission",
        unpopulated_id: "no metadata stored in the database",
    }


def test_command_includes_test_submissions_on_request(
    migrated_database_config_path: Path, tmp_path: Path, recipient_key_paths: tuple[Path, Path]
):
    """--include-test-submissions exports test submissions too."""
    public_key_path, private_key_path = recipient_key_paths
    test_id, submitted = _populate(migrated_database_config_path, tmp_path)
    output = tmp_path / "export.zip.c4gh"

    _invoke(migrated_database_config_path, *_export_args(output, public_key_path, "--include-test-submissions"))

    with _decrypt(output, private_key_path) as archive:
        assert json.loads(archive.read(f"{test_id}/metadata.json")) == submitted
        assert json.loads(archive.read("manifest.json"))["export"]["test_submissions_included"] is True


def test_command_refuses_to_overwrite(
    migrated_database_config_path: Path, tmp_path: Path, recipient_key_paths: tuple[Path, Path]
):
    """An existing output file aborts the command and is left as it is."""
    public_key_path, _ = recipient_key_paths
    output = tmp_path / "export.zip.c4gh"
    output.write_text("an earlier export")

    result = CliRunner().invoke(
        grzctl.cli.build_cli(),
        ["--config", str(migrated_database_config_path), *_export_args(output, public_key_path)],
    )

    assert result.exit_code != 0
    assert "Refusing to overwrite" in result.stderr
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
    assert "cannot be read" in result.stderr
    assert not output.exists()
