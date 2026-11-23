"""Tests for the early duplicate checks of ``grzctl download``.

The metadata is downloaded first. It is stored (``--populate``) and checked before any file is
transferred, so a submission whose tanG is taken, or whose case already has a QC-passed initial
submission, fails without its files being transferred and keeps the metadata it was rejected with.
"""

import datetime
import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import click.testing
import grzctl.cli
import pytest
import yaml
from grz_db.errors import AmbiguousCaseError
from grz_db.models.submission import FailureReasonEnum, SubmissionDb
from grz_pydantic_models.submission.metadata import GrzSubmissionMetadata
from grzctl.models.config import GrzctlConfig

SUBMITTER_ID = "260914050"
INBOX = "inbox"
UPLOAD_DATE = "2026-01-01"


def _metadata(test_metadata_path: Path, *, tan_g: str, date: str, submission_type: str = "initial") -> dict:
    """The test metadata as an initial submission of its case, with the given tanG and submission date."""
    raw = json.loads(test_metadata_path.read_text())
    raw["submission"]["submissionType"] = submission_type
    raw["submission"]["tanG"] = tan_g
    raw["submission"]["submissionDate"] = date
    return raw


def _id(raw: dict) -> str:
    return GrzSubmissionMetadata.model_validate(raw).submission_id


@pytest.fixture
def harness(tmp_path: Path, migrated_database_config: GrzctlConfig) -> SimpleNamespace:
    """The CLI on a migrated database, with the inbox of the test submitter configured."""
    data = migrated_database_config.model_dump(mode="json", exclude_none=True, context={"reveal_secrets": True})
    key_path = data["db"]["author"]["private_key_path"]
    data["leistungserbringer"] = {SUBMITTER_ID: {"inbox_buckets": {INBOX: {"private_key_path": str(key_path)}}}}
    config_path = tmp_path / "config.inbox.yaml"
    config_path.write_text(yaml.safe_dump(data))

    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()

    def invoke(*args: str) -> click.testing.Result:
        return runner.invoke(cli, ["--config", str(config_path), *args])

    def add(raw: dict, *, populate: bool, basic_qc_passed: bool = False) -> str:
        submission_id = _id(raw)
        metadata_file = tmp_path / f"{submission_id}.metadata.json"
        metadata_file.write_text(json.dumps(raw))
        result = invoke("db", "submission", "add", submission_id)
        assert result.exit_code == 0, result.output
        if populate:
            result = invoke(
                "db", "submission", "populate", submission_id, str(metadata_file),
                "--no-confirm", "--submission-date", UPLOAD_DATE,
            )  # fmt: skip
            assert result.exit_code == 0, result.output
        if basic_qc_passed:
            result = invoke("db", "submission", "modify", submission_id, "basic_qc_passed", "true")
            assert result.exit_code == 0, result.output
        return submission_id

    def show(submission_id: str) -> dict:
        result = invoke("db", "submission", "show", submission_id, "--json")
        assert result.exit_code == 0, result.output
        return json.loads(result.stdout)

    def download(submission_id: str, raw: dict, *args: str) -> SimpleNamespace:
        """Run ``download`` against a Worker that serves *raw* as the downloaded metadata.

        The S3 timestamp that populating needs is a fixed date, and the grz-cli version and the upload timing
        that download reads from the inbox are fixed too. Returns the result and what the fake
        worker did: whether it reached the file transfer, and the hook it was given.
        """
        output_dir = tmp_path / f"out-{submission_id}"
        output_dir.mkdir()
        metadata = GrzSubmissionMetadata.model_validate(raw)
        seen = SimpleNamespace(files_downloaded=False, on_metadata=None)

        def fake_download(s3_options, sid, force=False, metadata_version_check=None, on_metadata=None):
            seen.on_metadata = on_metadata
            if on_metadata is not None:
                on_metadata(metadata)
            seen.files_downloaded = True

        with (
            patch("grzctl.commands.download.Worker") as worker_cls,
            patch("grzctl.commands.download.init_s3_client"),
            patch("grzctl.commands.download.get_grzcli_version", return_value="3.0.0"),
            patch("grzctl.commands.download.get_upload_timing", return_value=None),
            patch(
                "grzctl.commands.download.get_metadata_upload_timestamp",
                return_value=datetime.datetime.fromisoformat(f"{UPLOAD_DATE}T12:00:00+00:00"),
            ),
        ):
            worker_cls.return_value.download.side_effect = fake_download
            result = invoke(
                "download", "--submission-id", submission_id, "--output-dir", str(output_dir),
                "--inbox", INBOX, *args,
            )  # fmt: skip
        return SimpleNamespace(result=result, seen=seen)

    return SimpleNamespace(invoke=invoke, add=add, show=show, download=download)


def _error_state(shown: dict) -> dict:
    return [s for s in shown["states"] if s["state"] == "Error"][-1]


def test_download_rejects_a_taken_tan_g_before_any_file(harness, test_metadata_path):
    """A tanG held by another submission fails while populating, before the files are transferred."""
    holder = _metadata(test_metadata_path, tan_g="a" * 64, date="2025-01-01")
    newcomer = _metadata(test_metadata_path, tan_g="a" * 64, date="2025-01-02")
    harness.add(holder, populate=True)
    newcomer_id = harness.add(newcomer, populate=False)

    outcome = harness.download(newcomer_id, newcomer)

    assert outcome.result.exit_code != 0
    assert not outcome.seen.files_downloaded
    error = _error_state(harness.show(newcomer_id))
    assert error["failure_reason"] == FailureReasonEnum.DUPLICATE_TANG.value


def test_download_rejects_a_duplicate_initial_before_any_file_and_keeps_its_metadata(harness, test_metadata_path):
    """A second initial of a case that has a QC-passed one fails basic QC, and has been populated by then."""
    first = _metadata(test_metadata_path, tan_g="a" * 64, date="2025-01-01")
    second = _metadata(test_metadata_path, tan_g="b" * 64, date="2025-01-02")  # same local case id, fresh tanG
    first_id = harness.add(first, populate=True, basic_qc_passed=True)
    second_id = harness.add(second, populate=False)

    outcome = harness.download(second_id, second)

    assert outcome.result.exit_code != 0
    assert not outcome.seen.files_downloaded
    shown = harness.show(second_id)
    error = _error_state(shown)
    assert error["failure_reason"] == FailureReasonEnum.DUPLICATE_INITIAL.value
    assert first_id in error["data"]["error"]
    assert shown["basic_qc_passed"] is False, "recorded as failing basic QC, as validate would"
    assert harness.show(first_id)["basic_qc_passed"] is True, "the case's QC-passed initial is left alone"
    # what populating stores, and the reports read: without it the rejected submission would be invisible to them
    assert shown["submission_type"] == "initial"
    assert shown["submission_uploaded_date"] == UPLOAD_DATE
    assert shown["tan_g"] == "b" * 64
    assert shown["submission_metadata"] is not None


def test_download_rejects_a_duplicate_initial_without_populating_when_told_not_to(harness, test_metadata_path):
    """With --no-populate the check still runs, and nothing is stored beyond what it records."""
    first = _metadata(test_metadata_path, tan_g="a" * 64, date="2025-01-01")
    second = _metadata(test_metadata_path, tan_g="b" * 64, date="2025-01-02")
    harness.add(first, populate=True, basic_qc_passed=True)
    second_id = harness.add(second, populate=True)  # linked to the case by an earlier populate

    outcome = harness.download(second_id, second, "--no-populate")

    assert outcome.result.exit_code != 0
    assert not outcome.seen.files_downloaded
    assert _error_state(harness.show(second_id))["failure_reason"] == FailureReasonEnum.DUPLICATE_INITIAL.value


def test_download_continues_when_nothing_is_taken_and_populates(harness, test_metadata_path):
    """A submission with a fresh tanG, and a case without a QC-passed initial, is stored and downloads its files."""
    other = _metadata(test_metadata_path, tan_g="a" * 64, date="2025-01-01")
    other["submission"]["localCaseId"] = "another-case"
    harness.add(other, populate=True, basic_qc_passed=True)
    fresh = _metadata(test_metadata_path, tan_g="b" * 64, date="2025-01-02")
    fresh_id = harness.add(fresh, populate=False)

    outcome = harness.download(fresh_id, fresh)

    assert outcome.result.exit_code == 0, outcome.result.output
    assert outcome.seen.files_downloaded
    shown = harness.show(fresh_id)
    assert shown["states"][-1]["state"] == "Downloaded"
    assert shown["tan_g"] == "b" * 64


def test_download_accepts_the_same_submission_again(harness, test_metadata_path):
    """Downloading a submission whose metadata is already stored does not collide with its own row."""
    raw = _metadata(test_metadata_path, tan_g="a" * 64, date="2025-01-01")
    submission_id = harness.add(raw, populate=True, basic_qc_passed=True)

    outcome = harness.download(submission_id, raw)

    assert outcome.result.exit_code == 0, outcome.result.output
    assert outcome.seen.files_downloaded


def test_download_without_update_db_makes_no_check(harness, test_metadata_path):
    """With --no-update-db the database is not opened, so there is nothing to check against."""
    raw = _metadata(test_metadata_path, tan_g="a" * 64, date="2025-01-01")
    submission_id = _id(raw)

    outcome = harness.download(submission_id, raw, "--no-update-db")

    assert outcome.result.exit_code == 0, outcome.result.output
    assert outcome.seen.on_metadata is None
    assert outcome.seen.files_downloaded


def test_download_continues_when_the_case_is_ambiguous(harness, test_metadata_path):
    """A case key that matches several cases does not stop the download, as it does not stop populate."""
    raw = _metadata(test_metadata_path, tan_g="a" * 64, date="2025-01-01")
    submission_id = harness.add(raw, populate=False)

    with patch.object(SubmissionDb, "assert_no_duplicate_initial", side_effect=AmbiguousCaseError("x", "y", [1, 2])):
        outcome = harness.download(submission_id, raw)

    assert outcome.result.exit_code == 0, outcome.result.output
    assert outcome.seen.files_downloaded
    assert harness.show(submission_id)["states"][-1]["state"] == "Downloaded"


def test_a_failure_to_record_the_failed_basic_qc_does_not_replace_the_duplicate(harness, test_metadata_path):
    """The failure reason is the duplicate, even if recording the failed basic QC itself fails."""
    first = _metadata(test_metadata_path, tan_g="a" * 64, date="2025-01-01")
    second = _metadata(test_metadata_path, tan_g="b" * 64, date="2025-01-02")
    harness.add(first, populate=True, basic_qc_passed=True)
    second_id = harness.add(second, populate=False)

    modify_submission = SubmissionDb.modify_submission

    def fail_recording_the_failed_basic_qc(self, submission_id, key, value, *args, **kwargs):
        if (key, value) == ("basic_qc_passed", "false"):
            raise RuntimeError("database is gone")
        return modify_submission(self, submission_id, key, value, *args, **kwargs)

    with patch.object(SubmissionDb, "modify_submission", fail_recording_the_failed_basic_qc):
        outcome = harness.download(second_id, second)

    assert outcome.result.exit_code != 0
    assert not outcome.seen.files_downloaded
    error = _error_state(harness.show(second_id))
    assert error["failure_reason"] == FailureReasonEnum.DUPLICATE_INITIAL.value
    assert "database is gone" not in error["data"]["error"]
