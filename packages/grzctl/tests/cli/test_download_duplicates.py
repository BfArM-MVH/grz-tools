"""Tests for the early duplicate checks of ``grzctl download``.

A submission whose tanG is taken, or whose case already has a QC-passed initial submission, has to
fail from its metadata alone: the metadata is downloaded first, so no file may be transferred.
"""

import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import click.testing
import grzctl.cli
import pytest
import yaml
from grz_db.errors import DuplicateInitialSubmissionError, DuplicateTanGError
from grz_db.models.submission import FailureReasonEnum, SubmissionDb
from grz_pydantic_models.submission.metadata import GrzSubmissionMetadata
from grzctl.commands.duplicate_checks import reject_duplicates
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

        Returns the result and what the fake worker did: whether it reached the file transfer.
        """
        output_dir = tmp_path / f"out-{submission_id}"
        output_dir.mkdir()
        metadata = GrzSubmissionMetadata.model_validate(raw)
        seen = SimpleNamespace(files_downloaded=False, metadata_check=None)

        def fake_download(s3_options, sid, force=False, metadata_version_check=None, metadata_check=None):
            seen.metadata_check = metadata_check
            if metadata_check is not None:
                metadata_check(metadata)
            seen.files_downloaded = True

        with patch("grzctl.commands.download.Worker") as worker_cls:
            worker_cls.return_value.download.side_effect = fake_download
            result = invoke(
                "download", "--submission-id", submission_id, "--output-dir", str(output_dir),
                "--inbox", INBOX, "--no-populate", *args,
            )  # fmt: skip
        return SimpleNamespace(result=result, seen=seen)

    return SimpleNamespace(invoke=invoke, add=add, show=show, download=download)


def _error_state(shown: dict) -> dict:
    return [s for s in shown["states"] if s["state"] == "Error"][-1]


def test_download_rejects_a_taken_tan_g_before_any_file(harness, test_metadata_path):
    """A tanG held by another submission fails the download before the files are transferred."""
    holder = _metadata(test_metadata_path, tan_g="a" * 64, date="2025-01-01")
    newcomer = _metadata(test_metadata_path, tan_g="a" * 64, date="2025-01-02")
    holder_id = harness.add(holder, populate=True)
    newcomer_id = harness.add(newcomer, populate=False)

    outcome = harness.download(newcomer_id, newcomer)

    assert outcome.result.exit_code != 0
    assert not outcome.seen.files_downloaded
    error = _error_state(harness.show(newcomer_id))
    assert error["failure_reason"] == FailureReasonEnum.DUPLICATE_TANG.value
    assert holder_id in error["data"]["error"]


def test_download_rejects_a_duplicate_initial_before_any_file(harness, test_metadata_path):
    """A second initial of a case that has a QC-passed one fails basic QC before the files are transferred."""
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


def test_download_continues_when_nothing_is_taken(harness, test_metadata_path):
    """A submission with a fresh tanG, and a case without a QC-passed initial, downloads its files."""
    other = _metadata(test_metadata_path, tan_g="a" * 64, date="2025-01-01")
    other["submission"]["localCaseId"] = "another-case"
    harness.add(other, populate=True, basic_qc_passed=True)
    fresh = _metadata(test_metadata_path, tan_g="b" * 64, date="2025-01-02")
    fresh_id = harness.add(fresh, populate=False)

    outcome = harness.download(fresh_id, fresh)

    assert outcome.result.exit_code == 0, outcome.result.output
    assert outcome.seen.files_downloaded
    assert harness.show(fresh_id)["states"][-1]["state"] == "Downloaded"


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
    assert outcome.seen.metadata_check is None
    assert outcome.seen.files_downloaded


def test_reject_duplicates_names_the_holder_and_keeps_to_the_rules(
    harness, migrated_database_config, test_metadata_path
):
    """The checks themselves: tanG first, a second initial of the case, and a followup that joins it."""
    holder = _metadata(test_metadata_path, tan_g="a" * 64, date="2025-01-01")
    harness.add(holder, populate=True, basic_qc_passed=True)
    db = SubmissionDb(db_url=migrated_database_config.db.database_url, author=None)

    def submission(**kwargs) -> GrzSubmissionMetadata:
        raw = _metadata(test_metadata_path, **kwargs)
        harness.add(raw, populate=False)  # the download has added the submission to the database already
        return GrzSubmissionMetadata.model_validate(raw)

    taken_tan_g = submission(tan_g="a" * 64, date="2025-01-02")
    with pytest.raises(DuplicateTanGError, match="already used by submission"):
        reject_duplicates(db, taken_tan_g.submission_id, taken_tan_g)

    second_initial = submission(tan_g="b" * 64, date="2025-01-03")
    with pytest.raises(DuplicateInitialSubmissionError):
        reject_duplicates(db, second_initial.submission_id, second_initial)

    follow_up = submission(tan_g="c" * 64, date="2025-01-04", submission_type="followup")
    reject_duplicates(db, follow_up.submission_id, follow_up)  # a followup joins the case, whatever it holds
