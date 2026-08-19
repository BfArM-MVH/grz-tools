"""Tests for the grzctl ``validate`` command."""

import json
import logging
from types import SimpleNamespace
from unittest.mock import patch

import click.testing
import grzctl.cli
import pytest
from grz_common.workers.worker import Worker
from grz_db.models.submission import SubmissionDb
from grz_pydantic_models.submission.metadata import GrzSubmissionMetadata


@pytest.fixture
def grzctl_config_path(tmp_path):
    import yaml

    config = {
        "leistungserbringer": {"000000000": {"inbox_buckets": {"inbox": {"private_key_path": "/dev/null"}}}},
        "archives": {
            "consented": {"s3": {"bucket": "consented"}, "public_key_path": "/dev/null"},
            "non_consented": {"s3": {"bucket": "non_consented"}, "public_key_path": "/dev/null"},
        },
        "db": {"database_url": "sqlite:///:memory:", "author": {"name": "test"}},
        "pruefbericht": {},
        "keys": {"grz_private_key_path": "/dev/null"},
        "identifiers": {"grz": "GRZT00000"},
    }
    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(config, f)
    return config_path


@pytest.mark.parametrize(
    ("flag", "expected_no_mmap"),
    [
        ("--mmap", False),
        ("--no-mmap", True),
        (None, True),  # default is no-mmap
    ],
)
def test_validate_forwards_mmap_to_worker(tmp_path, flag, expected_no_mmap, grzctl_config_path):
    """The grzctl ``validate`` command must forward the ``--mmap/--no-mmap`` flag to
    ``Worker.validate`` (as the inverted ``no_mmap`` keyword).

    Regression test for a previous wrapper implementation which invoked grz-cli's
    inner callback directly and silently dropped misnamed keywords into ``**kwargs``.
    """
    submission_dir = tmp_path / "submission"
    for sub in ("metadata", "files", "logs"):
        (submission_dir / sub).mkdir(parents=True)

    with patch("grzctl.commands.validate.Worker") as mock_worker_cls:
        mock_worker = mock_worker_cls.return_value
        mock_worker.parse_submission.return_value.metadata.content.submission_id = "S1"

        args = [
            "--config",
            str(grzctl_config_path),
            "validate",
            "--submission-dir",
            str(submission_dir),
            "--submitter-id",
            "000000000",
            "--no-update-db",
        ]
        if flag is not None:
            args.append(flag)

        runner = click.testing.CliRunner()
        cli = grzctl.cli.build_cli()
        result = runner.invoke(cli, args)

        assert result.exit_code == 0, result.output
        mock_worker.validate.assert_called_once()
        assert mock_worker.validate.call_args.kwargs["no_mmap"] is expected_no_mmap


#: The submitterId of the shipped example metadata these tests build their submissions from.
#: ``validate`` requires it on the command line and checks the metadata against it.
SUBMITTER_ID = "260914050"


def _set_up_case_with_qc_passed_initial(tmp_path, migrated_database_config_path, test_metadata_path) -> SimpleNamespace:
    """Register two competing initial submissions for one case; one has already passed basic QC.

    Returns the CLI harness plus the directory of the other (duplicate) initial submission,
    ready for ``validate``.
    """
    qc_passed_raw = json.loads(test_metadata_path.read_text())
    qc_passed_raw["submission"]["submissionType"] = "initial"
    duplicate_raw = json.loads(json.dumps(qc_passed_raw))
    duplicate_raw["submission"]["submissionDate"] = "2025-01-02"  # distinct id, same local case id
    duplicate_raw["submission"]["tanG"] = "f" * 64  # a re-upload comes with a fresh transaction number
    qc_passed_id = GrzSubmissionMetadata.model_validate(qc_passed_raw).submission_id
    duplicate_id = GrzSubmissionMetadata.model_validate(duplicate_raw).submission_id

    qc_passed_meta = tmp_path / "qc_passed.metadata.json"
    qc_passed_meta.write_text(json.dumps(qc_passed_raw))
    duplicate_meta = tmp_path / "duplicate.metadata.json"
    duplicate_meta.write_text(json.dumps(duplicate_raw))
    submission_dir = tmp_path / "duplicate-submission"
    for sub in ("metadata", "files", "logs", "encrypted_files"):
        (submission_dir / sub).mkdir(parents=True)
    (submission_dir / "metadata" / "metadata.json").write_text(json.dumps(duplicate_raw))

    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    config = str(migrated_database_config_path)

    def db(*args: str) -> click.testing.Result:
        return runner.invoke(cli, ["--config", config, "db", *args])

    for sid, meta in ((qc_passed_id, qc_passed_meta), (duplicate_id, duplicate_meta)):
        result = db("submission", "add", sid)
        assert result.exit_code == 0, result.output
        result = db("submission", "populate", sid, str(meta), "--no-confirm")
        assert result.exit_code == 0, result.output
    result = db("submission", "modify", qc_passed_id, "basic_qc_passed", "true")
    assert result.exit_code == 0, result.output

    return SimpleNamespace(
        runner=runner,
        cli=cli,
        config=config,
        db=db,
        qc_passed_id=qc_passed_id,
        duplicate_id=duplicate_id,
        submission_dir=submission_dir,
    )


def _assert_duplicate_failed_basic_qc(competing: SimpleNamespace) -> None:
    """The duplicate initial submission must carry ``basic_qc_passed=False`` and an ERROR log
    naming the case's QC-passed initial submission, which itself must be left unaffected.
    """
    result = competing.db("submission", "show", competing.duplicate_id, "--json")
    assert result.exit_code == 0, result.output
    duplicate = json.loads(result.stdout)
    assert duplicate["basic_qc_passed"] is False
    error = [s for s in duplicate["states"] if s["state"] == "Error"][-1]
    assert error["failure_reason"] == "duplicate_initial"
    assert competing.qc_passed_id in error["data"]["error"]
    result = competing.db("submission", "show", competing.qc_passed_id, "--json")
    assert json.loads(result.stdout)["basic_qc_passed"] is True


def test_validate_rejects_duplicate_initial_before_validating(
    tmp_path, monkeypatch, migrated_database_config_path, test_metadata_path
):
    """A duplicate initial submission fails basic QC before validation runs.

    The case already has an initial submission that passed basic QC, so the DB pre-check
    rejects the duplicate initial submission: ``basic_qc_passed`` is set to ``False``, an
    ERROR state log carries the ``duplicate_initial`` failure reason, and the command exits
    with a nonzero code, without spending validation effort.
    """
    competing = _set_up_case_with_qc_passed_initial(tmp_path, migrated_database_config_path, test_metadata_path)

    def _must_not_validate(self, **kwargs):
        raise AssertionError("validation must not run for a duplicate initial submission")

    monkeypatch.setattr(Worker, "validate", _must_not_validate)
    result = competing.runner.invoke(
        competing.cli,
        [
            "--config",
            competing.config,
            "validate",
            "--submission-dir",
            str(competing.submission_dir),
            "--submitter-id",
            SUBMITTER_ID,
            "--update-db",
        ],
    )
    assert result.exit_code != 0, "a failed basic QC must fail the validate step"
    _assert_duplicate_failed_basic_qc(competing)


def test_validate_records_duplicate_initial_detected_during_validation(
    tmp_path, monkeypatch, migrated_database_config_path, test_metadata_path
):
    """The write-time check catches a competing initial submission that already passed basic
    QC, once a blind pre-check has missed it.

    The pre-check is forced inconclusive, so the duplicate rejection is triggered instead by
    the write-time update, which records basic QC as failed the same way the pre-check would
    have.
    """
    competing = _set_up_case_with_qc_passed_initial(tmp_path, migrated_database_config_path, test_metadata_path)

    monkeypatch.setattr(SubmissionDb, "resolve_case", lambda self, *args, **kwargs: None)
    monkeypatch.setattr(Worker, "validate", lambda self, **kwargs: None)
    result = competing.runner.invoke(
        competing.cli,
        [
            "--config",
            competing.config,
            "validate",
            "--submission-dir",
            str(competing.submission_dir),
            "--submitter-id",
            SUBMITTER_ID,
            "--update-db",
        ],
    )
    assert result.exit_code != 0, "a failed basic QC must fail the validate step"
    _assert_duplicate_failed_basic_qc(competing)


def test_validate_without_update_db_warns_about_duplicate_initial(
    tmp_path, monkeypatch, caplog, migrated_database_config_path, test_metadata_path
):
    """With ``--no-update-db``, a duplicate initial submission produces a warning.

    Basic QC is not marked failed and the database is not written.
    """
    competing = _set_up_case_with_qc_passed_initial(tmp_path, migrated_database_config_path, test_metadata_path)

    monkeypatch.setattr(Worker, "validate", lambda self, **kwargs: None)
    with caplog.at_level(logging.WARNING, logger="grzctl.commands.validate"):
        result = competing.runner.invoke(
            competing.cli,
            [
                "--config",
                competing.config,
                "validate",
                "--submission-dir",
                str(competing.submission_dir),
                "--submitter-id",
                SUBMITTER_ID,
                "--no-update-db",
            ],
        )
    assert result.exit_code == 0, result.output
    assert "already has an initial submission" in caplog.text

    result = competing.db("submission", "show", competing.duplicate_id, "--json")
    assert result.exit_code == 0, result.output
    assert json.loads(result.stdout)["basic_qc_passed"] is None


def _submission_dir_with_initial_metadata(tmp_path, test_metadata_path) -> str:
    """A submission directory holding one ``initial`` metadata.json, ready for ``validate``."""
    raw = json.loads(test_metadata_path.read_text())
    raw["submission"]["submissionType"] = "initial"
    submission_dir = tmp_path / "submission"
    for sub in ("metadata", "files", "logs", "encrypted_files"):
        (submission_dir / sub).mkdir(parents=True)
    (submission_dir / "metadata" / "metadata.json").write_text(json.dumps(raw))
    return str(submission_dir)


def _validate_without_update_db(config_path, submission_dir) -> click.testing.Result:
    runner = click.testing.CliRunner()
    return runner.invoke(
        grzctl.cli.build_cli(),
        [
            "--config",
            str(config_path),
            "validate",
            "--submission-dir",
            submission_dir,
            "--submitter-id",
            SUBMITTER_ID,
            "--no-update-db",
        ],
    )


def test_validate_without_update_db_reports_a_check_it_could_not_run(
    tmp_path, monkeypatch, caplog, empty_database_config_path, test_metadata_path
):
    """A configured database that cannot answer must be reported, not passed over in silence.

    Nothing is written in this mode, so the command still succeeds. What must not happen is
    the operator reading the absence of a warning as "not a duplicate": here the database is
    configured but behind on its schema, so the question was never put to it.
    """
    submission_dir = _submission_dir_with_initial_metadata(tmp_path, test_metadata_path)
    monkeypatch.setattr(Worker, "validate", lambda self, **kwargs: None)

    with caplog.at_level(logging.WARNING, logger="grzctl.commands.validate"):
        result = _validate_without_update_db(empty_database_config_path, submission_dir)

    assert result.exit_code == 0, result.output
    warnings = [r.message for r in caplog.records if r.name == "grzctl.commands.validate"]
    assert len(warnings) == 1, warnings
    assert "Could not check whether this is a duplicate initial submission" in warnings[0]
    assert "not at latest schema" in warnings[0], "the reason the check could not run must be named"
