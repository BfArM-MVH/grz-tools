import shutil
from unittest.mock import MagicMock, patch

import click.testing
import pytest
import sqlalchemy
import yaml
from grz_db.models.submission import Submission, SubmissionStateEnum
from grzctl.cli import build_cli
from grzctl.models.config import DbConfig
from sqlalchemy.orm import selectinload
from sqlmodel import Session, select


@pytest.fixture
def full_config_path(
    tmp_path,
    db_config_content,
    s3_config_content,
    keys_config_content,
    pruefbericht_config_content,
):
    config_data = {}
    config_data.update(db_config_content)
    config_data.update(s3_config_content)
    config_data.update(keys_config_content)
    config_data.update(pruefbericht_config_content)

    if "author" in config_data.get("db", {}):
        config_data["db"]["author"]["private_key_passphrase"] = "test"

    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(config_data, f)

    runner = click.testing.CliRunner()
    cli = build_cli()
    result = runner.invoke(cli, ["db", "--config-file", str(config_path), "init"])
    assert result.exit_code == 0, f"DB Init failed: {result.output}"

    return config_path


@pytest.fixture
def db_engine(full_config_path):
    config = DbConfig.from_path(full_config_path)
    return sqlalchemy.create_engine(config.db.database_url)


@pytest.fixture
def test_metadata(tmp_path, submission_metadata):
    parsed_metadata = submission_metadata.content

    metadata_path = tmp_path / "metadata.json"
    metadata_path.write_text(parsed_metadata.model_dump_json(by_alias=True))

    return parsed_metadata, metadata_path


def setup_db_state(
    runner, cli, config_file, submission_id, metadata_path, initial_state: SubmissionStateEnum | None = None
):
    result = runner.invoke(cli, ["db", "--config-file", str(config_file), "submission", "add", submission_id])
    assert result.exit_code == 0, f"Setup add failed: {result.output}"

    if metadata_path and metadata_path.exists():
        result = runner.invoke(
            cli,
            [
                "db",
                "--config-file",
                str(config_file),
                "submission",
                "populate",
                submission_id,
                str(metadata_path),
                "--no-confirm",
            ],
        )
        assert result.exit_code == 0, f"Setup populate failed: {result.output}"

    if initial_state:
        result = runner.invoke(
            cli, ["db", "--config-file", str(config_file), "submission", "update", submission_id, initial_state.value]
        )
        assert result.exit_code == 0, f"Setup update state failed: {result.output}"


def get_state_history(engine, submission_id) -> list[SubmissionStateEnum]:
    with Session(engine) as session:
        statement = select(Submission).where(Submission.id == submission_id).options(selectinload(Submission.states))
        submission = session.exec(statement).first()
        if not submission:
            return []
        return [log.state for log in sorted(submission.states, key=lambda x: x.id)]


@pytest.mark.parametrize(
    "command_spec",
    [
        {
            "cmd": ["archive", "--submission-dir", "SUBMISSION_DIR"],
            "worker_patch": "grzctl.commands.archive.Worker",
            "id_source": "encrypted_submission",
            "initial_state": SubmissionStateEnum.ENCRYPTED,
            "intermediate_state": SubmissionStateEnum.ARCHIVING,
            "expected_state": SubmissionStateEnum.ARCHIVED,
        },
        {
            "cmd": ["download", "--submission-id", "SUBMISSION_ID", "--output-dir", "OUTPUT_DIR"],
            "worker_patch": "grzctl.commands.download.Worker",
            "id_source": "arg",
            "initial_state": SubmissionStateEnum.UPLOADED,
            "intermediate_state": SubmissionStateEnum.DOWNLOADING,
            "expected_state": SubmissionStateEnum.DOWNLOADED,
            "skip_populate": True,
        },
        {
            "cmd": ["decrypt", "--submission-dir", "SUBMISSION_DIR"],
            "worker_patch": "grzctl.commands.decrypt.Worker",
            "id_source": "encrypted_submission",
            "initial_state": SubmissionStateEnum.DOWNLOADED,
            "intermediate_state": SubmissionStateEnum.DECRYPTING,
            "expected_state": SubmissionStateEnum.DECRYPTED,
        },
        {
            "cmd": ["validate", "--submission-dir", "SUBMISSION_DIR"],
            "worker_patch": "grzctl.commands.validate.Worker",
            "extra_patch": "grzctl.commands.validate.validate_module.validate.callback",
            "id_source": "submission",
            "initial_state": SubmissionStateEnum.DECRYPTED,
            "intermediate_state": SubmissionStateEnum.VALIDATING,
            "expected_state": SubmissionStateEnum.VALIDATED,
        },
        {
            "cmd": ["encrypt", "--submission-dir", "SUBMISSION_DIR"],
            "worker_patch": "grzctl.commands.encrypt.Worker",
            "extra_patch": "grzctl.commands.encrypt.encrypt_module.encrypt.callback",
            "id_source": "submission",
            "initial_state": SubmissionStateEnum.VALIDATED,
            "intermediate_state": SubmissionStateEnum.ENCRYPTING,
            "expected_state": SubmissionStateEnum.ENCRYPTED,
        },
        {
            "cmd": ["clean", "--submission-id", "SUBMISSION_ID", "--yes-i-really-mean-it"],
            "extra_patch": "grzctl.commands.clean._clean_submission_from_bucket",
            "id_source": "arg",
            "initial_state": SubmissionStateEnum.QCED,
            "intermediate_state": SubmissionStateEnum.CLEANING,
            "expected_state": SubmissionStateEnum.CLEANED,
            "skip_populate": True,
        },
    ],
)
def test_db_wrappers(  # noqa: C901
    command_spec,
    db_engine,
    full_config_path,
    test_metadata,
    tmp_path,
):
    runner = click.testing.CliRunner()
    cli = build_cli()

    parsed_metadata, metadata_path_fixture = test_metadata
    submission_id = parsed_metadata.submission_id

    submission_dir = tmp_path / "submission"
    submission_dir.mkdir()
    for d in ["metadata", "files", "logs", "encrypted_files"]:
        (submission_dir / d).mkdir()
    output_dir = tmp_path / "output"

    if not command_spec.get("skip_populate", False):
        shutil.copy(metadata_path_fixture, submission_dir / "metadata" / "metadata.json")

    setup_db_state(
        runner,
        cli,
        full_config_path,
        submission_id,
        metadata_path_fixture if not command_spec.get("skip_populate") else None,
        initial_state=command_spec["initial_state"],
    )

    args = []
    for arg in command_spec["cmd"]:
        if arg == "SUBMISSION_DIR":
            args.append(str(submission_dir))
        elif arg == "OUTPUT_DIR":
            args.append(str(output_dir))
        elif arg == "SUBMISSION_ID":
            args.append(submission_id)
        else:
            args.append(arg)

    args.extend(["--config-file", str(full_config_path), "--update-db"])

    worker_patcher = None
    extra_patcher = None
    mock_worker = None
    mock_extra = None

    try:
        if "worker_patch" in command_spec:
            worker_patcher = patch(command_spec["worker_patch"])
            mock_worker_cls = worker_patcher.start()
            mock_worker = mock_worker_cls.return_value

            mock_submission = MagicMock()
            mock_submission.metadata.content.submission_id = submission_id
            mock_submission.submission_id = submission_id

            if command_spec["id_source"] == "submission":
                mock_worker.parse_submission.return_value = mock_submission
            elif command_spec["id_source"] == "encrypted_submission":
                mock_worker.parse_encrypted_submission.return_value = mock_submission

        if "extra_patch" in command_spec:
            extra_patcher = patch(command_spec["extra_patch"])
            mock_extra = extra_patcher.start()

        result = runner.invoke(cli, args)

        assert result.exit_code == 0, f"Command failed: {result.output}"

        if mock_extra:
            mock_extra.assert_called_once()
        elif mock_worker:
            method_name = command_spec["cmd"][0]
            getattr(mock_worker, method_name).assert_called_once()

        history = get_state_history(db_engine, submission_id)
        assert len(history) >= 2, f"History too short: {history}"
        assert history[-1] == command_spec["expected_state"]
        assert history[-2] == command_spec["intermediate_state"]
        if command_spec["initial_state"]:
            assert history[-3] == command_spec["initial_state"]

    finally:
        if extra_patcher:
            extra_patcher.stop()
        if worker_patcher:
            worker_patcher.stop()


def test_pruefbericht_wrapper(db_engine, full_config_path, test_metadata, tmp_path):
    runner = click.testing.CliRunner()
    cli = build_cli()

    parsed_metadata, metadata_path = test_metadata
    submission_id = parsed_metadata.submission_id

    setup_db_state(
        runner, cli, full_config_path, submission_id, metadata_path, initial_state=SubmissionStateEnum.ARCHIVED
    )

    with patch("grzctl.commands.pruefbericht._try_submit") as mock_submit:
        mock_submit.return_value = (None, "token")

        pb_path = tmp_path / "pruefbericht.json"
        pb_path.write_text("{}")

        with patch("grzctl.commands.pruefbericht.Pruefbericht") as MockPbModel:  # noqa: N806
            mock_pb = MockPbModel.model_validate_json.return_value
            mock_pb.submitted_case.tan = parsed_metadata.submission.tan_g

            args = [
                "pruefbericht",
                "submit",
                "--submission-id",
                submission_id,
                "--pruefbericht-file",
                str(pb_path),
                "--config-file",
                str(full_config_path),
                "--update-db",
            ]

            result = runner.invoke(cli, args)

            assert result.exit_code == 0, result.output
            mock_submit.assert_called_once()

            history = get_state_history(db_engine, submission_id)
            assert history[-1] == SubmissionStateEnum.REPORTED
            assert history[-2] == SubmissionStateEnum.REPORTING
            assert history[-3] == SubmissionStateEnum.ARCHIVED


def test_dbcontext_error_handling(db_engine, full_config_path, test_metadata, tmp_path):
    runner = click.testing.CliRunner()
    cli = build_cli()

    parsed_metadata, metadata_path = test_metadata
    submission_id = parsed_metadata.submission_id

    setup_db_state(runner, cli, full_config_path, submission_id, metadata_path, initial_state=SubmissionStateEnum.QCED)

    with patch("grzctl.commands.clean._clean_submission_from_bucket") as mock_clean:
        mock_clean.side_effect = RuntimeError("S3 Failure")

        args = [
            "clean",
            "--submission-id",
            submission_id,
            "--yes-i-really-mean-it",
            "--config-file",
            str(full_config_path),
            "--update-db",
        ]

        result = runner.invoke(cli, args)

        assert result.exit_code != 0

        history = get_state_history(db_engine, submission_id)
        assert history[-1] == SubmissionStateEnum.ERROR
        assert history[-2] == SubmissionStateEnum.CLEANING
        assert history[-3] == SubmissionStateEnum.QCED
