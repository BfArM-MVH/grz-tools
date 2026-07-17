import contextlib
import datetime
import shutil
from unittest.mock import MagicMock, patch

import click.testing
import pytest
import sqlalchemy
import yaml
from grz_db.errors import SubmissionNotFoundError
from grz_db.models.submission import Submission, SubmissionStateEnum
from grzctl.cli import build_cli
from grzctl.models.config import GrzctlConfig
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
    # Build a valid GrzctlConfig structure
    from tests.conftest import _grzctl_archives

    config_data = {
        "s3": {
            "inboxes": {
                "260914050": {
                    "inbox": {
                        "endpoint_url": "http://localhost:9000",
                        "private_key_path": "/dev/null",
                    }
                }
            }
        },
        "archives": _grzctl_archives(endpoint_url="http://localhost:9000"),
        "pruefbericht": {
            "authorization_url": "https://bfarm.localhost/token",
            "api_base_url": "https://bfarm.localhost/api/",
            "client_id": "pytest",
            "client_secret": "pysecret",
        },
        "identifiers": {"grz": "GRZK00007"},
    }
    config_data.update(db_config_content)
    config_data.update(keys_config_content)

    if "author" in config_data.get("db", {}):
        config_data["db"]["author"]["private_key_passphrase"] = "test"

    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(config_data, f)

    runner = click.testing.CliRunner()
    cli = build_cli()
    result = runner.invoke(cli, ["--config", str(config_path), "db", "init"])
    assert result.exit_code == 0, f"DB Init failed: {result.output}"

    return config_path


@pytest.fixture
def db_engine(full_config_path):
    config = GrzctlConfig.from_path(full_config_path)
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
    result = runner.invoke(cli, ["--config", str(config_file), "db", "submission", "add", submission_id])
    assert result.exit_code == 0, f"Setup add failed: {result.output}"

    if metadata_path and metadata_path.exists():
        result = runner.invoke(
            cli,
            [
                "--config",
                str(config_file),
                "db",
                "submission",
                "populate",
                submission_id,
                str(metadata_path),
                "--no-confirm",
            ],
            catch_exceptions=False,
        )
        assert result.exit_code == 0, f"Setup populate failed: {result.output}"

    if initial_state:
        result = runner.invoke(
            cli,
            ["--config", str(config_file), "db", "submission", "update", submission_id, initial_state.value],
            catch_exceptions=False,
        )
        assert result.exit_code == 0, f"Setup update state failed: {result.output}"


def get_state_history(engine, submission_id) -> list[SubmissionStateEnum]:
    with Session(engine) as session:
        statement = select(Submission).where(Submission.id == submission_id).options(selectinload(Submission.states))
        submission = session.exec(statement).first()
        if not submission:
            return []
        return [log.state for log in sorted(submission.states, key=lambda x: x.id)]


@contextlib.contextmanager
def mock_command(command_spec, submission_id):
    """
    Patch the Worker class and/or extra callable described by *command_spec*, then yield
    ``(mock_worker_instance, mock_extra)`` so callers can make assertions on them.

    Uses ExitStack so every patch is registered individually and all are torn down
    automatically when the ``with`` block exits, regardless of exceptions.
    """
    with contextlib.ExitStack() as stack:
        mock_worker = None
        mock_extra = None

        if "worker_patch" in command_spec:
            # Patch the Worker *class* so that Worker(...) returns our mock instance.
            mock_worker_cls = stack.enter_context(patch(command_spec["worker_patch"]))
            mock_worker = mock_worker_cls.return_value

            # Whichever parse_* method the command calls must return a submission
            # whose ID matches the one we put in the DB.
            mock_submission = MagicMock()
            mock_submission.metadata.content.submission_id = submission_id
            mock_submission.submission_id = submission_id

            if command_spec["id_source"] == "submission":
                mock_worker.parse_submission.return_value = mock_submission
            elif command_spec["id_source"] == "encrypted_submission":
                mock_worker.parse_encrypted_submission.return_value = mock_submission

        if "extra_patch" in command_spec:
            # Some commands invoke a standalone callable (e.g. validate.callback,
            # _clean_submission_from_bucket) instead of a Worker method.
            mock_extra = stack.enter_context(patch(command_spec["extra_patch"]))

        yield mock_worker, mock_extra


def build_args(
    command_spec,
    submission_dir=None,
    output_dir=None,
    submission_id=None,
    config_path=None,
    extra_flags=(),
):
    """Translate a command_spec cmd list (with placeholders) into a concrete argv list."""
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
    if config_path:
        args.insert(0, str(config_path))
        args.insert(0, "--config")
    args.extend(extra_flags)
    return args


@pytest.mark.parametrize(
    "command_spec",
    [
        {
            "cmd": ["download", "--submission-id", "SUBMISSION_ID", "--output-dir", "OUTPUT_DIR"],
            "worker_patch": "grzctl.commands.download.Worker",
            "id_source": "arg",
            "initial_state": None,
            "intermediate_state": SubmissionStateEnum.DOWNLOADING,
            "expected_state": SubmissionStateEnum.DOWNLOADED,
            "skip_populate": True,
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
            "worker_patch": "grzctl.commands.cli_wrappers.Worker",
            "extra_patch": "grzctl.commands.cli_wrappers.validate_module.validate.callback",
            "id_source": "submission",
            "initial_state": SubmissionStateEnum.DECRYPTED,
            "intermediate_state": SubmissionStateEnum.VALIDATING,
            "expected_state": SubmissionStateEnum.VALIDATED,
        },
        {
            "cmd": ["encrypt", "--submission-dir", "SUBMISSION_DIR"],
            "worker_patch": "grzctl.commands.cli_wrappers.Worker",
            "extra_patch": "grzctl.commands.cli_wrappers.encrypt_module.encrypt.callback",
            "id_source": "submission",
            "initial_state": SubmissionStateEnum.VALIDATED,
            "intermediate_state": SubmissionStateEnum.ENCRYPTING,
            "expected_state": SubmissionStateEnum.ENCRYPTED,
        },
        {
            "cmd": ["archive", "--submission-dir", "SUBMISSION_DIR"],
            "worker_patch": "grzctl.commands.archive.Worker",
            "id_source": "encrypted_submission",
            "initial_state": SubmissionStateEnum.ENCRYPTED,
            "intermediate_state": SubmissionStateEnum.ARCHIVING,
            "expected_state": SubmissionStateEnum.ARCHIVED,
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
def test_db_wrappers(
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

    args = build_args(
        command_spec,
        submission_dir=submission_dir,
        output_dir=output_dir,
        submission_id=submission_id,
        config_path=full_config_path,
        extra_flags=["--update-db"],
    )

    with mock_command(command_spec, submission_id) as (mock_worker, mock_extra):
        # download.py's --populate path calls `worker.parse_submission().metadata.content`
        # plus `get_metadata_upload_timestamp(...)`. Configure the mocked Worker and patch
        # the S3 helpers so the populate runs end-to-end against the real db. Harmless for
        # commands that never call these.
        if mock_worker is not None:
            mock_worker.parse_submission.return_value.metadata.content = parsed_metadata
        with (
            patch("grzctl.commands.download.init_s3_client") as mock_init_s3,
            patch("grzctl.commands.download.get_metadata_upload_timestamp") as mock_get_ts,
        ):
            mock_init_s3.return_value = MagicMock()
            mock_get_ts.return_value.date.return_value = datetime.date.today()
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


@pytest.mark.parametrize(
    "command_spec",
    [
        {
            "cmd": ["decrypt", "--submission-dir", "SUBMISSION_DIR"],
            "worker_patch": "grzctl.commands.decrypt.Worker",
            "id_source": "encrypted_submission",
        },
        {
            "cmd": ["validate", "--submission-dir", "SUBMISSION_DIR"],
            "worker_patch": "grzctl.commands.cli_wrappers.Worker",
            "extra_patch": "grzctl.commands.cli_wrappers.validate_module.validate.callback",
            "id_source": "submission",
        },
        {
            "cmd": ["encrypt", "--submission-dir", "SUBMISSION_DIR"],
            "worker_patch": "grzctl.commands.cli_wrappers.Worker",
            "extra_patch": "grzctl.commands.cli_wrappers.encrypt_module.encrypt.callback",
            "id_source": "submission",
        },
        {
            "cmd": ["archive", "--submission-dir", "SUBMISSION_DIR"],
            "worker_patch": "grzctl.commands.archive.Worker",
            "id_source": "encrypted_submission",
        },
        {
            "cmd": ["clean", "--submission-id", "SUBMISSION_ID", "--yes-i-really-mean-it"],
            "extra_patch": "grzctl.commands.clean._clean_submission_from_bucket",
            "id_source": "arg",
        },
    ],
)
def test_db_wrappers_submission_not_in_db(
    command_spec,
    db_engine,
    full_config_path,
    test_metadata,
    tmp_path,
):
    """Test that commands fail with SubmissionNotFoundError when the submission has not been added to the DB yet."""
    runner = click.testing.CliRunner()
    cli = build_cli()

    parsed_metadata, _ = test_metadata
    submission_id = parsed_metadata.submission_id

    submission_dir = tmp_path / "submission"
    submission_dir.mkdir()
    for d in ["metadata", "files", "logs", "encrypted_files"]:
        (submission_dir / d).mkdir()
    output_dir = tmp_path / "output"

    # Deliberately skip adding the submission to the DB

    args = build_args(
        command_spec,
        submission_dir=submission_dir,
        output_dir=output_dir,
        submission_id=submission_id,
        config_path=full_config_path,
        extra_flags=["--update-db"],
    )

    with mock_command(command_spec, submission_id):
        result = runner.invoke(cli, args)

        assert result.exit_code != 0, (
            f"Expected command '{command_spec['cmd'][0]}' to fail when submission is not in DB, "
            f"but it succeeded. Output: {result.output}"
        )
        assert isinstance(result.exception, SubmissionNotFoundError), (
            f"Expected SubmissionNotFoundError, got {type(result.exception)}: {result.exception}"
        )


@pytest.mark.parametrize(
    "command_spec",
    [
        {
            "cmd": ["decrypt", "--submission-dir", "SUBMISSION_DIR"],
            "worker_patch": "grzctl.commands.decrypt.Worker",
            "id_source": "encrypted_submission",
            "wrong_state": SubmissionStateEnum.ENCRYPTED,  # expected: DOWNLOADED
            "intermediate_state": SubmissionStateEnum.DECRYPTING,
            "expected_state": SubmissionStateEnum.DECRYPTED,
        },
        {
            "cmd": ["validate", "--submission-dir", "SUBMISSION_DIR"],
            "worker_patch": "grzctl.commands.cli_wrappers.Worker",
            "extra_patch": "grzctl.commands.cli_wrappers.validate_module.validate.callback",
            "id_source": "submission",
            "wrong_state": SubmissionStateEnum.ENCRYPTED,  # expected: DECRYPTED
            "intermediate_state": SubmissionStateEnum.VALIDATING,
            "expected_state": SubmissionStateEnum.VALIDATED,
        },
        {
            "cmd": ["encrypt", "--submission-dir", "SUBMISSION_DIR"],
            "worker_patch": "grzctl.commands.cli_wrappers.Worker",
            "extra_patch": "grzctl.commands.cli_wrappers.encrypt_module.encrypt.callback",
            "id_source": "submission",
            "wrong_state": SubmissionStateEnum.DOWNLOADED,  # expected: VALIDATED
            "intermediate_state": SubmissionStateEnum.ENCRYPTING,
            "expected_state": SubmissionStateEnum.ENCRYPTED,
        },
        {
            "cmd": ["archive", "--submission-dir", "SUBMISSION_DIR"],
            "worker_patch": "grzctl.commands.archive.Worker",
            "id_source": "encrypted_submission",
            "wrong_state": SubmissionStateEnum.DOWNLOADED,  # expected: ENCRYPTED
            "intermediate_state": SubmissionStateEnum.ARCHIVING,
            "expected_state": SubmissionStateEnum.ARCHIVED,
        },
    ],
)
def test_db_wrappers_wrong_initial_state(
    command_spec,
    db_engine,
    full_config_path,
    test_metadata,
    tmp_path,
):
    """
    Test that commands still run (with a warning) when the submission is in an unexpected state.
    The DbContext only logs a warning for mismatched states and does not block execution.
    The state transition (intermediate → expected) should still be recorded.
    """
    runner = click.testing.CliRunner()
    cli = build_cli()

    parsed_metadata, metadata_path_fixture = test_metadata
    submission_id = parsed_metadata.submission_id

    submission_dir = tmp_path / "submission"
    submission_dir.mkdir()
    for d in ["metadata", "files", "logs", "encrypted_files"]:
        (submission_dir / d).mkdir()
    shutil.copy(metadata_path_fixture, submission_dir / "metadata" / "metadata.json")

    # Add submission with a wrong/unexpected initial state
    setup_db_state(
        runner,
        cli,
        full_config_path,
        submission_id,
        metadata_path_fixture,
        initial_state=command_spec["wrong_state"],
    )

    args = build_args(
        command_spec,
        submission_dir=submission_dir,
        submission_id=submission_id,
        config_path=full_config_path,
        extra_flags=["--update-db"],
    )

    with mock_command(command_spec, submission_id):
        result = runner.invoke(cli, args)

        # The command should still succeed — wrong state is only a warning, not a hard failure
        assert result.exit_code == 0, (
            f"Command '{command_spec['cmd'][0]}' unexpectedly failed with wrong initial state. Output: {result.output}"
        )

        history = get_state_history(db_engine, submission_id)
        assert history[-1] == command_spec["expected_state"], f"Unexpected final state: {history}"
        assert history[-2] == command_spec["intermediate_state"], f"Unexpected intermediate state: {history}"
        # The wrong state should be present earlier in history
        assert command_spec["wrong_state"] in history, f"Wrong state not found in history: {history}"


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
                "--config",
                str(full_config_path),
                "pruefbericht",
                "submit",
                "--submission-id",
                submission_id,
                "--pruefbericht-file",
                str(pb_path),
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
            "--config",
            str(full_config_path),
            "clean",
            "--submission-id",
            submission_id,
            "--yes-i-really-mean-it",
            "--update-db",
        ]

        result = runner.invoke(cli, args)

        assert result.exit_code != 0

        history = get_state_history(db_engine, submission_id)
        assert history[-1] == SubmissionStateEnum.ERROR
        assert history[-2] == SubmissionStateEnum.CLEANING
        assert history[-3] == SubmissionStateEnum.QCED


@pytest.mark.parametrize(
    ("valid_metadata", "expected_basic_qc_passed"),
    [(True, True), (False, None)],
)
def test_validation_basic_qc_passed_update(
    valid_metadata, expected_basic_qc_passed, db_engine, full_config_path, test_metadata, tmp_path
):
    runner = click.testing.CliRunner()
    cli = build_cli()

    parsed_metadata, metadata_path_fixture = test_metadata
    submission_id = parsed_metadata.submission_id

    submission_dir = tmp_path / "submission"
    submission_dir.mkdir()
    for d in ["metadata", "files", "logs", "encrypted_files"]:
        (submission_dir / d).mkdir()
    shutil.copy(metadata_path_fixture, submission_dir / "metadata" / "metadata.json")

    setup_db_state(
        runner,
        cli,
        full_config_path,
        submission_id,
        metadata_path_fixture,
        initial_state=SubmissionStateEnum.DECRYPTED,
    )

    # Mock the validate.callback
    with (
        patch("grzctl.commands.cli_wrappers.validate_module.validate.callback") as mock_validate_callback,
    ):
        # Fail validation on purpose for negative case
        if not valid_metadata:
            mock_validate_callback.side_effect = Exception("validation failed")

        validate_args = [
            "--config",
            str(full_config_path),
            "validate",
            "--submission-dir",
            str(submission_dir),
            "--update-db",
        ]

        # Test to see if basic_qc_passed is updated to true on successful validation
        runner.invoke(cli, validate_args)

    # Create Submission from Db and check
    with Session(db_engine) as db_session:
        submission_from_db = db_session.exec(select(Submission).where(Submission.id == submission_id)).first()
    assert submission_from_db.basic_qc_passed is expected_basic_qc_passed
