"""
Tests for the ``grzctl db case`` subcommand group.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

import click.testing
import grzctl.cli
import pytest
from grz_pydantic_models.submission.metadata import GrzSubmissionMetadata


def _invoke(cli, config_path: Path, *args: str, input: str | None = None) -> click.testing.Result:
    """Invoke ``grzctl db ... `` with the given config file and arguments."""
    runner = click.testing.CliRunner()
    return runner.invoke(cli, ["--config", str(config_path), "db", *args], input=input)


def _create_case(
    cli, config_path: Path, submitter_id: str, local_case_id: str, psn: str | None = None
) -> click.testing.Result:
    args = ["case", "create", submitter_id, local_case_id]
    if psn is not None:
        args += ["--psn", psn]
    return _invoke(cli, config_path, *args)


def _list_cases(cli, config_path: Path) -> list[dict]:
    result = _invoke(cli, config_path, "case", "list", "--json")
    assert result.exit_code == 0, result.stderr
    return json.loads(result.stdout)


def _as_initial_submission(test_metadata_path: Path, tmp_path: Path) -> Path:
    """Write a copy of *test_metadata_path* typed ``initial``, so populating it creates a case."""
    raw = json.loads(test_metadata_path.read_text())
    raw["submission"]["submissionType"] = "initial"
    metadata_path = tmp_path / "initial.metadata.json"
    metadata_path.write_text(json.dumps(raw))
    return metadata_path


def _add_and_populate(cli, config_path: Path, test_metadata_path: Path) -> str:
    """Add + populate the example submission (test-type, so no case is created). Returns the submission id."""
    submission_id = GrzSubmissionMetadata.model_validate_json(test_metadata_path.read_text()).submission_id

    result_add = _invoke(cli, config_path, "submission", "add", submission_id)
    assert result_add.exit_code == 0, result_add.stderr

    result_populate = _invoke(
        cli, config_path, "submission", "populate", submission_id, str(test_metadata_path), "--no-confirm"
    )
    assert result_populate.exit_code == 0, result_populate.stderr

    return submission_id


def test_case_create_duplicate_pair_rejected(migrated_database_config_path: Path):
    """A submitter's local case ID identifies one patient, so it identifies one case."""
    cli = grzctl.cli.build_cli()

    result_create = _create_case(cli, migrated_database_config_path, "123456789", "case-A")
    assert result_create.exit_code == 0, result_create.stderr

    # the created case id should be reported (message goes to stderr)
    cases = _list_cases(cli, migrated_database_config_path)
    assert len(cases) == 1
    assert str(cases[0]["id"]) in result_create.stderr

    result_duplicate = _create_case(cli, migrated_database_config_path, "123456789", "case-A")

    assert result_duplicate.exit_code != 0
    assert "already exists" in result_duplicate.stderr
    assert len(_list_cases(cli, migrated_database_config_path)) == 1


def test_case_create_duplicate_psn_rejected(migrated_database_config_path: Path):
    """Uniqueness is enforced on --psn: reusing a psn makes the second create fail."""
    cli = grzctl.cli.build_cli()

    result_first = _create_case(cli, migrated_database_config_path, "123456789", "case-A", psn="RKI-1")
    assert result_first.exit_code == 0, result_first.stderr

    # a different (submitter_id, local_case_id) but the same psn must fail
    result_second = _create_case(cli, migrated_database_config_path, "987654321", "case-B", psn="RKI-1")
    assert result_second.exit_code != 0

    # only the first case was created
    assert len(_list_cases(cli, migrated_database_config_path)) == 1


def test_case_create_invalid_submitter_id_rejected(migrated_database_config_path: Path):
    """A malformed submitter ID fails ``Case``'s pydantic validation, not an uncaught crash.

    ``submitter_id`` is pattern-constrained (``^[0-9]{9}$``) and ``Case`` validates on assignment,
    so this is a ``ValidationError`` (a ``ValueError`` subclass), distinct from the db-level
    ``DuplicateCaseError``/``DuplicatePsnError`` this command already catches.
    """
    cli = grzctl.cli.build_cli()

    result = _create_case(cli, migrated_database_config_path, "not-a-submitter-id", "case-X")

    assert result.exit_code != 0
    # a caught, reported error, not an exception that escaped to the runner uncaught
    assert "Error:" in result.stderr
    assert "Traceback" not in result.stderr
    assert len(_list_cases(cli, migrated_database_config_path)) == 0


def test_case_list_and_show_json(migrated_database_config_path: Path):
    """``case list --json`` and ``case show --json`` expose the expected fields."""
    cli = grzctl.cli.build_cli()

    result_create = _create_case(cli, migrated_database_config_path, "123456789", "case-B")
    assert result_create.exit_code == 0, result_create.stderr

    cases = _list_cases(cli, migrated_database_config_path)
    assert len(cases) == 1
    case = cases[0]
    for key in ("id", "psn", "submitter_id", "local_case_id", "submission_count"):
        assert key in case
    assert case["submitter_id"] == "123456789"
    assert case["local_case_id"] == "case-B"
    assert case["psn"] is None
    assert case["submission_count"] == 0

    case_id = case["id"]
    result_show = _invoke(cli, migrated_database_config_path, "case", "show", str(case_id), "--json")
    assert result_show.exit_code == 0, result_show.stderr
    shown = json.loads(result_show.stdout)
    for key in ("id", "psn", "submitter_id", "local_case_id", "submissions"):
        assert key in shown
    assert shown["id"] == case_id
    assert shown["submitter_id"] == "123456789"
    assert shown["local_case_id"] == "case-B"
    assert shown["psn"] is None
    assert shown["submissions"] == []


def test_case_modify_psn(migrated_database_config_path: Path):
    """``case modify <id> psn`` updates the PSN; duplicates and unknown keys are rejected."""
    cli = grzctl.cli.build_cli()

    result_create = _create_case(cli, migrated_database_config_path, "123456789", "case-C")
    assert result_create.exit_code == 0, result_create.stderr
    case_id = _list_cases(cli, migrated_database_config_path)[0]["id"]

    result_modify = _invoke(cli, migrated_database_config_path, "case", "modify", str(case_id), "psn", "RKI-1")
    assert result_modify.exit_code == 0, result_modify.stderr

    result_show = _invoke(cli, migrated_database_config_path, "case", "show", str(case_id), "--json")
    assert result_show.exit_code == 0, result_show.stderr
    assert json.loads(result_show.stdout)["psn"] == "RKI-1"

    # submitter_id is mutable (only ``id`` is immutable)
    result_modify_submitter = _invoke(
        cli, migrated_database_config_path, "case", "modify", str(case_id), "submitter_id", "999999999"
    )
    assert result_modify_submitter.exit_code == 0, result_modify_submitter.stderr
    result_show_submitter = _invoke(cli, migrated_database_config_path, "case", "show", str(case_id), "--json")
    assert json.loads(result_show_submitter.stdout)["submitter_id"] == "999999999"

    # modifying psn to a value already used by another case is rejected
    result_other = _create_case(cli, migrated_database_config_path, "123456789", "case-D", psn="RKI-2")
    assert result_other.exit_code == 0, result_other.stderr
    result_modify_dup = _invoke(cli, migrated_database_config_path, "case", "modify", str(case_id), "psn", "RKI-2")
    assert result_modify_dup.exit_code != 0

    # an unknown modify KEY is rejected by the click Choice
    result_modify_bad_key = _invoke(
        cli, migrated_database_config_path, "case", "modify", str(case_id), "not_a_real_key", "value"
    )
    assert result_modify_bad_key.exit_code != 0


def test_case_modify_duplicate_pair_rejected(migrated_database_config_path: Path):
    """Modifying ``local_case_id`` into another case's ``(submitter_id, local_case_id)`` pair is
    rejected cleanly, not an uncaught crash.

    The pair is enforced by the partial unique index ``ux_cases_submitter_local_case``, so this
    raises ``DuplicateCaseError``, distinct from the ``CaseNotFoundError``/``DuplicatePsnError``/
    ``ValueError`` this command already catches.
    """
    cli = grzctl.cli.build_cli()

    result_first = _create_case(cli, migrated_database_config_path, "123456789", "case-A")
    assert result_first.exit_code == 0, result_first.stderr
    result_second = _create_case(cli, migrated_database_config_path, "123456789", "case-B")
    assert result_second.exit_code == 0, result_second.stderr
    second_case_id = next(
        c["id"] for c in _list_cases(cli, migrated_database_config_path) if c["local_case_id"] == "case-B"
    )

    result_modify = _invoke(
        cli, migrated_database_config_path, "case", "modify", str(second_case_id), "local_case_id", "case-A"
    )

    assert result_modify.exit_code != 0
    # a caught, reported error, not an exception that escaped to the runner uncaught
    assert "Error:" in result_modify.stderr
    assert "Traceback" not in result_modify.stderr

    # the change must not have gone through
    unchanged = next(
        c["local_case_id"] for c in _list_cases(cli, migrated_database_config_path) if c["id"] == second_case_id
    )
    assert unchanged == "case-B"


def test_case_show_unknown_exits_nonzero(migrated_database_config_path: Path):
    """Showing a nonexistent case id exits non-zero."""
    cli = grzctl.cli.build_cli()
    result = _invoke(cli, migrated_database_config_path, "case", "show", "999999")
    assert result.exit_code != 0


def test_case_delete_empty(migrated_database_config_path: Path):
    """An empty case can be deleted, after which it can no longer be shown."""
    cli = grzctl.cli.build_cli()

    result_create = _create_case(cli, migrated_database_config_path, "123456789", "case-D")
    assert result_create.exit_code == 0, result_create.stderr
    case_id = _list_cases(cli, migrated_database_config_path)[0]["id"]

    result_delete = _invoke(cli, migrated_database_config_path, "case", "delete", str(case_id))
    assert result_delete.exit_code == 0, result_delete.stderr

    result_show = _invoke(cli, migrated_database_config_path, "case", "show", str(case_id))
    assert result_show.exit_code != 0


def test_populate_creates_no_case_for_test_submission(migrated_database_config_path: Path, test_metadata_path: Path):
    """Test submissions are never case-tracked: populating one leaves the case table empty."""
    cli = grzctl.cli.build_cli()

    _add_and_populate(cli, migrated_database_config_path, test_metadata_path)

    assert _list_cases(cli, migrated_database_config_path) == []


def test_case_delete_refuses_when_linked(migrated_database_config_path: Path):
    """Deleting a case that still has linked submissions is refused."""
    cli = grzctl.cli.build_cli()

    submission_id = "123456789_2025-01-01_0000000a"
    result_add = _invoke(cli, migrated_database_config_path, "submission", "add", submission_id)
    assert result_add.exit_code == 0, result_add.stderr
    # a case only takes a submission whose type is known and is not a test submission
    result_type = _invoke(
        cli, migrated_database_config_path, "submission", "modify", submission_id, "submission_type", "initial"
    )
    assert result_type.exit_code == 0, result_type.stderr

    # a submission can be linked to a case created by hand, not only one populate resolved
    result_create = _create_case(cli, migrated_database_config_path, "123456789", "delete-target")
    assert result_create.exit_code == 0, result_create.stderr
    linked_case = _list_cases(cli, migrated_database_config_path)[0]
    result_relink = _invoke(cli, migrated_database_config_path, "case", "relink", submission_id, str(linked_case["id"]))
    assert result_relink.exit_code == 0, result_relink.stderr

    result_delete = _invoke(cli, migrated_database_config_path, "case", "delete", str(linked_case["id"]))
    assert result_delete.exit_code != 0


def test_case_relink(migrated_database_config_path: Path):
    """A submission can be relinked to a different case."""
    cli = grzctl.cli.build_cli()

    submission_id = "123456789_2025-01-01_0000000a"
    result_add = _invoke(cli, migrated_database_config_path, "submission", "add", submission_id)
    assert result_add.exit_code == 0, result_add.stderr
    # a case only takes a submission whose type is known and is not a test submission
    result_type = _invoke(
        cli, migrated_database_config_path, "submission", "modify", submission_id, "submission_type", "initial"
    )
    assert result_type.exit_code == 0, result_type.stderr

    result_first = _create_case(cli, migrated_database_config_path, "123456789", "first-case")
    assert result_first.exit_code == 0, result_first.stderr
    original_case_id = _list_cases(cli, migrated_database_config_path)[0]["id"]
    result_link = _invoke(cli, migrated_database_config_path, "case", "relink", submission_id, str(original_case_id))
    assert result_link.exit_code == 0, result_link.stderr

    # create a second, empty case to relink the submission into
    result_create = _create_case(cli, migrated_database_config_path, "123456789", "relink-target")
    assert result_create.exit_code == 0, result_create.stderr
    new_case_id = next(c["id"] for c in _list_cases(cli, migrated_database_config_path) if c["id"] != original_case_id)

    result_relink = _invoke(cli, migrated_database_config_path, "case", "relink", submission_id, str(new_case_id))
    assert result_relink.exit_code == 0, result_relink.stderr

    # the submission now shows up under the new case
    result_show = _invoke(cli, migrated_database_config_path, "case", "show", str(new_case_id), "--json")
    assert result_show.exit_code == 0, result_show.stderr
    shown = json.loads(result_show.stdout)
    assert submission_id in {s["id"] for s in shown["submissions"]}
    # basic QC state is what tells the case's QC-passed initial submission apart from any
    # duplicate initial submissions of the same case, so `case show` must carry it
    assert all("basic_qc_passed" in s for s in shown["submissions"])


def test_populate_declined_confirmation_creates_no_case(
    migrated_database_config_path: Path, test_metadata_path: Path, tmp_path: Path
):
    """Answering 'no' at the populate confirmation must leave the database unchanged, cases included.

    The metadata is retyped ``initial`` first: the shipped example is a test submission, which is
    never case-tracked, so against it the assertion would hold whatever the answer was.
    """
    cli = grzctl.cli.build_cli()
    metadata_path = _as_initial_submission(test_metadata_path, tmp_path)
    submission_id = GrzSubmissionMetadata.model_validate_json(metadata_path.read_text()).submission_id

    result_add = _invoke(cli, migrated_database_config_path, "submission", "add", submission_id)
    assert result_add.exit_code == 0, result_add.stderr

    result = _invoke(
        cli, migrated_database_config_path, "submission", "populate", submission_id, str(metadata_path), input="n\n"
    )
    assert result.exit_code == 0, result.stderr
    assert "Database populated successfully" not in result.stderr
    assert _list_cases(cli, migrated_database_config_path) == []

    # accepting the same populate does create one, so the decline is what the assertion above rests on
    result_accepted = _invoke(
        cli, migrated_database_config_path, "submission", "populate", submission_id, str(metadata_path), "--no-confirm"
    )
    assert result_accepted.exit_code == 0, result_accepted.stderr
    assert len(_list_cases(cli, migrated_database_config_path)) == 1


@pytest.mark.parametrize("command", [("case", "list"), ("case", "show", "{case_id}")])
def test_case_json_survives_a_coloured_environment(migrated_database_config_path, command):
    """``--json`` output must be JSON for a consumer piping it, whatever the shell exports.

    rich decides on ANSI colour when a Console is built, which for the CLI modules happens at
    import, and it honours ``FORCE_COLOR`` even when stdout is a pipe. This runs out of process
    because the test suite pins ``TTY_COMPATIBLE=0`` at import for every in-process test: an
    in-process run would inherit that pin and never see escapes a real user's shell can
    produce.
    """
    cli = grzctl.cli.build_cli()
    assert _create_case(cli, migrated_database_config_path, "111111111", "case-colour").exit_code == 0
    case_id = _list_cases(cli, migrated_database_config_path)[0]["id"]

    environment = {**os.environ, "FORCE_COLOR": "1"}
    environment.pop("TTY_COMPATIBLE", None)
    completed = subprocess.run(
        [
            str(Path(sys.executable).parent / "grzctl"),
            "--config",
            str(migrated_database_config_path),
            "db",
            *(part.format(case_id=case_id) for part in command),
            "--json",
        ],
        capture_output=True,
        text=True,
        env=environment,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    assert "\x1b" not in completed.stdout, "ANSI escapes reached a consumer parsing this"
    json.loads(completed.stdout)  # raises if anything but JSON was written
