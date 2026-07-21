"""
Tests for the ``grzctl db case`` subcommand group.
"""

import json
from pathlib import Path

import click.testing
import grzctl.cli
from grz_pydantic_models.submission.metadata import GrzSubmissionMetadata


def _invoke(cli, config_path: Path, *args: str) -> click.testing.Result:
    """Invoke ``grzctl db ... `` with the given config file and arguments."""
    runner = click.testing.CliRunner()
    return runner.invoke(cli, ["db", "--config-file", str(config_path), *args])


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


def _add_and_populate(cli, config_path: Path, test_metadata_path: Path) -> str:
    """Add + populate the example submission, auto-creating its case. Returns the submission id."""
    submission_id = GrzSubmissionMetadata.model_validate_json(test_metadata_path.read_text()).submission_id

    result_add = _invoke(cli, config_path, "submission", "add", submission_id)
    assert result_add.exit_code == 0, result_add.stderr

    result_populate = _invoke(
        cli, config_path, "submission", "populate", submission_id, str(test_metadata_path), "--no-confirm"
    )
    assert result_populate.exit_code == 0, result_populate.stderr

    return submission_id


def test_case_create_duplicate_pair_allowed(blank_database_config_path: Path):
    """Two cases with the same (submitter_id, local_case_id) both succeed; that pair is not unique."""
    cli = grzctl.cli.build_cli()

    result_create = _create_case(cli, blank_database_config_path, "123456789", "case-A")
    assert result_create.exit_code == 0, result_create.stderr

    # the created case id should be reported (message goes to stderr)
    cases = _list_cases(cli, blank_database_config_path)
    assert len(cases) == 1
    assert str(cases[0]["id"]) in result_create.stderr

    # re-creating the same (submitter_id, local_case_id) pair now also succeeds
    result_duplicate = _create_case(cli, blank_database_config_path, "123456789", "case-A")
    assert result_duplicate.exit_code == 0, result_duplicate.stderr

    # ... resulting in two independent cases
    assert len(_list_cases(cli, blank_database_config_path)) == 2


def test_case_create_duplicate_psn_rejected(blank_database_config_path: Path):
    """Uniqueness is enforced on --psn: reusing a psn makes the second create fail."""
    cli = grzctl.cli.build_cli()

    result_first = _create_case(cli, blank_database_config_path, "123456789", "case-A", psn="RKI-1")
    assert result_first.exit_code == 0, result_first.stderr

    # a different (submitter_id, local_case_id) but the same psn must fail
    result_second = _create_case(cli, blank_database_config_path, "987654321", "case-B", psn="RKI-1")
    assert result_second.exit_code != 0

    # only the first case was created
    assert len(_list_cases(cli, blank_database_config_path)) == 1


def test_case_list_and_show_json(blank_database_config_path: Path):
    """``case list --json`` and ``case show --json`` expose the expected fields."""
    cli = grzctl.cli.build_cli()

    result_create = _create_case(cli, blank_database_config_path, "123456789", "case-B")
    assert result_create.exit_code == 0, result_create.stderr

    cases = _list_cases(cli, blank_database_config_path)
    assert len(cases) == 1
    case = cases[0]
    for key in ("id", "psn", "submitter_id", "local_case_id", "submission_count"):
        assert key in case
    assert case["submitter_id"] == "123456789"
    assert case["local_case_id"] == "case-B"
    assert case["psn"] is None
    assert case["submission_count"] == 0

    case_id = case["id"]
    result_show = _invoke(cli, blank_database_config_path, "case", "show", str(case_id), "--json")
    assert result_show.exit_code == 0, result_show.stderr
    shown = json.loads(result_show.stdout)
    for key in ("id", "psn", "submitter_id", "local_case_id", "submissions"):
        assert key in shown
    assert shown["id"] == case_id
    assert shown["submitter_id"] == "123456789"
    assert shown["local_case_id"] == "case-B"
    assert shown["psn"] is None
    assert shown["submissions"] == []


def test_case_modify_psn(blank_database_config_path: Path):
    """``case modify <id> psn`` updates the PSN; duplicates and unknown keys are rejected."""
    cli = grzctl.cli.build_cli()

    result_create = _create_case(cli, blank_database_config_path, "123456789", "case-C")
    assert result_create.exit_code == 0, result_create.stderr
    case_id = _list_cases(cli, blank_database_config_path)[0]["id"]

    result_modify = _invoke(cli, blank_database_config_path, "case", "modify", str(case_id), "psn", "RKI-1")
    assert result_modify.exit_code == 0, result_modify.stderr

    result_show = _invoke(cli, blank_database_config_path, "case", "show", str(case_id), "--json")
    assert result_show.exit_code == 0, result_show.stderr
    assert json.loads(result_show.stdout)["psn"] == "RKI-1"

    # submitter_id is now mutable (only ``id`` is immutable)
    result_modify_submitter = _invoke(
        cli, blank_database_config_path, "case", "modify", str(case_id), "submitter_id", "999999999"
    )
    assert result_modify_submitter.exit_code == 0, result_modify_submitter.stderr
    result_show_submitter = _invoke(cli, blank_database_config_path, "case", "show", str(case_id), "--json")
    assert json.loads(result_show_submitter.stdout)["submitter_id"] == "999999999"

    # modifying psn to a value already used by another case is rejected
    result_other = _create_case(cli, blank_database_config_path, "123456789", "case-D", psn="RKI-2")
    assert result_other.exit_code == 0, result_other.stderr
    result_modify_dup = _invoke(cli, blank_database_config_path, "case", "modify", str(case_id), "psn", "RKI-2")
    assert result_modify_dup.exit_code != 0

    # an unknown modify KEY is rejected by the click Choice
    result_modify_bad_key = _invoke(
        cli, blank_database_config_path, "case", "modify", str(case_id), "not_a_real_key", "value"
    )
    assert result_modify_bad_key.exit_code != 0


def test_case_show_unknown_exits_nonzero(blank_database_config_path: Path):
    """Showing a nonexistent case id exits non-zero."""
    cli = grzctl.cli.build_cli()
    result = _invoke(cli, blank_database_config_path, "case", "show", "999999")
    assert result.exit_code != 0


def test_case_delete_empty(blank_database_config_path: Path):
    """An empty case can be deleted, after which it can no longer be shown."""
    cli = grzctl.cli.build_cli()

    result_create = _create_case(cli, blank_database_config_path, "123456789", "case-D")
    assert result_create.exit_code == 0, result_create.stderr
    case_id = _list_cases(cli, blank_database_config_path)[0]["id"]

    result_delete = _invoke(cli, blank_database_config_path, "case", "delete", str(case_id))
    assert result_delete.exit_code == 0, result_delete.stderr

    result_show = _invoke(cli, blank_database_config_path, "case", "show", str(case_id))
    assert result_show.exit_code != 0


def test_case_delete_refuses_when_linked(blank_database_config_path: Path, test_metadata_path: Path):
    """Deleting a case that still has linked submissions is refused."""
    cli = grzctl.cli.build_cli()

    _add_and_populate(cli, blank_database_config_path, test_metadata_path)

    # populate auto-creates exactly one case for the submission
    cases = _list_cases(cli, blank_database_config_path)
    assert len(cases) == 1
    linked_case = cases[0]
    assert linked_case["submission_count"] == 1

    result_delete = _invoke(cli, blank_database_config_path, "case", "delete", str(linked_case["id"]))
    assert result_delete.exit_code != 0


def test_case_relink(blank_database_config_path: Path, test_metadata_path: Path):
    """A submission can be relinked to a different case."""
    cli = grzctl.cli.build_cli()

    submission_id = _add_and_populate(cli, blank_database_config_path, test_metadata_path)

    cases = _list_cases(cli, blank_database_config_path)
    assert len(cases) == 1
    original_case_id = cases[0]["id"]

    # create a second, empty case to relink the (test-type) submission into
    result_create = _create_case(cli, blank_database_config_path, "123456789", "relink-target")
    assert result_create.exit_code == 0, result_create.stderr
    new_case_id = next(c["id"] for c in _list_cases(cli, blank_database_config_path) if c["id"] != original_case_id)

    result_relink = _invoke(cli, blank_database_config_path, "case", "relink", submission_id, str(new_case_id))
    assert result_relink.exit_code == 0, result_relink.stderr

    # the submission now shows up under the new case
    result_show = _invoke(cli, blank_database_config_path, "case", "show", str(new_case_id), "--json")
    assert result_show.exit_code == 0, result_show.stderr
    shown = json.loads(result_show.stdout)
    assert submission_id in {s["id"] for s in shown["submissions"]}
