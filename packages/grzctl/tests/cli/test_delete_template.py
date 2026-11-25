import datetime

import click.testing
import grzctl.cli
import pytest
import yaml
from grzctl.commands.delete.template import DELETION_REQUEST_TEMPLATE
from grzctl.models.deletion_request import DeletionRequest, deletion_request_problems

TAN_G_A = "3f9a" + "0" * 56 + "0c21"
TAN_G_B = "8b1e" + "1" * 56 + "00d4"


def _filled_in(**overrides) -> dict:
    """The template filled in as from a request form, with *overrides* applied."""
    request = yaml.safe_load(DELETION_REQUEST_TEMPLATE)
    request.update(
        instructed_at=datetime.date(2026, 11, 21),
        submitter_id="260914050",
        requester_name="Dr. Erika Mustermann",
        requester_email="erika.mustermann@example.org",
        requested_at=datetime.date(2026, 11, 20),
        reason="Widerruf der Einwilligung zur Teilnahme am Modellvorhaben\n",
        data_node_id="GRZM00006",
        tan_g=[TAN_G_A, TAN_G_B],
    )
    request.update(overrides)
    return request


def test_command_prints_the_template():
    result = click.testing.CliRunner().invoke(grzctl.cli.build_cli(), ["delete", "template"])

    assert result.exit_code == 0, result.output
    assert result.output == DELETION_REQUEST_TEMPLATE


def test_unfilled_template_reports_each_field_to_fill_in_once():
    problems = deletion_request_problems(yaml.safe_load(DELETION_REQUEST_TEMPLATE))

    assert [problem.split(":")[0] for problem in problems] == [
        "instructed_at",
        "submitter_id",
        "requester_name",
        "requester_email",
        "requested_at",
        "reason",
        "data_node_id",
        "tan_g",
    ]


def test_filled_in_template_has_no_problems():
    request = _filled_in()

    assert deletion_request_problems(request) == []
    assert DeletionRequest.model_validate(request).tan_g == [TAN_G_A, TAN_G_B]


@pytest.mark.parametrize(
    ("overrides", "location"),
    [
        pytest.param({"tan_g": [TAN_G_A[:-1]]}, "tan_g.0", id="short-tan-g"),
        pytest.param({"tan_g": [TAN_G_A, TAN_G_A.upper()]}, "tan_g", id="repeated-tan-g"),
        pytest.param({"tan_g": []}, "tan_g", id="no-tan-g"),
        pytest.param({"submitter_id": "26091405"}, "submitter_id", id="short-le-id"),
        pytest.param({"data_node_id": "KDKM00001"}, "data_node_id", id="kdk-id"),
        pytest.param({"requester_email": "erika"}, "requester_email", id="bad-email"),
        pytest.param({"requested_at": "20.11.2026"}, "requested_at", id="german-date"),
        pytest.param({"unknown": "x"}, "unknown", id="unknown-field"),
    ],
)
def test_reports_invalid_values(overrides: dict, location: str):
    problems = deletion_request_problems(_filled_in(**overrides))

    assert [problem.split(":")[0] for problem in problems] == [location]
