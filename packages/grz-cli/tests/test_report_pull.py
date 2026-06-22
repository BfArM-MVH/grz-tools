"""Tests for `grz-cli report pull`."""

import gzip
import io
import json
from datetime import UTC, datetime
from pathlib import Path

import click.testing
from grz_cli.cli import build_cli


def _status_payload_json() -> str:
    payload = {
        "generated_at": datetime(2026, 6, 2, 12, 0, tzinfo=UTC).isoformat(),
        "grz_id": "GRZX00000",
        "le_id": "123456789",
        "submissions": [
            {
                "submission_id": "123456789_2025-09-15_d0f805c5",
                "submission_date": "2025-09-15",
                "latest_state": "Reported",
                "latest_state_at": datetime(2026, 6, 2, 11, 59, tzinfo=UTC).isoformat(),
                "failure_reason": None,
                "basic_qc_passed": True,
                "detailed_qc_passed": None,
                "has_detailed_qc_report": False,
                "reported_date": "2026-06-01",
            },
            {
                "submission_id": "123456789_2025-10-15_aaaaaaaa",
                "submission_date": "2025-10-15",
                "latest_state": "Error",
                "latest_state_at": datetime(2026, 6, 2, 10, 0, tzinfo=UTC).isoformat(),
                "failure_reason": "validation_error",
                "basic_qc_passed": False,
                "detailed_qc_passed": False,
                "has_detailed_qc_report": True,
                "reported_date": None,
            },
        ],
    }
    return json.dumps(payload)


class _FakeS3Client:
    def __init__(self, body: bytes):
        self._body = body

    def get_object(self, Bucket: str, Key: str):  # noqa: N803
        assert Key == ".status.json.gz"
        return {
            "Body": io.BytesIO(self._body),
            "ContentEncoding": "gzip",
        }


def test_report_pull_prints_decoded_json(s3_config_path: Path, monkeypatch):
    payload = _status_payload_json().encode("utf-8")
    client = _FakeS3Client(gzip.compress(payload))
    monkeypatch.setattr("grz_cli.commands.report.init_s3_client", lambda _opts: client)

    runner = click.testing.CliRunner()
    result = runner.invoke(
        build_cli(),
        [
            "--config-file",
            str(s3_config_path),
            "report",
            "pull",
            "--json",
        ],
        catch_exceptions=False,
    )

    assert result.exit_code == 0, result.output
    decoded = json.loads(result.output)
    assert decoded["le_id"] == "123456789"
    assert decoded["submissions"][0]["submission_id"] == "123456789_2025-09-15_d0f805c5"


def test_report_pull_can_write_output_file(s3_config_path: Path, tmp_path: Path, monkeypatch):
    payload = _status_payload_json().encode("utf-8")
    client = _FakeS3Client(gzip.compress(payload))
    monkeypatch.setattr("grz_cli.commands.report.init_s3_client", lambda _opts: client)

    output_path = tmp_path / "status.json"
    runner = click.testing.CliRunner()
    result = runner.invoke(
        build_cli(),
        [
            "--config-file",
            str(s3_config_path),
            "report",
            "pull",
            "--json",
            "--output",
            str(output_path),
        ],
        catch_exceptions=False,
    )

    assert result.exit_code == 0, result.output
    assert output_path.exists()
    decoded = json.loads(output_path.read_text(encoding="utf-8"))
    assert decoded["grz_id"] == "GRZX00000"


def test_report_pull_defaults_to_summary_output(s3_config_path: Path, monkeypatch):
    payload = _status_payload_json().encode("utf-8")
    client = _FakeS3Client(gzip.compress(payload))
    monkeypatch.setattr("grz_cli.commands.report.init_s3_client", lambda _opts: client)

    runner = click.testing.CliRunner()
    result = runner.invoke(
        build_cli(),
        [
            "--config-file",
            str(s3_config_path),
            "report",
            "pull",
        ],
        catch_exceptions=False,
    )

    assert result.exit_code == 0, result.output
    assert "Submission Status Summary" in result.output
    assert "Total Submissions" in result.output and "2" in result.output
    assert "State: Reported" in result.output and "1" in result.output
    assert "State: Error" in result.output and "1" in result.output
    assert "Failure: validation_error" in result.output and "1" in result.output


def test_report_pull_fails_on_mismatching_configured_le_id(s3_config_path: Path, monkeypatch):
    payload = _status_payload_json().replace("123456789", "987654321", 1).encode("utf-8")
    client = _FakeS3Client(gzip.compress(payload))
    monkeypatch.setattr("grz_cli.commands.report.init_s3_client", lambda _opts: client)

    runner = click.testing.CliRunner()
    result = runner.invoke(
        build_cli(),
        [
            "--config-file",
            str(s3_config_path),
            "report",
            "pull",
            "--json",
        ],
    )

    assert result.exit_code != 0
    assert "Configured identifiers.le does not match report LE ID" in result.output
