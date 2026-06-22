"""Commands for LE-facing reports."""

import gzip
from pathlib import Path
from typing import Any

import click
import grz_common.cli as grzcli
import rich.console
import rich.table
from grz_common.transfer import init_s3_client
from grz_pydantic_models.status import SubmissionStatusReport

from ..models.config import ReportPullConfig

_STATUS_KEY = ".status.json.gz"


def _build_summary_table(status_report: SubmissionStatusReport) -> rich.table.Table:
    table = rich.table.Table(title="Submission Status Summary")
    table.add_column("Metric", style="cyan")
    table.add_column("Value", justify="right")

    submissions = status_report.submissions
    latest_state_counts: dict[str, int] = {}
    failure_reason_counts: dict[str, int] = {}

    for submission in submissions:
        latest_state = submission.latest_state or "unknown"
        latest_state_counts[latest_state] = latest_state_counts.get(latest_state, 0) + 1

        if submission.failure_reason:
            failure_reason_counts[submission.failure_reason] = (
                failure_reason_counts.get(submission.failure_reason, 0) + 1
            )

    table.add_row("GRZ ID", status_report.grz_id)
    table.add_row("LE ID", status_report.le_id)
    table.add_row("Generated At", status_report.generated_at.isoformat())
    table.add_row("Total Submissions", str(len(submissions)))
    table.add_row(
        "Basic QC Passed",
        str(sum(1 for s in submissions if s.basic_qc_passed is True)),
    )
    table.add_row(
        "Detailed QC Passed",
        str(sum(1 for s in submissions if s.detailed_qc_passed is True)),
    )
    table.add_row(
        "Detailed QC Failed",
        str(sum(1 for s in submissions if s.detailed_qc_passed is False)),
    )
    table.add_row(
        "Reported Submissions",
        str(sum(1 for s in submissions if s.reported_date is not None)),
    )

    if latest_state_counts:
        for state, count in sorted(latest_state_counts.items()):
            table.add_row(f"State: {state}", str(count))

    if failure_reason_counts:
        for reason, count in sorted(failure_reason_counts.items()):
            table.add_row(f"Failure: {reason}", str(count))

    return table


@click.group()
def report():
    """Read a status report from the configured bucket.

    These commands use the S3 settings from the config.
    """


@report.command("pull")
@grzcli.configuration
@grzcli.output_json
@click.option(
    "--output",
    "output_path",
    type=click.Path(dir_okay=False, writable=True, resolve_path=True, path_type=Path),
    default=None,
    help="Optional path to write the output.",
)
def pull(configuration: dict[str, Any], output_json: bool, output_path: Path | None, **kwargs):
    """Download the status report from S3.

    Use ``--json`` for full machine-readable output; otherwise a compact summary is shown.
    """
    config = ReportPullConfig.model_validate(configuration)
    s3_client = init_s3_client(config.s3)

    response = s3_client.get_object(Bucket=config.s3.bucket, Key=_STATUS_KEY)
    payload_bytes = response["Body"].read()

    content_encoding = response.get("ContentEncoding")
    if content_encoding == "gzip" or _STATUS_KEY.endswith(".gz"):
        payload_bytes = gzip.decompress(payload_bytes)

    decoded_payload = payload_bytes.decode("utf-8")
    status_report = SubmissionStatusReport.model_validate_json(decoded_payload)
    configured_le_id = config.identifiers.le
    if status_report.le_id != configured_le_id:
        raise click.UsageError(
            f"Configured identifiers.le does not match report LE ID ({configured_le_id} != {status_report.le_id}). "
        )

    if output_json:
        rendered_output = status_report.model_dump_json(indent=2)
    else:
        summary_console = rich.console.Console(record=True)
        summary_console.print(_build_summary_table(status_report))
        rendered_output = summary_console.export_text()

    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(rendered_output, encoding="utf-8")
        click.echo(f"Wrote status report to {output_path}")
        return

    click.echo(rendered_output)
