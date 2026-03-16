"""Command for submitting Prüfberichte."""

import datetime
import logging
from typing import Any

import click
import grz_common.cli as grzcli
import requests
from grz_common.workers.submission import Submission
from grz_db.models.submission import SubmissionDb, SubmissionStateEnum
from grz_pydantic_models.pruefbericht.v0 import LibraryType as PruefberichtLibraryType
from grz_pydantic_models.pruefbericht.v0 import Pruefbericht, SubmittedCase
from grz_pydantic_models.submission.metadata.v1 import REDACTED_TAN, GrzSubmissionMetadata, Relation
from pydantic_core import to_jsonable_python

from ..dbcontext import DbContext
from ..models.config import DbConfig, PruefberichtConfig

log = logging.getLogger(__name__)
fail_or_pass = click.option(
    "--fail/--pass", "failed", help="Fail an otherwise valid submission (e.g. failed internal QC)"
)


def _get_new_token(auth_url: str, client_id: str, client_secret: str) -> tuple[str, datetime.datetime]:
    log.info("Refreshing access token...")

    response = requests.post(
        auth_url,
        headers={"Content-Type": "application/x-www-form-urlencoded"},
        data={"grant_type": "client_credentials", "client_id": client_id, "client_secret": client_secret},
        timeout=60,
    )

    if response.status_code != requests.codes.ok:
        log.error("There was a problem refreshing the access token")
        response.raise_for_status()

    response_json = response.json()
    token = response_json["access_token"]
    expires_in = response_json["expires_in"]
    # take off a second to provide at least a minimal safety margin
    expires_at = datetime.datetime.now() + datetime.timedelta(seconds=expires_in - 1)

    log.info("Successfully obtained a new access token.")
    return token, expires_at


def _submit_pruefbericht(base_url: str, token: str, pruefbericht: Pruefbericht):
    log.info("Submitting Prüfbericht...")

    response = requests.post(
        base_url.rstrip("/") + "/upload",
        headers={"Authorization": f"bearer {token}"},
        json=to_jsonable_python(pruefbericht),
        timeout=60,
    )

    if response.status_code != requests.codes.ok:
        log.warning("There was a problem submitting the Prüfbericht.")
        response.raise_for_status()


def _get_most_expensive_library_type(library_types: set[str]) -> PruefberichtLibraryType:
    return PruefberichtLibraryType.most_expensive(library_types)


def get_pruefbericht_library_type(metadata: GrzSubmissionMetadata) -> PruefberichtLibraryType:
    """
    Determine the singular representative library type of a submission to submit with the Prüfbericht.
    This should be library type of the index patient with the highest reimbursement value.
    """
    index_patient = metadata.index_donor
    index_patient_submission_library_types = {str(datum.library_type) for datum in index_patient.lab_data}
    return _get_most_expensive_library_type(index_patient_submission_library_types)


def _generate_pruefbericht_from_metadata(metadata: GrzSubmissionMetadata, failed: bool) -> Pruefbericht:
    return Pruefbericht(
        SubmittedCase=SubmittedCase(
            submissionDate=metadata.submission.submission_date,
            submissionType=metadata.submission.submission_type,
            tan=metadata.submission.tan_g,
            submitterId=metadata.submission.submitter_id,
            dataNodeId=metadata.submission.genomic_data_center_id,
            diseaseType=metadata.submission.disease_type,
            dataCategory="genomic",
            libraryType=get_pruefbericht_library_type(metadata),
            coverageType=metadata.submission.coverage_type,
            dataQualityCheckPassed=not failed,
        )
    )


def _generate_pruefbericht_from_database(
    submission_id: str, configuration: dict[str, Any], failed: bool
) -> Pruefbericht:
    """Generate Prüfbericht by fetching submission data from the database."""
    config = DbConfig.model_validate(configuration)

    if not config.db.database_url:
        raise ValueError("database_url must be provided in configuration to fetch from database")

    db = SubmissionDb(db_url=str(config.db.database_url), author=None, debug=False)
    submission = db.get_submission(submission_id)

    if submission is None:
        raise ValueError(f"Submission with ID '{submission_id}' not found in database")

    # Check if submission has the required fields populated
    required_fields = [
        "submission_date",
        "submission_type",
        "tan_g",
        "submitter_id",
        "data_node_id",
        "disease_type",
        "coverage_type",
    ]

    missing_fields = [field for field in required_fields if getattr(submission, field) is None]
    if missing_fields:
        raise ValueError(f"Submission {submission_id} is missing required fields: {', '.join(missing_fields)}")

    # Get donors to determine library types
    donors = db.get_donors(submission_id)
    if not donors:
        raise ValueError(f"No donors found for submission {submission_id}")

    # Find index donor
    index_donor = next((d for d in donors if d.relation == Relation.index_), None)
    if index_donor is None:
        raise ValueError(f"No index donor found for submission {submission_id}")

    # Convert database library_types to strings and determine most expensive type
    index_donor_library_types = {str(lt.value if hasattr(lt, "value") else lt) for lt in index_donor.library_types}

    library_type = _get_most_expensive_library_type(index_donor_library_types)

    # Generate the Prüfbericht
    return Pruefbericht(
        SubmittedCase=SubmittedCase(
            submissionDate=submission.submission_date,
            submissionType=submission.submission_type,
            tan=submission.tan_g,
            submitterId=submission.submitter_id,
            dataNodeId=submission.data_node_id,
            diseaseType=submission.disease_type,
            dataCategory="genomic",
            libraryType=library_type,
            coverageType=submission.coverage_type,
            dataQualityCheckPassed=not failed,
        )
    )


@click.group()
def pruefbericht():
    """Generate and submit Prüfberichte."""


@pruefbericht.group()
def generate():
    """Generate a Prüfbericht JSON from submission metadata."""


@generate.command("from-submission-dir")
@click.argument(
    "submission_dir",
    metavar="PATH",
    type=grzcli.DIR_R_E,
    required=True,
)
@fail_or_pass
def from_submission_dir(submission_dir, failed):
    """Generate Prüfbericht from submission directory.

    This is equivalent to `from-metadata ${submission_dir}/metadata/metadata.json`.
    """
    submission = Submission(metadata_dir=f"{submission_dir}/metadata", files_dir=f"{submission_dir}/files")
    metadata = submission.metadata.content
    pruefbericht = _generate_pruefbericht_from_metadata(metadata, failed)
    click.echo(pruefbericht.model_dump_json(indent=None, by_alias=True))


@generate.command("from-metadata")
@click.argument("metadata_file", type=click.Path(exists=True))
@fail_or_pass
def from_metadata(metadata_file, failed):
    """Generate Prüfbericht from metadata.json"""
    with open(metadata_file) as f:
        metadata = GrzSubmissionMetadata.model_validate_json(f.read())
    pruefbericht = _generate_pruefbericht_from_metadata(metadata, failed)
    click.echo(pruefbericht.model_dump_json(indent=None, by_alias=True))


@generate.command("from-database")
@grzcli.submission_id
@grzcli.configuration
@fail_or_pass
def from_database(submission_id, configuration, failed, config_file=None):
    """Generate Prüfbericht from database using submission ID."""
    try:
        pruefbericht = _generate_pruefbericht_from_database(submission_id, configuration, failed)
        click.echo(pruefbericht.model_dump_json(indent=None, by_alias=True))
    except ValueError as e:
        raise click.ClickException(str(e)) from e


@pruefbericht.command()
@click.option("--pruefbericht-file", type=click.Path(exists=True), required=True, help="Path to pruefbericht file")
@grzcli.submission_id
@grzcli.configuration
@click.option(
    "--token", help="Access token to try instead of requesting a new one.", envvar="GRZ_PRUEFBERICHT_ACCESS_TOKEN"
)
@click.option("--print-token", is_flag=True, help="Print obtained access token to stdout.")
@click.option(
    "--allow-redacted-tan-g",
    help="Allow submission of a Prüfbericht with a redacted TAN.",
    is_flag=True,
)
@grzcli.update_db
def submit(  # noqa: PLR0913
    configuration: dict[str, Any],
    pruefbericht_file,
    submission_id,
    token,
    print_token,
    allow_redacted_tan_g,
    update_db,
    **kwargs,
):
    """Submit a Prüfbericht JSON to BfArM."""
    config = PruefberichtConfig.model_validate(configuration)

    with open(pruefbericht_file) as f:
        pruefbericht = Pruefbericht.model_validate_json(f.read())

    if (auth_url := config.pruefbericht.authorization_url) is None:
        raise ValueError("pruefbericht.auth_url must be provided to submit Prüfberichte")
    if (client_id := config.pruefbericht.client_id) is None:
        raise ValueError("pruefbericht.client_id must be provided to submit Prüfberichte")
    if (client_secret := config.pruefbericht.client_secret) is None:
        raise ValueError("pruefbericht.client_secret must be provided to submit Prüfberichte")
    if (api_base_url := config.pruefbericht.api_base_url) is None:
        raise ValueError("pruefbericht.api_base_url must be provided to submit Prüfberichte")

    if pruefbericht.submitted_case.tan == REDACTED_TAN and not allow_redacted_tan_g:
        raise ValueError("Refusing to submit a Prüfbericht with a redacted TAN")

    with DbContext(
        configuration=configuration,
        submission_id=submission_id,
        start_state=SubmissionStateEnum.REPORTING,
        end_state=SubmissionStateEnum.REPORTED,
        enabled=update_db,
    ):
        expiry, token = _try_submit(
            pruefbericht=pruefbericht,
            api_base_url=str(api_base_url),
            auth_url=str(auth_url),
            client_id=client_id,
            client_secret=client_secret,
            token=token,
        )

    log.info("Prüfbericht submitted successfully.")

    if expiry and print_token:
        log.info(f"New token expires at {expiry.isoformat()}")
        click.echo(token)


def _try_submit(  # noqa: PLR0913
    pruefbericht: Pruefbericht, api_base_url: str, auth_url: str, client_id: str, client_secret: str, token: str
) -> tuple[Any, Any]:
    if token:
        # replace newlines in token if accidentally present from pasting
        token = token.replace("\n", "")
        expiry = None
    else:
        token, expiry = _get_new_token(
            auth_url=auth_url,
            client_id=client_id,
            client_secret=client_secret,
        )

    try:
        _submit_pruefbericht(base_url=api_base_url, token=token, pruefbericht=pruefbericht)
    except requests.HTTPError as error:
        if error.response.status_code == requests.codes.unauthorized:
            # get a new token and try again
            log.warning("Provided token has expired. Attempting to refresh.")
            token, expiry = _get_new_token(
                auth_url=auth_url,
                client_id=client_id,
                client_secret=client_secret,
            )
            _submit_pruefbericht(base_url=api_base_url, token=token, pruefbericht=pruefbericht)
        else:
            log.error("Encountered an irrecoverable error while submitting the Prüfbericht!")
            raise error
    return expiry, token
