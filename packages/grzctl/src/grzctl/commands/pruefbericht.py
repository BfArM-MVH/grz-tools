"""Command for submitting Prüfberichte."""

import datetime
import logging
from http import HTTPStatus
from typing import Any

import click
import grz_common.cli as grzcli
import requests
from grz_common.exceptions import (
    ConfigurationError,
    GrzError,
    NetworkError,
    PruefberichtGenerationError,
    PruefberichtRejectedError,
)
from grz_common.workers.submission import Submission
from grz_db.models.submission import SubmissionDb, SubmissionStateEnum
from grz_pydantic_models.pruefbericht.v0 import LibraryType as PruefberichtLibraryType
from grz_pydantic_models.pruefbericht.v0 import Pruefbericht, SubmittedCase
from grz_pydantic_models.submission.metadata.v1 import REDACTED_TAN, GrzSubmissionMetadata, Relation
from pydantic_core import to_jsonable_python

from ..commands import grzctl_configuration
from ..dbcontext import DbContext
from ..models.config import GrzctlConfig
from ..models.pruefbericht import PruefberichtModel

log = logging.getLogger(__name__)
fail_or_pass = click.option(
    "--fail/--pass", "failed", help="Fail an otherwise valid submission (e.g. failed internal QC)"
)


def _http_error(error: requests.RequestException, client_error: type[GrzError]) -> GrzError:
    """Classify a failed HTTP request as the failure it stands for.

    :param error: What ``requests`` raised. Its message names the URL of the request.
    :param client_error: The class for a client error other than refused credentials.
    :returns: A :class:`NetworkError` if the request did not get through or the server answered with a
        server error, a :class:`ConfigurationError` if the server refuses the credentials, and
        ``client_error`` otherwise.
    """
    status = error.response.status_code if error.response is not None else None
    if status is None or status >= HTTPStatus.INTERNAL_SERVER_ERROR:
        error_class: type[GrzError] = NetworkError
    elif status in {HTTPStatus.UNAUTHORIZED, HTTPStatus.FORBIDDEN}:
        error_class = ConfigurationError
    else:
        error_class = client_error
    return error_class(str(error))


def _get_new_token(auth_url: str, client_id: str, client_secret: str) -> tuple[str, datetime.datetime]:
    log.info("Refreshing access token...")

    try:
        response = requests.post(
            auth_url,
            headers={"Content-Type": "application/x-www-form-urlencoded"},
            data={"grant_type": "client_credentials", "client_id": client_id, "client_secret": client_secret},
            timeout=60,
        )
        if response.status_code != HTTPStatus.OK:
            log.error("There was a problem refreshing the access token")
            response.raise_for_status()
    except requests.RequestException as e:
        # the token request sends only the client credentials, so BfArM can reject nothing else
        raise _http_error(e, ConfigurationError) from e

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

    if response.status_code != HTTPStatus.OK:
        log.warning("There was a problem submitting the Prüfbericht.")
        response.raise_for_status()


def get_pruefbericht_library_type(metadata: GrzSubmissionMetadata) -> PruefberichtLibraryType:
    """
    Determine the singular representative library type of a submission to submit with the Prüfbericht.
    This should be library type of the index patient with the highest reimbursement value.
    """
    index_patient = metadata.index_donor
    index_patient_submission_library_types = {str(datum.library_type) for datum in index_patient.lab_data}
    return PruefberichtLibraryType.most_expensive(index_patient_submission_library_types)


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


def _generate_pruefbericht_from_database(submission_id: str, configuration: GrzctlConfig, failed: bool) -> Pruefbericht:
    """Generate Prüfbericht by fetching submission data from the database.

    :raises PruefberichtGenerationError: If the database lacks the submission or data that the Prüfbericht needs.
    """
    db = configuration.db

    db_service = SubmissionDb(db_url=str(db.database_url), author=None, debug=False)
    submission = db_service.get_submission(submission_id)

    if submission is None:
        raise PruefberichtGenerationError(f"Submission with ID '{submission_id}' not found in database")

    # Check if submission has the required fields populated
    required_fields = [
        "submission_uploaded_date",
        "submission_type",
        "tan_g",
        "submitter_id",
        "data_node_id",
        "disease_type",
        "coverage_type",
    ]

    missing_fields = [field for field in required_fields if getattr(submission, field) is None]
    if missing_fields:
        raise PruefberichtGenerationError(
            f"Submission {submission_id} is missing required fields: {', '.join(missing_fields)}"
        )

    # Get donors to determine library types
    donors = db_service.get_donors(submission_id)
    if not donors:
        raise PruefberichtGenerationError(f"No donors found for submission {submission_id}")

    # Find index donor
    index_donor = next((d for d in donors if d.relation == Relation.index_), None)
    if index_donor is None:
        raise PruefberichtGenerationError(f"No index donor found for submission {submission_id}")

    # Convert database library_types to strings and determine most expensive type
    index_donor_library_types = {str(lt.value if hasattr(lt, "value") else lt) for lt in index_donor.library_types}

    try:
        library_type = PruefberichtLibraryType.most_expensive(index_donor_library_types)

        # Generate the Prüfbericht
        return Pruefbericht(
            SubmittedCase=SubmittedCase(
                submissionDate=submission.submission_uploaded_date,
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
    except ValueError as e:
        raise PruefberichtGenerationError(f"The Prüfbericht of {submission_id} cannot be generated: {e}") from e


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
@grzctl_configuration
@fail_or_pass
def from_database(submission_id, configuration: GrzctlConfig, failed):
    """Generate Prüfbericht from database using submission ID."""
    try:
        pruefbericht = _generate_pruefbericht_from_database(submission_id, configuration, failed)
        click.echo(pruefbericht.model_dump_json(indent=None, by_alias=True))
    except PruefberichtGenerationError as e:
        raise click.ClickException(str(e)) from e


@pruefbericht.command()
@click.option("--pruefbericht-file", type=click.Path(exists=True), required=True, help="Path to pruefbericht file")
@grzcli.submission_id
@grzctl_configuration
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
def submit(  # noqa: PLR0913, PLR0917
    configuration: GrzctlConfig,
    pruefbericht_file,
    submission_id,
    token,
    print_token,
    allow_redacted_tan_g,
    update_db,
    **kwargs,
):
    """Submit a Prüfbericht JSON to BfArM."""
    with open(pruefbericht_file) as f:
        pruefbericht = Pruefbericht.model_validate_json(f.read())

    auth_url, client_id, client_secret, api_base_url = _get_submission_credentials(configuration.pruefbericht)

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
            api_base_url=api_base_url,
            auth_url=auth_url,
            client_id=client_id,
            client_secret=client_secret,
            token=token,
        )

    log.info("Prüfbericht submitted successfully.")

    if expiry and print_token:
        log.info(f"New token expires at {expiry.isoformat()}")
        click.echo(token)


def _get_submission_credentials(pb: PruefberichtModel) -> tuple[str, str, str, str]:
    """Return ``(auth_url, client_id, client_secret, api_base_url)``, or raise if one is not configured."""
    if (auth_url := pb.authorization_url) is None:
        raise ConfigurationError("pruefbericht.authorization_url must be provided to submit Prüfberichte")
    if (client_id := pb.client_id) is None:
        raise ConfigurationError("pruefbericht.client_id must be provided to submit Prüfberichte")
    if (client_secret := pb.client_secret) is None:
        raise ConfigurationError("pruefbericht.client_secret must be provided to submit Prüfberichte")
    if (api_base_url := pb.api_base_url) is None:
        raise ConfigurationError("pruefbericht.api_base_url must be provided to submit Prüfberichte")
    return str(auth_url), client_id, client_secret.get_secret_value(), str(api_base_url)


def _try_submit(  # noqa: PLR0913, PLR0917
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
        try:
            _submit_pruefbericht(base_url=api_base_url, token=token, pruefbericht=pruefbericht)
        except requests.HTTPError as error:
            if error.response is None or error.response.status_code != HTTPStatus.UNAUTHORIZED:
                raise
            # get a new token and try again
            log.warning("Provided token has expired. Attempting to refresh.")
            token, expiry = _get_new_token(
                auth_url=auth_url,
                client_id=client_id,
                client_secret=client_secret,
            )
            _submit_pruefbericht(base_url=api_base_url, token=token, pruefbericht=pruefbericht)
    except requests.RequestException as e:
        log.error("Submitting the Prüfbericht failed.")
        raise _http_error(e, PruefberichtRejectedError) from e
    return expiry, token
