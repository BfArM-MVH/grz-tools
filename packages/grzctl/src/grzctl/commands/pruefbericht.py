"""Command for submitting Prüfberichte."""

import datetime
import json
import logging
import sys

import click
import requests
from grz_common.cli import config_file, output_json, submission_dir
from grz_common.workers.submission import Submission

# pyrefly: ignore
from grz_pydantic_models.pruefbericht import LibraryType, Pruefbericht, SubmittedCase

# pyrefly: ignore
from grz_pydantic_models.submission.metadata.v1 import (
    GenomicStudySubtype,
    GrzSubmissionMetadata,
    Relation,
    SequenceSubtype,
)
from pydantic_core import to_jsonable_python

from ..models.config import PruefberichtConfig

log = logging.getLogger(__name__)


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


def _get_library_type(metadata: GrzSubmissionMetadata) -> LibraryType:
    index_patients = [donor for donor in metadata.donors if donor.relation == Relation.index_]
    if len(index_patients) != 1:
        # TODO handle in grz-pydantic-model validation?
        raise ValueError("Multiple index patients in one submission")
    index_patient = index_patients[0]

    # assume somatic is tumor in context of genomDE
    data_tumor = filter(lambda datum: datum.sequence_subtype == SequenceSubtype.somatic, index_patient.lab_data)
    data_germline = filter(lambda datum: datum.sequence_subtype == SequenceSubtype.germline, index_patient.lab_data)

    unique_tumor_library_types = {datum.library_type for datum in data_tumor}
    if len(unique_tumor_library_types) > 1:
        raise ValueError("Multiple different library types detected for tumor data")

    unique_germline_library_types = {datum.library_type for datum in data_germline}
    if len(unique_germline_library_types) > 1:
        raise ValueError("Multiple different library types detected for germline data")

    match metadata.submission.genomic_study_subtype:
        case GenomicStudySubtype.tumor_only:
            if not unique_tumor_library_types:
                raise ValueError("No somatic library types detected in tumor submission")
            library_type = LibraryType(str(next(iter(unique_tumor_library_types))))
        case GenomicStudySubtype.germline_only:
            if not unique_germline_library_types:
                raise ValueError("No germline library types detected in germline submission")
            library_type = LibraryType(str(next(iter(unique_germline_library_types))))
        case GenomicStudySubtype.tumor_germline:
            unique_library_types = unique_tumor_library_types | unique_germline_library_types
            if len(unique_library_types) != 1:
                # TODO: is there a better solution to this?
                raise ValueError("None or multiple different library types detected for tumor+germline sample")
            library_type = LibraryType(str(next(iter(unique_library_types))))
        case _:
            raise ValueError(f"Unknown genomic study subtype: {metadata.submission.genomic_study_subtype}")

    return library_type


@click.command()
@config_file
@submission_dir
@output_json
@click.option("--fail/--pass", "failed", help="Fail an otherwise valid submission (e.g. failed internal QC)")
@click.option(
    "--token", help="Access token to try instead of requesting a new one.", envvar="GRZ_PRUEFBERICHT_ACCESS_TOKEN"
)
def pruefbericht(config_file, submission_dir, output_json, failed, token):
    """
    Submit a Prüfbericht to BfArM.
    """
    config = PruefberichtConfig.from_path(config_file)

    if config.pruefbericht.authorization_url is None:
        raise ValueError("pruefbericht.auth_url must be provided to submit Prüfberichte")
    if config.pruefbericht.client_id is None:
        raise ValueError("pruefbericht.client_id must be provided to submit Prüfberichte")
    if config.pruefbericht.client_secret is None:
        raise ValueError("pruefbericht.client_secret must be provided to submit Prüfberichte")

    submission = Submission(metadata_dir=f"{submission_dir}/metadata", files_dir=f"{submission_dir}/files")

    metadata = submission.metadata.content
    pruefbericht = Pruefbericht(
        SubmittedCase=SubmittedCase(
            submissionDate=metadata.submission.submission_date,
            submissionType=metadata.submission.submission_type,
            tan=metadata.submission.tan_g,
            submitterId=metadata.submission.submitter_id,
            dataNodeId=metadata.submission.genomic_data_center_id,
            diseaseType=metadata.submission.disease_type,
            dataCategory="genomic",
            libraryType=_get_library_type(metadata),
            coverageType=metadata.submission.coverage_type,
            dataQualityCheckPassed=not failed,
        )
    )

    if token:
        # replace newlines in token if accidentally present from pasting
        token = token.replace("\n", "")
        expiry = None
    else:
        token, expiry = _get_new_token(
            auth_url=str(config.pruefbericht.authorization_url),
            client_id=config.pruefbericht.client_id,
            client_secret=config.pruefbericht.client_secret,
        )

    try:
        _submit_pruefbericht(base_url=str(config.pruefbericht.api_base_url), token=token, pruefbericht=pruefbericht)
    except requests.HTTPError as error:
        if error.response.status_code == requests.codes.unauthorized:
            # get a new token and try again
            log.warning("Provided token has expired. Attempting to refresh.")
            token, expiry = _get_new_token(
                auth_url=str(config.pruefbericht.authorization_url),
                client_id=config.pruefbericht.client_id,
                client_secret=config.pruefbericht.client_secret,
            )
            _submit_pruefbericht(base_url=str(config.pruefbericht.api_base_url), token=token, pruefbericht=pruefbericht)
        else:
            log.error("Encountered an irrecoverable error while submitting the Prüfbericht!")
            raise error

    log.info("Prüfbericht submitted successfully.")

    if output_json and expiry:
        json.dump({"token": token, "expires": expiry.isoformat()}, sys.stdout)
    elif expiry:
        log.info(f"New token expires at {expiry.isoformat()}")
        print(token)
