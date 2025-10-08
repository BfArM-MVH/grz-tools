import shutil
from operator import itemgetter
from os import PathLike
from typing import Literal

import humanfriendly
import yaml
from grz_pydantic_models.submission.metadata import GrzSubmissionMetadata
from snakemake.io import InputFiles, Wildcards

cfg_path = lambda dpath: lambda wildcards: (
    lookup(dpath=dpath, within=config)(wildcards)
    if "{" in dpath
    else lookup(dpath=dpath, within=config)
)

ALL_INBOX_PAIRS = [
    (submitter, inbox)
    for submitter, inboxes in config["config_paths"]["inbox"].items()
    for inbox in inboxes.keys()
]

## RULE INPUT FUNCTIONS


def all_pending_submissions(wildcards):
    batch_file = checkpoints.select_submissions.get().output["submissions_batch"]
    with open(batch_file) as f:
        targets = [line.strip() for line in f if line.strip()]
    return targets


def get_submission_to_sync(wildcards):
    if hasattr(wildcards, "submission_id"):
        return [rules.filter_single_submission.output.filtered_submission]
    else:
        return expand(
            "results/scan_inbox/{submitter}/{inbox}/submissions.json",
            zip,
            submitter=map(itemgetter(0), ALL_INBOX_PAIRS),
            inbox=map(itemgetter(1), ALL_INBOX_PAIRS),
        )


def get_cleanup_prerequisite(wildcards):
    """
    Determines the prerequisite for the 'clean' rule.
    - If validation was successful, it depends on the downstream QC/reporting rules.
    - If validation failed, it depends directly on the validation flag.
    """
    validation_flag_file = checkpoints.validate.get(
        submitter_id=wildcards.submitter_id,
        inbox=wildcards.inbox,
        submission_id=wildcards.submission_id,
    ).output.validation_flag

    with open(validation_flag_file) as f:
        is_valid = f.read().strip() == "true"

    if is_valid:
        if wildcards.qc_status == "with_qc":
            return rules.process_qc_results.output.marker
        else:
            return rules.submit_pruefbericht.output.answer
    else:
        return validation_flag_file


def get_final_submission_target(wildcards):
    validation_flag_file = checkpoints.validate.get(
        submitter_id=wildcards.submitter_id,
        inbox=wildcards.inbox,
        submission_id=wildcards.submission_id,
    ).output.validation_flag

    with open(validation_flag_file) as f:
        is_valid = f.read().strip() == "true"

    if is_valid:
        return rules.finalize_success.output.target
    else:
        return rules.finalize_fail.output.target


def get_failed_finalize_inputs(wildcards):
    """Input function for the failure endpoint to conditionally add cleanup."""
    inputs = {
        "db_config_path": cfg_path("config_paths/db")(wildcards),
        "validation_errors": rules.validate.output.validation_errors,
    }
    if config.get("on-failed-validation", "do-nothing") == "cleanup":
        inputs["cleanup_done"] = rules.clean.output.clean_results
    return inputs


## PARAMETER / HELPER FUNCTIONS


def perhaps_temp(f):
    if config.get("temp-outputs", False):
        return f
    else:
        return temp(f)


def get_endpoint_url(wildcards: Wildcards, input: InputFiles) -> str:
    """
    Get the endpoint URL from the inbox config file.

    Args:
        wildcards: Unused
        input: InputFiles with attribute `inbox_config_path` pointing to the inbox config file.

    Returns:
        A string with describing an endpoint URL, e.g.
        "https://localhost:port".
    """
    with open(input.inbox_config_path) as f:
        endpoint_url = yaml.safe_load(f)["s3"]["endpoint_url"].rstrip("/")
        return endpoint_url


def get_s3_bucket(wildcards: Wildcards, input: InputFiles) -> str:
    """
    Get the S3 bucket name from the inbox config file.
    """
    with open(input.inbox_config_path) as f:
        bucket = yaml.safe_load(f)["s3"]["bucket"]
        return bucket


def get_s3_metadata_key(wildcards: Wildcards) -> str:
    """
    Get the key pointing to the metadata.json file in S3.

    Args:
        wildcards: Wildcards with attributes `inbox` and `submission_id`

    Returns:
        A string with the S3 key pointing to the metadata.json file,
        i.e., "{inbox}/{submission_id}/metadata/metadata.json"
    """
    return f"{wildcards.submission_id}/metadata/metadata.json"


def register_s3_access_key(
    wildcards: Wildcards, input: InputFiles
) -> Literal["success"]:
    """
    Export the S3 access key in the environment as AWS_ACCESS_KEY_ID.

    Try looking up environment variable `GRZ_S3__ACCESS_KEY` first,
    then look up the access key in the inbox config file.

    Args:
        wildcards: Unused
        input: InputFiles with attribute `inbox_config_path` pointing to the inbox config file.

    Returns:
        "success" if the access key was registered, else raises ValueError.

    Raises:
        ValueError: If no S3 access key is found.
    """
    if access_key := os.environ.get("GRZ_S3__ACCESS_KEY", ""):
        os.environ["AWS_ACCESS_KEY_ID"] = access_key
        return "success"

    with open(input.inbox_config_path) as f:
        access_key = yaml.safe_load(f).get("s3", {}).get("access_key", "")
        if not access_key:
            raise ValueError("No S3 access_key found.")
        os.environ["AWS_ACCESS_KEY_ID"] = access_key
        return "success"


def register_s3_secret(wildcards: Wildcards, input: InputFiles) -> Literal["success"]:
    """
    Export the S3 secret in the environment as AWS_SECRET_ACCESS_KEY.

    Try looking up environment variable `GRZ_S3__SECRET` first,
    then look up the secret in the inbox config file.

    Args:
        wildcards: Unused
        input: InputFiles with attribute `inbox_config_path` pointing to the inbox config file.

    Returns:
        "success" if the access key was registered, else raises ValueError.

    Raises:
        ValueError: If no S3 secret is found.
    """
    if secret := os.environ.get("GRZ_S3__SECRET", ""):
        os.environ["AWS_SECRET_ACCESS_KEY"] = secret
        return "success"

    with open(input.inbox_config_path) as f:
        secret = yaml.safe_load(f).get("s3", {}).get("secret", "")
        if not secret:
            raise ValueError("No S3 secret found.")
        os.environ["AWS_SECRET_ACCESS_KEY"] = secret
        return "success"


def get_qc_workflow_revision(wildcards: Wildcards) -> str:
    return config["qc"]["revision"]


## RESOURCE ESTIMATION FUNCTIONS


def parse_metadata(metadata_path: PathLike | str) -> GrzSubmissionMetadata:
    with open(metadata_path) as f:
        json_str = f.read()
        metadata: GrzSubmissionMetadata = GrzSubmissionMetadata.model_validate_json(
            json_str
        )
    return metadata


def dataset_size(metadata: GrzSubmissionMetadata) -> int:
    bytes = 0
    for donor in metadata.donors:
        for lab_datum in donor.lab_data:
            if sequence_data := lab_datum.sequence_data:
                for file in sequence_data.files:
                    bytes += file.file_size_in_bytes
    return bytes


def estimate_download_size(wildcards: Wildcards, input: InputFiles) -> str:
    """
    Estimate the total size of the files to be downloaded.

    Sums all file_size_in_bytes from all sequence_data.files.

    Args:
        wildcards: Unused.
        input: InputFiles with attribute `metadata` pointing to the metadata.json file.

    Returns:
        A humanfriendly string to be used with snakemake's "disk"
        resource (not "disk_mb"!)
    """
    bytes = dataset_size(parse_metadata(input.metadata))
    return humanfriendly.format_size(bytes)


def estimate_decrypt_size(wildcards: Wildcards, input: InputFiles) -> str:
    """
    Estimate the total size of the files to be decrypted (plus the original files).

    Args:
        wildcards: Unused.
        input: InputFiles with attribute `data` pointing to the data directory and `metadata` pointing to the metadata.json file.

    Returns:
        A humanfriendly string to be used with snakemake's "disk"
        resource (not "disk_mb"!)
    """
    bytes = dataset_size(parse_metadata(input.metadata))
    # TODO: check: returning double the size here because we double the files?
    return humanfriendly.format_size(2 * bytes)


def estimate_re_encrypt_size(wildcards: Wildcards, input: InputFiles) -> str:
    """
    Estimate the total size of the files to be re-encrypted (plus the original files).

    Args:
        wildcards: Unused.
        input: InputFiles with attribute `data` pointing to the data directory and `metadata` pointing to the metadata.json file.

    Returns:
        A humanfriendly string to be used with snakemake's "disk"
        resource (not "disk_mb"!)
    """
    bytes = dataset_size(parse_metadata(input.metadata))
    # TODO: check: returning double the size here because we double the files?
    return humanfriendly.format_size(2 * bytes)
