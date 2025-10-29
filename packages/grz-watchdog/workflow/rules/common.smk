import datetime
import math
import os
import random
import shutil
from collections.abc import Callable
from operator import itemgetter
from os import PathLike
from typing import Literal

import humanfriendly
import yaml
from grz_cli.commands.validate import validate
from grz_db.models.submission import SubmissionDb
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


def cleanup_stale_temp_outputs():
    """
    For use with the `onerror` hook: Make sure to remove relevant temp files (scan_inbox and sync_database outputs),
    so subsequent runs are forced to re-run scan_inbox and sync_database (if they need their outputs).

    For example, this covers the following case, where:
     - invalid submission A is processed, the workflow fails, is terminated (but the temp files aren't removed automatically)
     - a subsequent valid submission B is processed, which needs to re-run scan_inbox and sync_database, otherwise submission B can't be found because the old scan_inbox / sync_database outputs may only include submission A.
    """
    from snakemake.logging import logger

    logger.info("Checking for stale temporary state files...")
    relevant_temp_outputs = [
        str(rules.scan_inbox.output.submissions),
        str(rules.sync_database.output.submissions),
    ]
    if hasattr(rules, "daemon_keepalive"):
        relevant_temp_outputs.append(str(rules.daemon_keepalive.output.marker))

    paths_to_check = set()

    for template in relevant_temp_outputs:
        for submitter_id, inbox in ALL_INBOX_PAIRS:
            path = template.format(submitter_id=submitter_id, inbox=inbox)
            paths_to_check.add(path)

    for path in paths_to_check:
        path = Path(path)
        if path.exists():
            logger.info(f"Removing stale file: {path}")
            try:
                if path.is_dir():
                    shutil.rmtree(path)
                else:
                    path.unlink()
            except OSError as e:
                logger.error(f"Failed to remove {path}: {e}")
    logger.debug("Cleanup stale temp outputs done")


## RULE INPUT FUNCTIONS


def all_pending_submissions(wildcards: Wildcards):
    batch_file = checkpoints.select_submissions.get().output["submissions_batch"]
    with open(batch_file) as f:
        targets = [line.strip() for line in f if line.strip()]
    return targets


def get_submission_to_sync(wildcards: Wildcards):
    if hasattr(wildcards, "submission_id"):
        return [rules.filter_single_submission.output.filtered_submission]
    else:
        return expand(
            "<results>/scan_inbox/{submitter}/{inbox}/submissions.json",
            zip,
            submitter=map(itemgetter(0), ALL_INBOX_PAIRS),
            inbox=map(itemgetter(1), ALL_INBOX_PAIRS),
        )


def get_cleanup_prerequisite(wildcards: Wildcards):
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


def should_run_qc(
    db_config_path: str, submitter_id: str, target_percentage: float, salt: str
) -> bool:
    """
    Determines if a validated submission should undergo QC.
    """
    with open(db_config_path) as f:
        db_url = yaml.safe_load(f)["db"]["database_url"]
    db = SubmissionDb(db_url=db_url, author=None)

    today = datetime.date.today()
    num_validated, num_qcing, num_since_last_qcing = db.get_monthly_qc_stats(
        submitter_id=submitter_id, year=today.year, month=today.month
    )

    # always qc the first (validated) submission each month
    if num_validated == 0:
        return True

    # for example, for a target_percentage of 2%, pick 1 submission at random from 50 submissions
    target_fraction = target_percentage / 100.0
    block = math.floor(1 / target_fraction)
    block_index = num_qcing
    current_index_in_block = num_since_last_qcing

    # ensure the seed is consistent across submitter, date, block_index and (secret) salt
    seed = f"{submitter_id}-{today.year}-{today.month}-{block_index}-{salt}"
    rng = random.Random(seed)  # noqa: S311
    target_index_in_block = rng.randint(0, block - 1)

    if current_index_in_block == target_index_in_block:
        return True

    # if we somehow have exceeded the block size, return True, otherwise False
    return current_index_in_block >= block - 1


def get_final_submission_target(wildcards: Wildcards):
    validation_flag_file = checkpoints.validate.get(
        submitter_id=wildcards.submitter_id,
        inbox=wildcards.inbox,
        submission_id=wildcards.submission_id,
    ).output.validation_flag

    with open(validation_flag_file) as f:
        is_valid = f.read().strip() == "true"

    if not is_valid:
        return rules.finalize_fail.output.target.format(
            submitter_id=wildcards.submitter_id,
            inbox=wildcards.inbox,
            submission_id=wildcards.submission_id,
            qc_status="without_qc",
        )

    strategy = config["qc"]["selection_strategy"]
    run_qc = should_run_qc(
        db_config_path=cfg_path("config_paths/db")(wildcards),
        submitter_id=wildcards.submitter_id,
        target_percentage=strategy["target_percentage"],
        salt=strategy["salt"],
    )

    qc_status = "with_qc" if run_qc else "without_qc"

    return rules.finalize_success.output.target.format(
        submitter_id=wildcards.submitter_id,
        inbox=wildcards.inbox,
        submission_id=wildcards.submission_id,
        qc_status=qc_status,
    )


def get_successful_finalize_inputs(wildcards: Wildcards):
    """
    Gather inputs for a successfully validated submission.
    Conditionally include the QC result marker based on the qc_status wildcard.
    """
    inputs = {
        "archived_marker": rules.archive.output.marker,
        "pruefbericht_answer": rules.submit_pruefbericht.output.answer,
        "pruefbericht": rules.generate_pruefbericht.output.pruefbericht,
        "clean_results": rules.clean.output.clean_results,
        "db_config_path": cfg_path("config_paths/db")(wildcards),
    }

    if wildcards.qc_status == "with_qc":
        inputs["qc_processed_marker"] = rules.process_qc_results.output.marker
    return inputs


def get_failed_finalize_inputs(wildcards: Wildcards):
    """Input function for the failure endpoint to conditionally add cleanup."""
    inputs = {
        "db_config_path": cfg_path("config_paths/db")(wildcards),
        "validation_errors": rules.validate.output.validation_errors,
    }
    if config.get("on-failed-validation", "do-nothing") == "cleanup":
        inputs["clean_results"] = rules.clean.output.clean_results
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


def get_qc_workflow_references_directory() -> str:
    return directory(config["qc"].get("reference_directory", "<resources>/references"))


def get_prepare_qc_nextflow_extra_params(wildcards: Wildcards, input, output):
    work_dir = os.path.abspath(output.work_dir)
    extra = config.get("qc", {}).get("prepare-qc", {}).get("extra", "")
    return f"-resume -work-dir {work_dir} {extra}"


def get_run_qc_nextflow_extra_params(wildcards: Wildcards, input, output):
    work_dir = os.path.abspath(output.work_dir)
    extra = config.get("qc", {}).get("run-qc", {}).get("extra", "")
    return f"-resume -work-dir {work_dir} {extra}"


def get_prepare_qc_nextflow_configs(wildcards: Wildcards):
    config_paths = config.get("qc", {}).get("prepare-qc", {}).get("configs", [])
    if not config_paths:
        return ""
    return " ".join([f"-c {os.path.abspath(p)}" for p in config_paths])


def get_prepare_qc_nextflow_profiles(wildcards: Wildcards):
    profiles = config.get("qc", {}).get("prepare-qc", {}).get("profiles", ["conda"])
    return ",".join(profiles)


def get_run_qc_nextflow_configs(wildcards: Wildcards):
    config_paths = config.get("qc", {}).get("run-qc", {}).get("configs", [])
    if not config_paths:
        return ""
    return " ".join([f"-c {os.path.abspath(p)}" for p in config_paths])


def get_run_qc_nextflow_profiles(wildcards: Wildcards):
    profiles = config.get("qc", {}).get("run-qc", {}).get("profiles", ["conda"])
    return ",".join(profiles)


def get_target_qc_percentage(wildcards: Wildcards) -> float | None:
    """
    If automatic qc selection strategy is enabled, returns the target qc percentage.
    Otherwise, returns None.
    """
    selection_strategy = config["qc"].get("selection_strategy", {})
    if bool(selection_strategy.get("enabled", False)):
        return selection_strategy.get("target_qc_percentage", 2.0)
    else:
        return None


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


# if snakemake runtime is given as an int, will be interpreted as _minutes_
def estimated_processing_time_in_minutes(
    data_size_in_bytes: int, processing_speed_in_mb_per_s: int
) -> int:
    if processing_speed_in_mb_per_s <= 0:
        raise ValueError("processing_speed_in_mb_per_s must be positive")
    bytes_per_second = processing_speed_in_mb_per_s * 1024**2
    bytes_per_minute = float(bytes_per_second) / 60
    minutes = int(math.ceil(data_size_in_bytes / bytes_per_minute))
    return minutes


def estimate_runtime(what: str) -> Callable[[Wildcards, InputFiles], int]:
    mb_per_s = config["estimates"]["speed"][what]

    def inner(wildcards: Wildcards, input: InputFiles) -> int:
        return estimated_processing_time_in_minutes(
            dataset_size(parse_metadata(input.metadata)), mb_per_s
        )

    return inner


for step in ("download", "decrypt", "validate", "encrypt", "archive", "qc"):
    locals()[f"estimate_{step}_runtime"] = estimate_runtime(step)


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


def estimate_qc_disk(wildcards: Wildcards, input: InputFiles) -> str:
    bytes = dataset_size(parse_metadata(input.metadata))
    return humanfriendly.format_size(bytes * 6)
