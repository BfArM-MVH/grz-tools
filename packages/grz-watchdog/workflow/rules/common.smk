import shutil
from typing import Literal

import humanfriendly
import yaml
from grz_pydantic_models.submission.metadata import GrzSubmissionMetadata
from snakemake.io import Wildcards, InputFiles


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
    with open(input.metadata, "rt") as f:
        json_str = f.read()
        metadata: GrzSubmissionMetadata = GrzSubmissionMetadata.model_validate_json(
            json_str
        )
    bytes = 0
    for donor in metadata.donors:
        for lab_datum in donor.lab_data:
            if sequence_data := lab_datum.sequence_data:
                for file in sequence_data.files:
                    bytes += file.file_size_in_bytes
    return humanfriendly.format_size(bytes)


def estimate_decrypt_size(wildcards: Wildcards, input: InputFiles) -> str:
    """
    Estimate the total size of the files to be decrypted (plus the original files).

    Delegates to shutil.disk_usage on the downloaded/encrypted files.

    Args:
        wildcards: Unused.
        input: InputFiles with attribute `data` pointing to the data directory.

    Returns:
        A humanfriendly string to be used with snakemake's "disk"
        resource (not "disk_mb"!)
    """
    _total_bytes, used_bytes, _free_bytes = shutil.disk_usage(Path(input.data))
    # TODO: check: returning double the size here because we double the files?
    return humanfriendly.format_size(2 * used_bytes)


def estimate_re_encrypt_size(wildcards: Wildcards, input: InputFiles) -> str:
    """
    Estimate the total size of the files to be re-encrypted (plus the original files).

    Delegates to shutil.disk_usage on the downloaded/encrypted files.

    Args:
        wildcards: Unused.
        input: InputFiles with attribute `data` pointing to the data directory.

    Returns:
        A humanfriendly string to be used with snakemake's "disk"
        resource (not "disk_mb"!)
    """
    _total_bytes, used_bytes, _free_bytes = shutil.disk_usage(Path(input.data))
    # TODO: check: returning double the size here because we double the files?
    return humanfriendly.format_size(2 * used_bytes)


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
    with open(input.inbox_config_path, "rt") as f:
        endpoint_url = yaml.safe_load(f)["s3"]["endpoint_url"].rstrip("/")
        return endpoint_url


def get_s3_metadata_key(wildcards: Wildcards) -> str:
    """
    Get the key pointing to the metadata.json file in S3.

    Args:
        wildcards: Wildcards with attributes `inbox` and `submission_id`

    Returns:
        A string with the S3 key pointing to the metadata.json file,
        i.e., "{inbox}/{submission_id}/metadata/metadata.json"
    """

    return f"{wildcards.inbox}/{wildcards.submission_id}/metadata/metadata.json"


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

    with open(input.inbox_config_path, "rt") as f:
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

    with open(input.inbox_config_path, "rt") as f:
        secret = yaml.safe_load(f).get("s3", {}).get("secret", "")
        if not secret:
            raise ValueError("No S3 secret found.")
        os.environ["AWS_SECRET_ACCESS_KEY"] = secret
        return "success"
