"""Command for cleaning a submission from the S3 inbox."""

import logging
import sys

import click
from grz_common.cli import config_file, submission_id
from grz_common.transfer import init_s3_resource

from ..models.config import CleanConfig

log = logging.getLogger(__name__)


@click.command()
@submission_id
@config_file
@click.option("--yes-i-really-mean-it", is_flag=True)
def clean(submission_id, config_file, yes_i_really_mean_it: bool):
    """
    Remove all files of a submission from the S3 inbox.
    """
    config = CleanConfig.from_path(config_file)
    bucket_name = config.s3.bucket

    if yes_i_really_mean_it or click.confirm(
        f"Are you SURE you want to delete the submission from the inbox bucket ({bucket_name})?",
        default=False,
        show_default=True,
    ):
        prefix = submission_id
        prefix = prefix + "/" if not prefix.endswith("/") else prefix

        resource = init_s3_resource(config.s3)
        bucket = resource.Bucket(bucket_name)
        log.info(f"Deleting {prefix} from {bucket_name} …")

        responses = bucket.objects.filter(Prefix=prefix).delete()
        if not responses:
            sys.exit(f"No objects with prefix '{prefix}' in bucket '{bucket_name}' found for deletion.")

        successfully_deleted_keys = []
        errors_encountered = []
        objects_found = 0

        # responses is a list of dicts reporting the result of the deletion API call
        # because the API does things in batches
        for response in responses:
            deleted_batch = response.get("Deleted", [])
            for deleted_obj in deleted_batch:
                key = deleted_obj.get("Key")
                if key:
                    successfully_deleted_keys.append(key)

            errors_batch = response.get("Errors", [])
            for error_obj in errors_batch:
                errors_encountered.append(
                    {
                        "Key": error_obj.get("Key", "N/A"),
                        "Code": error_obj.get("Code", "N/A"),
                        "Message": error_obj.get("Message", "N/A"),
                    }
                )

            objects_found += len(deleted_batch) + len(errors_batch)

        log.info(f"Total objects attempted to delete: {objects_found}")
        log.info(f"Successfully deleted: {len(successfully_deleted_keys)} objects.")
        log.info(f"Failed to delete: {len(errors_encountered)} objects.")

        for error in errors_encountered:
            log.error(f"  - Key: {error['Key']}, Code: {error['Code']}, Message: {error['Message']}")

        if errors_encountered:
            sys.exit(f"Errors encountered while deleting objects from bucket '{bucket_name}'. See log for details.")

        log.info(f"Deleted '{prefix}' from '{bucket_name}'.")
