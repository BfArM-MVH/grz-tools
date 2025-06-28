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

    if not submission_id:
        sys.exit("No submission ID provided. Please specify a submission ID to clean.")

    if yes_i_really_mean_it or click.confirm(
        f"Are you SURE you want to delete the submission '{submission_id}' from the bucket '{bucket_name}'?",
        default=False,
        show_default=True,
    ):
        prefix = submission_id
        prefix = prefix + "/" if not prefix.endswith("/") else prefix

        resource = init_s3_resource(config.s3)
        bucket = resource.Bucket(bucket_name)
        log.info(f"Cleaning '{prefix}' from '{bucket_name}' …")
        # add a marker at start of cleaning to
        #  1.) ensure user can upload the "cleaned" marker at the end _before_ we start deleting things
        #  2.) detect incomplete cleans if needed
        bucket.put_object(Body=b"", Key=f"{submission_id}/cleaning")

        # keep metadata.json to prevent future re-uploads
        keys_to_keep = {f"{submission_id}/metadata/metadata.json", f"{submission_id}/cleaning"}
        num_deleted = 0
        for obj in bucket.objects.filter(Prefix=prefix):
            if obj.key not in keys_to_keep:
                _ = obj.delete()
                num_deleted += 1
        if not num_deleted:
            sys.exit(f"No objects with prefix '{prefix}' in bucket '{bucket_name}' found for deletion.")

        log.info(f"Successfully deleted {num_deleted} objects.")

        # redact metadata.json since it contains tanG + localCaseId
        bucket.put_object(Body=b"", Key=f"{submission_id}/metadata/metadata.json")

        # mark that we've cleaned this submission
        bucket.put_object(Body=b"", Key=f"{submission_id}/cleaned")
        bucket.Object(f"{submission_id}/cleaning").delete()

        log.info(f"Cleaned '{prefix}' from '{bucket_name}'.")
