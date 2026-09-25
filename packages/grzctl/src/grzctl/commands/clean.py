"""Command for cleaning a submission from the S3 inbox."""

import logging
import sys

import click
import grz_common.cli as grzcli
from grz_common.transfer import init_s3_resource, s3_errors
from grz_db.models.submission import SubmissionStateEnum

from ..commands import grzctl_configuration, inbox_option
from ..dbcontext import DbContext
from ..models.config import GrzctlConfig
from .db.cli import get_submission_db_instance
from .inbox_resolution import require_inbox

log = logging.getLogger(__name__)


@click.command()
@grzctl_configuration
@grzcli.submission_id
@click.option("--yes-i-really-mean-it", is_flag=True)
@grzcli.update_db
@inbox_option
def clean(
    configuration: GrzctlConfig,
    submission_id: str,
    yes_i_really_mean_it: bool,
    update_db: bool,
    inbox_name: str,
    **kwargs,
):
    """
    Remove all files of a submission from the S3 inbox.
    """
    submitter_id = submission_id.split("_", maxsplit=1)[0]
    resolved_inbox = require_inbox(
        configuration,
        submitter_id=submitter_id,
        submission_id=submission_id,
        inbox_name=inbox_name,
        db_service=get_submission_db_instance(db_url=configuration.db.database_url) if update_db else None,
    )
    s3_options = configuration.inbox_target(submitter_id=submitter_id, inbox_name=resolved_inbox).s3
    bucket_name = s3_options.bucket
    inbox_desc = f"'{resolved_inbox}' (bucket '{bucket_name}')" if resolved_inbox != bucket_name else f"'{bucket_name}'"

    if not submission_id:
        sys.exit("No submission ID provided. Please specify a submission ID to clean.")

    if yes_i_really_mean_it or click.confirm(
        f"Are you SURE you want to delete the submission '{submission_id}' from inbox {inbox_desc}?",
        default=False,
        show_default=True,
    ):
        with DbContext(
            configuration=configuration,
            submission_id=submission_id,
            start_state=SubmissionStateEnum.CLEANING,
            end_state=SubmissionStateEnum.CLEANED,
            enabled=update_db,
        ):
            _clean_submission_from_bucket(bucket_name, s3_options, submission_id, inbox_desc)


def _clean_submission_from_bucket(bucket_name: str, s3_options, submission_id: str, inbox_desc: str):
    prefix = submission_id
    prefix = prefix + "/" if not prefix.endswith("/") else prefix

    resource = init_s3_resource(s3_options)
    bucket = resource.Bucket(bucket_name)
    log.info(f"Cleaning '{prefix}' from inbox {inbox_desc} …")
    with s3_errors(f"Cleaning {submission_id} from inbox {inbox_desc}"):
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
            sys.exit(f"No objects with prefix '{prefix}' in inbox {inbox_desc} found for deletion.")

        log.info(f"Successfully deleted {num_deleted} objects.")

        # redact metadata.json since it contains tanG + localCaseId
        bucket.put_object(Body=b"", Key=f"{submission_id}/metadata/metadata.json")

        # mark that we've cleaned this submission
        bucket.put_object(Body=b"", Key=f"{submission_id}/cleaned")
        bucket.Object(f"{submission_id}/cleaning").delete()

    log.info(f"Cleaned '{prefix}' from inbox {inbox_desc}.")
