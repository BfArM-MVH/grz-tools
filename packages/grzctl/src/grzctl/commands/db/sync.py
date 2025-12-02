"""
Logic for synchronizing submissions from an inbox to the database.
"""

import logging

from grz_common.workers.download import InboxSubmissionState, InboxSubmissionSummary
from grz_db.models.author import Author
from grz_db.models.submission import SubmissionDb, SubmissionStateEnum

log = logging.getLogger(__name__)


def sync_submissions(
    db: SubmissionDb,
    submissions: list[InboxSubmissionSummary],
    author: Author,
) -> None:
    """
    Synchronizes a list of submissions with the database, i.e., registers submissions in the database and sets their
    state to "uploading" or "uploaded" respectively.

    Args:
        db: The database service instance.
        submissions: List of submission summaries from the S3 inbox.
        author: The author to sign the state changes with.
    """
    if not submissions:
        return

    all_db_submissions = db.list_submissions(limit=None)
    current_db_states: dict[str, str] = {}
    for sub in all_db_submissions:
        latest = sub.get_latest_state()
        if latest:
            current_db_states[sub.id] = latest.state.value.lower()
        else:
            current_db_states[sub.id] = "missing"

    for submission in submissions:
        submission_id = submission.submission_id
        s3_state = submission.state

        target_db_state = _determine_target_state(s3_state)

        # we only update to "uploading" or "uploaded", everything else is not of relevance here.
        if not target_db_state:
            log.debug(f"Skipping submission {submission_id} with S3 state {s3_state}.")
            continue

        current_db_state = current_db_states.get(submission_id)

        match (current_db_state, target_db_state):
            # if the submission is new (i.e., not in DB yet), we add it and set the state based on S3 status.
            case (None, target) if target in (SubmissionStateEnum.UPLOADING, SubmissionStateEnum.UPLOADED):
                log.info(f"Submission {submission_id} is new. Adding to DB.")
                try:
                    db.add_submission(submission_id)
                    _update_state(db, submission_id, target, author)
                except Exception as e:
                    log.error(f"Failed to add/update new submission {submission_id}: {e}")

            # if the submission is "uploading" in DB, but "complete" in S3, update to "uploaded"
            case ("uploading", SubmissionStateEnum.UPLOADED):
                log.info(f"Updating state for {submission_id} from 'uploading' to 'uploaded'.")
                try:
                    _update_state(db, submission_id, SubmissionStateEnum.UPLOADED, author)
                except Exception as e:
                    log.error(f"Failed to update submission {submission_id}: {e}")

            # … otherwise just skip.
            case _:
                pass


def _determine_target_state(s3_state: InboxSubmissionState) -> SubmissionStateEnum | None:
    match s3_state:
        case InboxSubmissionState.COMPLETE:
            return SubmissionStateEnum.UPLOADED
        case InboxSubmissionState.INCOMPLETE:
            return SubmissionStateEnum.UPLOADING
        case _:
            return None


def _update_state(db: SubmissionDb, submission_id: str, state: SubmissionStateEnum, author: Author):
    db.update_submission_state(submission_id, state)
