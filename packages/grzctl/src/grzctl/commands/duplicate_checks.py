"""Checks of a submission's metadata against the submissions already in the database.

Both checks need nothing but the metadata and the database, so they can run as soon as the
metadata is known.
"""

import logging
from typing import TYPE_CHECKING

from grz_db.errors import DuplicateInitialSubmissionError
from grz_pydantic_models.submission.metadata import GrzSubmissionMetadata

if TYPE_CHECKING:
    from grz_db.models.submission import SubmissionDb

log = logging.getLogger(__name__)


def check_duplicate_initial(
    db: "SubmissionDb", metadata: GrzSubmissionMetadata, submission_id: str | None = None
) -> None:
    """Raise :class:`DuplicateInitialSubmissionError` if this initial submission is a duplicate.

    The rule and the index enforcing it belong to the database layer, which answers this in
    one transaction; all that is left here is unpacking the metadata. Validation is what the
    answer buys: the database would reject the submission anyway, but only once basic QC is
    being recorded, by which point the effort has been spent.

    A resolution failure is left to propagate. Reporting "not a duplicate" when the
    question could not be answered is how a duplicate would pass basic QC unnoticed.
    Nothing later in this function asks again, so raising here is the only chance to notice.

    :param db: The submission database.
    :param metadata: Metadata of the submission to check.
    :param submission_id: ID the submission is stored under. Defaults to the one the metadata derives.
    """
    submission = metadata.submission
    db.assert_no_duplicate_initial(
        submission_id or metadata.submission_id,
        submitter_id=submission.submitter_id,
        local_case_id=submission.local_case_id,
        submission_type=submission.submission_type,
    )


def reject_duplicates(db: "SubmissionDb", submission_id: str, metadata: GrzSubmissionMetadata) -> None:
    """Fail a submission whose tanG or case is taken, before its files are transferred.

    The checks repeat what populating and validation enforce later. They are no replacement: a
    competing initial submission can still pass basic QC in between, and the database indices
    stay the final authority. Asking early only spares the download and decryption of a
    submission that is going to be rejected anyway.

    A duplicate initial submission is recorded as failing basic QC, exactly as validation does.

    :param db: The submission database.
    :param submission_id: ID the submission is stored under.
    :param metadata: Metadata of the submission, as read from the inbox.
    :raises DuplicateTanGError: if another submission holds the tanG.
    :raises DuplicateInitialSubmissionError: if the case already has a QC-passed initial submission.
    """
    db.assert_no_duplicate_tan_g(submission_id, metadata.submission.tan_g)
    try:
        check_duplicate_initial(db, metadata, submission_id)
    except DuplicateInitialSubmissionError as e:
        log.warning(f"{e} Failing basic QC for submission '{submission_id}' before downloading its files.")
        _ = db.modify_submission(submission_id, "basic_qc_passed", "false")
        raise
