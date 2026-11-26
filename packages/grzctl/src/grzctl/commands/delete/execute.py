"""Delete one submission of a BfArM deletion request from storage and scrub it in the database."""

from dataclasses import asdict, dataclass
from typing import Any

from grz_common.models.s3 import S3Options
from grz_db.models.submission import (
    ChangeRequestEnum,
    RequestRawContentType,
    Submission,
    SubmissionDb,
    SubmissionStateEnum,
)

from ... import get_versions
from ...models.config import GrzctlConfig
from ...models.deletion_request import DeletionRequest
from ..inbox_resolution import require_inbox
from .storage import Holdings, archive_holdings, delete_from_archive, delete_from_inbox, inbox_holdings


@dataclass(frozen=True)
class StorageTargets:
    """The buckets that may hold a submission.

    :param archives: Each archive, by its role.
    :param inbox: The submission's inbox, by its name.
    """

    archives: dict[str, S3Options]
    inbox: tuple[str, S3Options]

    @property
    def buckets(self) -> dict[str, S3Options]:
        """All of them, each by its role."""
        return {**self.archives, self.inbox[0]: self.inbox[1]}


class SubmissionDeletionFailedError(Exception):
    """Deleting a submission stopped before it finished; the submission stays in state ``Deleting``."""


def storage_targets(
    configuration: GrzctlConfig, db: SubmissionDb, submission_id: str, submitter_id: str
) -> StorageTargets:
    """Find the buckets that may hold a submission: both archives and its inbox.

    :param configuration: The grzctl configuration.
    :param db: The submission database, which records the submission's inbox.
    :param submission_id: ID of the submission.
    :param submitter_id: The LE that submitted it.
    :raises click.ClickException: if the submission's inbox cannot be told.
    """
    inbox_name = require_inbox(
        configuration,
        submitter_id=submitter_id,
        submission_id=submission_id,
        db_service=db,
        scan=True,
        hint="Record it with grzctl db submission modify, or delete the submission from its inbox by hand.",
    )
    return StorageTargets(
        archives={
            "consented archive": configuration.archives.consented.s3,
            "non-consented archive": configuration.archives.non_consented.s3,
        },
        inbox=(f"inbox '{inbox_name}'", configuration.inbox_target(submitter_id, inbox_name).s3),
    )


def grzctl_versions() -> dict[str, str]:
    """The versions of grzctl and its packages, to record with a state."""
    return {name: version or "unknown" for name, version in get_versions().items()}


def stored_holdings(targets: StorageTargets, submission_id: str) -> dict[str, Holdings]:
    """What each bucket of *targets* holds of a submission that deleting it would remove, by the bucket's role."""
    holdings = {
        location: archive_holdings(s3_options, submission_id) for location, s3_options in targets.archives.items()
    }
    inbox_location, inbox_options = targets.inbox
    holdings[inbox_location] = inbox_holdings(inbox_options, submission_id)
    return holdings


def finish_deletion(
    db: SubmissionDb, submission_id: str, targets: StorageTargets, note: dict[str, Any] | None = None
) -> dict[str, Any]:
    """Delete a submission whose deletion has started from storage, scrub it and record it as ``Deleted``.

    Deletes the submission from the archives that hold it and from its inbox, checks that nothing
    is left, purges its personal data and sets the state ``Deleted`` with what was deleted.
    Stopping anywhere leaves the state as it was, and calling this again continues.

    :param db: The submission database, with an author to sign the entries.
    :param submission_id: ID of the submission, which must have a ``Deleting`` or ``Deleted`` state.
    :param targets: Its buckets, from :func:`storage_targets`.
    :param note: Recorded with the ``Deleted`` state too, e.g. why the deletion was repeated.
    :returns: What was deleted, as recorded with the ``Deleted`` state.
    :raises SubmissionDeletionFailedError: if a bucket still holds some of the submission.
    :raises StorageDeletionError: if S3 refuses to delete some objects.
    """
    storage = [
        deletion
        for location, s3_options in targets.archives.items()
        if (deletion := delete_from_archive(s3_options, location, submission_id)) is not None
    ]
    inbox_location, inbox_options = targets.inbox
    if (deletion := delete_from_inbox(inbox_options, inbox_location, submission_id)) is not None:
        storage.append(deletion)

    if left := {location: held for location, held in stored_holdings(targets, submission_id).items() if held.objects}:
        raise SubmissionDeletionFailedError(
            "; ".join(f"the {location} still holds {held.objects} objects" for location, held in left.items())
        )

    facts = {
        **(note or {}),
        "storage": [asdict(deletion) for deletion in storage],
        "database": asdict(db.purge_personal_data(submission_id)),
    }
    db.update_submission_state(
        submission_id, SubmissionStateEnum.DELETED, data=facts, grzctl_versions=grzctl_versions()
    )
    return facts


def delete_submission(  # noqa: PLR0913
    db: SubmissionDb,
    submission: Submission,
    targets: StorageTargets,
    request: DeletionRequest,
    request_form: tuple[bytes, RequestRawContentType],
    *,
    record_request: bool,
) -> dict[str, Any]:
    """Carry out the deletion request for one submission.

    Records the request as a ``Delete`` change request (unless *record_request* is false, for a
    deletion that is continued), sets the state ``Deleting`` and finishes the deletion with
    :func:`finish_deletion`. Stopping anywhere leaves the submission in ``Deleting``, and calling
    this again continues.

    :param db: The submission database, with an author to sign the entries.
    :param submission: The submission.
    :param targets: Its buckets, from :func:`storage_targets`.
    :param request: The deletion request.
    :param request_form: The request form and its type, kept with the change request.
    :param record_request: Whether to record the request; false when continuing a deletion,
        whose request is recorded already.
    :returns: What was deleted, as recorded with the ``Deleted`` state.
    :raises SubmissionDeletionFailedError: if a bucket still holds some of the submission.
    :raises StorageDeletionError: if S3 refuses to delete some objects.
    """
    if record_request:
        db.add_change_request(
            submission.id,
            ChangeRequestEnum.DELETE,
            requester_name=request.requester_name,
            requester_email=request.requester_email,
            requested_at=request.requested_at,
            request_raw_content=request_form[0],
            request_raw_content_type=request_form[1],
            data={
                "instructed_at": request.instructed_at.isoformat(),
                "submitter_id": request.submitter_id,
                "data_node_id": request.data_node_id,
                "reason": request.reason,
                "comment": request.comment,
            },
        )
    db.update_submission_state(submission.id, SubmissionStateEnum.DELETING, grzctl_versions=grzctl_versions())
    return finish_deletion(db, submission.id, targets)
