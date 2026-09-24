import logging
from collections.abc import Iterator
from functools import cached_property
from typing import Any

import grz_common.exceptions as grzexc
from grz_db.errors import DuplicateInitialSubmissionError, DuplicateTanGError, SubmissionNotFoundError
from grz_db.models.author import Author
from grz_db.models.submission import FailureReasonEnum, SubmissionDb, SubmissionStateEnum

from . import get_versions
from .commands.db.cli import get_submission_db_instance
from .models.config import GrzctlConfig

log = logging.getLogger(__name__)

_FAILURE_REASONS: dict[type[BaseException], FailureReasonEnum] = {
    grzexc.MissingSubmissionFileError: FailureReasonEnum.FILE_NOT_FOUND,
    grzexc.SubmissionValidationError: FailureReasonEnum.VALIDATION_ERROR,
    grzexc.DecryptionError: FailureReasonEnum.DECRYPTION_ERROR,
    grzexc.DuplicateUploadError: FailureReasonEnum.DUPLICATE_TANG,
    DuplicateTanGError: FailureReasonEnum.DUPLICATE_TANG,
    DuplicateInitialSubmissionError: FailureReasonEnum.DUPLICATE_INITIAL,
    grzexc.IncompleteSubmissionError: FailureReasonEnum.INCOMPLETE_SUBMISSION,
    grzexc.SubmissionCleanedError: FailureReasonEnum.SUBMISSION_CLEANED,
    KeyboardInterrupt: FailureReasonEnum.INTERRUPTED,
    grzexc.ConfigurationError: FailureReasonEnum.CONFIGURATION_ERROR,
    grzexc.TransferError: FailureReasonEnum.TRANSFER_ERROR,
    grzexc.EncryptionError: FailureReasonEnum.ENCRYPTION_ERROR,
    grzexc.DetailedQCError: FailureReasonEnum.DETAILED_QC_ERROR,
    grzexc.PruefberichtGenerationError: FailureReasonEnum.PRUEFBERICHT_GENERATION_ERROR,
    grzexc.PruefberichtRejectedError: FailureReasonEnum.PRUEFBERICHT_REJECTED,
}
"""The failure reason of each expected error. Any other exception records ``unknown``."""


def _causes(error: BaseException | None) -> Iterator[BaseException]:
    """Yield ``error``, then its ``__cause__``, then the cause's ``__cause__``, and so on."""
    seen: set[int] = set()
    # a cause can form a cycle, as in ``raise e from e``
    while error is not None and id(error) not in seen:
        seen.add(id(error))
        yield error
        error = error.__cause__


def _classify(error: BaseException | None) -> tuple[FailureReasonEnum, BaseException | None]:
    """Find the failure reason of ``error``, and the exception that decides it.

    The first exception along the causes of ``error`` whose type is mapped decides the reason.

    :returns: The failure reason and the deciding exception, or ``unknown`` and ``None``.
    """
    for exc in _causes(error):
        for exc_class, failure_reason in _FAILURE_REASONS.items():
            if isinstance(exc, exc_class):
                return failure_reason, exc
    return FailureReasonEnum.UNKNOWN, None


class DbContext:
    """
    Context manager that brackets a long-running operation with DB state transitions.

    **Lifecycle:**

    1. *Enter*: connects to the DB, validates prerequisites (see below), then
       transitions the submission to ``start_state``.
    2. *Body*: the caller performs the actual work.
    3. *Exit (success)*: transitions the submission to ``end_state``.
       *Exit (exception)*: transitions the submission to ``ERROR`` and stores
       ``{"error": "<message>"}`` in the state log; the exception is then re-raised.

    **Prerequisite validation:**

    Before setting ``start_state``, the current state of the submission is compared
    against ``expected_prior_states``:

    - If the current state **matches** one of the expected prior states, the
      transition proceeds normally.
    - If the current state **does not match**, a warning is logged but the
      transition still proceeds (no hard failure).
    - If the submission **does not exist** in the DB:

      - for the entry states, ``PROCESSING`` and ``UPLOADING``, the submission is
        automatically created and the transition proceeds;
      - otherwise: ``SubmissionNotFoundError`` is raised immediately.

    Errors raised inside ``__enter__`` (other than ``SubmissionNotFoundError``) are
    wrapped in ``RuntimeError`` so callers always receive a consistent exception type.

    Example::

        with DbContext(
            config,
            submission_id,
            start_state=SubmissionStateEnum.DOWNLOADING,
            end_state=SubmissionStateEnum.DOWNLOADED,
        ):
            do_work_here()

    :param configuration: Nested dictionary that must contain a ``"db"`` key matching
        the ``DbConfig`` model (database URL, optional author credentials, …).
    :param submission_id: Submission ID to operate on.
    :param start_state: State written to the DB when entering the context.
    :param end_state: State written to the DB when exiting the context successfully.
    :param enabled: Set to ``False`` to skip all DB interactions (useful when no DB
        is configured).
    """

    _SUBMISSION_ENTRY_STATES = frozenset({SubmissionStateEnum.PROCESSING, SubmissionStateEnum.UPLOADING})
    """States at which a brand-new submission may be created.

    These are explicit because ``PROCESSING`` is not the enum member ``UPLOADING`` precedes.
    """

    def __init__(
        self,
        configuration: dict[str, Any] | GrzctlConfig,
        submission_id: str,
        start_state: SubmissionStateEnum,
        end_state: SubmissionStateEnum,
        enabled: bool = True,
    ):
        self.configuration = configuration
        self.submission_id = submission_id
        self.start_state = start_state
        self.end_state = end_state
        self.enabled = enabled
        self.db: SubmissionDb | None = None

    @cached_property
    def config(self) -> GrzctlConfig:
        """Parse and cache the GrzctlConfig from the raw configuration dict or object."""
        if isinstance(self.configuration, GrzctlConfig):
            return self.configuration
        return GrzctlConfig.from_configuration(self.configuration)

    @cached_property
    def grzctl_versions(self) -> dict[str, str]:
        """Get version information."""
        versions = get_versions()
        return {k: v or "unknown" for k, v in versions.items()}

    @cached_property
    def expected_prior_states(self) -> set[SubmissionStateEnum | None]:
        """Return the states the submission may be in before transitioning to ``start_state``.

        The entry states start a new submission, so they expect no prior state at all.
        Every other transition expects the previous ``SubmissionStateEnum`` member,
        whose order mirrors the pipeline order.
        """
        if self.start_state in self._SUBMISSION_ENTRY_STATES:
            return {None}
        members = list(SubmissionStateEnum)
        start_index = members.index(self.start_state)
        return {members[start_index - 1]}

    def __enter__(self):
        """Initializes DB connection, checks prerequisites, and sets the initial state."""
        if not self.enabled:
            return self

        try:
            db_config = self.config.db

            self.db = get_submission_db_instance(db_config.database_url, author=self.author)

            # Check if the state transition is valid.
            self._check_prerequisites()

            log.debug(f"Updating submission {self.submission_id} state to {self.start_state.name}")
            self.db.update_submission_state(self.submission_id, self.start_state, grzctl_versions=self.grzctl_versions)

        except SubmissionNotFoundError:
            raise
        except Exception as e:
            raise RuntimeError("Failed to connect to DB") from e

        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_val: BaseException | None,
        exc_tb: Any,
    ) -> bool:
        """Exit the database context.

        Commits the transaction if no exception occurred, otherwise rolls back.

        :returns: ``False``, so any exception is propagated.
        """
        if not self.db:
            return False

        if exc_type:
            error_state = SubmissionStateEnum.ERROR
            failure_reason, deciding = _classify(exc_val)
            recorded = deciding if deciding is not None else exc_val
            # an interruption carries no message, so its type names it
            data: dict[str, Any] = {"error": str(recorded) or type(recorded).__name__}
            log.error(f"Operation failed for {self.submission_id}. Updating DB to {error_state.name}.")
            try:
                self.db.update_submission_state(
                    self.submission_id,
                    error_state,
                    failure_reason=failure_reason,
                    data=data,
                    grzctl_versions=self.grzctl_versions,
                )
            except Exception as db_exc:
                log.error(f"Failed to write error state to DB: {db_exc}")
            return False

        else:
            log.info(f"Operation successful. Updating DB to {self.end_state.name}.")
            try:
                self.db.update_submission_state(
                    self.submission_id, self.end_state, grzctl_versions=self.grzctl_versions
                )
            except Exception as db_exc:
                log.error(f"Failed to write success state to DB: {db_exc}")

        return True

    @property
    def author(self) -> Author:
        return self.config.db.signing_author

    def _map_exception_to_failure_reason(
        self, exc_type: type[BaseException], exc_val: BaseException | None
    ) -> FailureReasonEnum:
        """Map an exception to its failure reason, see :func:`_classify`."""
        return _classify(exc_val)[0]

    def _check_prerequisites(self):
        """
        Checks if the state transition is valid.
        """
        submission = self.db.get_submission(self.submission_id)
        if submission is None:
            if None in self.expected_prior_states:
                # submission missing in DB is expected; add submission to DB
                submission = self.db.add_submission(self.submission_id)
            else:
                raise SubmissionNotFoundError(self.submission_id)

        latest_state_log = submission.get_latest_state()
        current_state = latest_state_log.state if latest_state_log else None

        if current_state not in self.expected_prior_states:
            log.warning(
                f"Submission {self.submission_id} is currently in state '{current_state}'. "
                f"Expected any of '{self.expected_prior_states}' before updating to '{self.start_state.name}'."
            )

        # The entry states expect no prior state ({None}), and no state-log entry has
        # state None, so the history check below would always warn for a new submission.
        if None in self.expected_prior_states:
            return

        history = submission.states
        found_in_history = any(entry.state in self.expected_prior_states for entry in history)

        if not found_in_history:
            log.warning(
                f"Submission {self.submission_id} is being updated to '{self.start_state.name}' "
                f"but state history does not contain any of '{self.expected_prior_states}'."
            )
