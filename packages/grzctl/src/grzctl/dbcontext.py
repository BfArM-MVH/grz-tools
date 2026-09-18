import logging
from functools import cached_property
from typing import Any

from grz_common.exceptions import (
    DecryptionError,
    DetailedQCError,
    EncryptionError,
    IncompleteSubmissionError,
    NetworkError,
    UploadError,
)
from grz_common.pipeline.components import DataValidationError
from grz_db.errors import DuplicateInitialSubmissionError, DuplicateTanGError, SubmissionNotFoundError
from grz_db.models.author import Author
from grz_db.models.submission import FailureReasonEnum, SubmissionDb, SubmissionStateEnum
from pydantic import ValidationError

from . import get_versions
from .commands.db.cli import get_submission_db_instance
from .models.config import GrzctlConfig

log = logging.getLogger(__name__)


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

      - for the entry states (``PROCESSING`` via ``grzctl process`` and
        ``UPLOADING`` via ``grzctl upload``) the submission is automatically
        created and the transition proceeds;
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

    ``grzctl process`` starts the streaming pipeline, ``grzctl upload`` the manual step-by-step
    flow. These are explicit because ``PROCESSING`` is not the enum member ``UPLOADING`` precedes.
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
            error_message = str(exc_val)
            error_state = SubmissionStateEnum.ERROR
            failure_reason = self._map_exception_to_failure_reason(exc_type, exc_val)  # new
            log.error(f"Operation failed for {self.submission_id}. Updating DB to {error_state.name}.")
            try:
                self.db.update_submission_state(
                    self.submission_id,
                    error_state,
                    failure_reason=failure_reason,
                    data={"error": error_message},
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
        # cached on the configuration, so the several contexts of one run share a single unlocked key
        return self.config.db.signing_author

    def _map_exception_to_failure_reason(
        self, exc_type: type[BaseException], exc_val: BaseException | None
    ) -> FailureReasonEnum:
        """Maps an exception to the closest FailureReasonEnum value.

        Checks ``exc_val``, then its ``__cause__``, then the cause's ``__cause__``, and so on.
        The first exception whose type is mapped decides the failure reason.
        """
        exception_map: dict[type[BaseException], FailureReasonEnum] = {
            FileNotFoundError: FailureReasonEnum.FILE_NOT_FOUND,
            ValidationError: FailureReasonEnum.VALIDATION_ERROR,
            DataValidationError: FailureReasonEnum.VALIDATION_ERROR,
            DecryptionError: FailureReasonEnum.DECRYPTION_ERROR,
            EncryptionError: FailureReasonEnum.ENCRYPTION_ERROR,
            NetworkError: FailureReasonEnum.NETWORK_ERROR,
            UploadError: FailureReasonEnum.UPLOAD_ERROR,
            DetailedQCError: FailureReasonEnum.DETAILED_QC_ERROR,
            DuplicateTanGError: FailureReasonEnum.DUPLICATE_TANG,
            DuplicateInitialSubmissionError: FailureReasonEnum.DUPLICATE_INITIAL,
            IncompleteSubmissionError: FailureReasonEnum.INCOMPLETE_SUBMISSION,
        }
        seen: set[int] = set()
        exc = exc_val
        # a cause can form a cycle, as in ``raise e from e``
        while exc is not None and id(exc) not in seen:
            for exc_class, failure_reason in exception_map.items():
                if isinstance(exc, exc_class):
                    return failure_reason
            seen.add(id(exc))
            exc = exc.__cause__
        return FailureReasonEnum.UNKNOWN

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
