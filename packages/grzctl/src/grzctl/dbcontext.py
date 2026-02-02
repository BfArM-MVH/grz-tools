import logging
from typing import Any

from grz_db.models.submission import SubmissionDb, SubmissionStateEnum
from pydantic import ValidationError

from .commands.db.cli import get_submission_db_instance
from .models.config import DbConfig

log = logging.getLogger(__name__)


class DbContext:
    """
    Context manager to handle automatic DB state updates for state updates.

    Usage:
        with DbContext(config, submission_id, start_state=SubmissionStateEnum.DOWNLOADING, end_state=SubmissionStateEnum.DOWNLOADED) as db:
            do_work_here()

    If an exception occurs within the block, the state is updated to ERROR,
    and the error message is added in the `data` blob, i.e., `{"error": str(error)}`.
    If the block finishes successfully, the state is updated to success_state.
    """

    def __init__(
        self,
        configuration: dict[str, Any],
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

    def __enter__(self):
        """Initializes DB connection, checks prerequisites, and sets the initial state."""
        if not self.enabled:
            return self

        try:
            db_config = DbConfig.model_validate(self.configuration).db
            self.db = get_submission_db_instance(db_config.db_url, author=db_config.author)

            if self.db:
                self._check_prerequisites()

                log.debug(f"Updating submission {self.submission_id} state to {self.start_state.name}")
                self.db.update_submission_state(self.submission_id, self.start_state)

        except (ValidationError, KeyError) as e:
            log.warning(f"DB Configuration invalid or missing. State updates skipped. ({e})")
            self.db = None
        except Exception as e:
            log.error(f"Failed to connect to DB: {e}. State updates skipped.")
            self.db = None

        return self

    def __exit__(self, exc_type: type[BaseException] | None, exc_val: BaseException | None, exc_tb: Any) -> bool:
        """Handles success or failure state updates."""
        if not self.db:
            return False

        if exc_type:
            error_message = str(exc_val)
            error_state = SubmissionStateEnum.ERROR
            log.error(f"Operation failed for {self.submission_id}. Updating DB to {error_state.name}.")
            try:
                self.db.update_submission_state(self.submission_id, error_state, data={"error": error_message})
            except Exception as db_exc:
                log.error(f"Failed to write error state to DB: {db_exc}")

            return False

        else:
            log.info(f"Operation successful. Updating DB to {self.end_state.name}.")
            try:
                self.db.update_submission_state(self.submission_id, self.end_state)
            except Exception as db_exc:
                log.error(f"Failed to write success state to DB: {db_exc}")

        return True

    def _check_prerequisites(self):
        """
        Checks if the state transition is valid.
        """
        if not self.db:
            return

        # determine expected prior state
        members = list(SubmissionStateEnum)
        start_index = members.index(self.start_state)
        if start_index == 0:
            # first state in the enum, no prior state expected
            return

        expected_prior_state = members[start_index - 1]
        expected_prior_state = str(expected_prior_state).casefold()

        submission = self.db.get_submission(self.submission_id)
        if not submission:
            return

        latest_state_log = submission.get_latest_state()
        current_state = latest_state_log.state if latest_state_log else None
        current_state = str(current_state).casefold() if current_state else None

        if current_state != expected_prior_state:
            log.warning(
                f"Submission {self.submission_id} is currently in state '{current_state}'. "
                f"Expected '{expected_prior_state}' before updating to '{self.start_state.name}'."
            )

        history = submission.states
        found_in_history = any(str(entry.state).casefold() == expected_prior_state for entry in history)

        if not found_in_history:
            log.warning(
                f"Submission {self.submission_id} is being updated to '{self.start_state.name}' "
                f"but state history does not contain '{expected_prior_state}'."
            )
