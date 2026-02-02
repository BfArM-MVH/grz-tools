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
        """Initializes DB connection and sets the initial state."""
        if not self.enabled:
            return self

        try:
            db_config = DbConfig.model_validate(self.configuration).db
            self.db = get_submission_db_instance(db_config.db_url, author=db_config.author)

            if self.db:
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
