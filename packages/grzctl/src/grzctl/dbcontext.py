import logging
from collections.abc import Iterable
from functools import cached_property
from pathlib import Path
from typing import Any

from grz_db.errors import SubmissionNotFoundError
from grz_db.models.author import Author
from grz_db.models.submission import SubmissionDb, SubmissionStateEnum
from pydantic import ValidationError

from .commands.db.cli import get_submission_db_instance
from .models.config import DbConfig

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

      - and ``None`` is in ``expected_prior_states``: the submission is
        automatically created and the transition proceeds.
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
    :param expected_prior_states: Iterable of states considered valid preconditions.
        Pass ``None`` (the default) to derive the expected prior state automatically
          from the state preceding ``start_state`` in ``SubmissionStateEnum``.
        Include ``None`` as an *element* to allow the submission to be absent from
          the DB; it will then be created on the fly.
    :param enabled: Set to ``False`` to skip all DB interactions (useful when no DB
        is configured).
    """

    def __init__(  # noqa: PLR0913
        self,
        configuration: dict[str, Any],
        submission_id: str,
        start_state: SubmissionStateEnum,
        end_state: SubmissionStateEnum,
        expected_prior_states: Iterable[SubmissionStateEnum | None] | None = None,
        enabled: bool = True,
    ):
        self.configuration = configuration
        self.submission_id = submission_id
        self.start_state = start_state
        self.end_state = end_state
        self.enabled = enabled
        self._expected_prior_states = set(expected_prior_states) if expected_prior_states else None
        self.db: SubmissionDb | None = None

    @property
    def expected_prior_states(self) -> set[SubmissionStateEnum | None]:
        if self._expected_prior_states is None:
            # determine expected prior state based on order of enums
            members = list(SubmissionStateEnum)
            start_index = members.index(self.start_state)

            if start_index == 0:
                # first state in the enum, no prior state expected
                return {None}
            else:
                # return previous state in the enum as expected prior state
                return {members[start_index - 1]}
        else:
            return self._expected_prior_states

    def __enter__(self):
        """Initializes DB connection, checks prerequisites, and sets the initial state."""
        if not self.enabled:
            return self

        try:
            db_config = DbConfig.model_validate(self.configuration).db

            self.db = get_submission_db_instance(db_config.database_url, author=self.author)

            # Check if the state transition is valid.
            self._check_prerequisites()

            log.debug(f"Updating submission {self.submission_id} state to {self.start_state.name}")
            self.db.update_submission_state(self.submission_id, self.start_state)

        except SubmissionNotFoundError:
            raise
        except (ValidationError, KeyError) as e:
            raise RuntimeError("DB Configuration invalid or missing") from e
        except Exception as e:
            raise RuntimeError("Failed to connect to DB") from e

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

    @cached_property
    def author(self) -> Author:
        db_config = DbConfig.model_validate(self.configuration).db

        if not db_config.author:
            raise ValueError("Author configuration is missing")

        if db_config.author.private_key_path is None:
            raise ValueError("Author private key path is required but was None")

        key_path = Path(db_config.author.private_key_path)
        if not key_path.exists():
            raise FileNotFoundError(f"Author private key not found at: {key_path}")

        return Author(
            name=db_config.author.name,
            private_key_bytes=key_path.read_bytes(),
            private_key_passphrase=db_config.author.private_key_passphrase,
        )

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
                f"Expected any of '{self._expected_prior_states}' before updating to '{self.start_state.name}'."
            )

        history = submission.states
        found_in_history = any(entry.state in self.expected_prior_states for entry in history)

        if not found_in_history:
            log.warning(
                f"Submission {self.submission_id} is being updated to '{self.start_state.name}' "
                f"but state history does not contain any of '{self.expected_prior_states}'."
            )
