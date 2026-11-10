import logging
import subprocess
from unittest import mock
from unittest.mock import MagicMock

import grz_common.exceptions as grzexc
import pytest
from grz_common.interrupt import TerminateInterrupt
from grz_common.pipeline.components import DataValidationError, UploadIntegrityError
from grz_common.pipeline.context import FileError
from grz_common.workers.submission import SubmissionMetadata
from grz_db.errors import DuplicateInitialSubmissionError, DuplicateTanGError, SubmissionNotFoundError
from grz_db.models.submission import RETIRED_FAILURE_REASONS, FailureReasonEnum, SubmissionStateEnum
from grzctl.dbcontext import PROCESS_TRANSITIONS, STEP_TRANSITIONS, DbContext, FilesFailedError, Transitions
from pydantic import ValidationError


@pytest.fixture
def db_context() -> DbContext:
    """Create a DbContext instance without triggering __init__ DB connections."""
    ctx = DbContext.__new__(DbContext)
    ctx.enabled = False
    ctx.db = None
    return ctx


@pytest.fixture
def mock_db():
    return MagicMock()


@pytest.fixture
def ctx(mock_db):
    context = DbContext(
        configuration={},
        submission_id="123_2025-01-01_00000000",
        start_state=SubmissionStateEnum.ENCRYPTING,
        end_state=SubmissionStateEnum.ENCRYPTED,
        transitions=STEP_TRANSITIONS,
        enabled=True,
    )
    context.db = mock_db  # bypass __enter__
    return context


class TestMapExceptionToFailureReason:
    @pytest.mark.parametrize(
        "exception,expected",
        [
            (FileNotFoundError("missing file"), FailureReasonEnum.UNKNOWN),
            (grzexc.DecryptionError("failed"), FailureReasonEnum.DECRYPTION_ERROR),
            (grzexc.EncryptionError("failed"), FailureReasonEnum.ENCRYPTION_ERROR),
            (grzexc.MissingSubmissionFileError("failed"), FailureReasonEnum.FILE_NOT_FOUND),
            (grzexc.DownloadError("failed"), FailureReasonEnum.TRANSFER_ERROR),
            (grzexc.MissingObjectError("failed"), FailureReasonEnum.TRANSFER_ERROR),
            (grzexc.DuplicateUploadError("failed"), FailureReasonEnum.DUPLICATE_TANG),
            (grzexc.NetworkError("failed"), FailureReasonEnum.TRANSFER_ERROR),
            (grzexc.UploadError("failed"), FailureReasonEnum.TRANSFER_ERROR),
            (UploadIntegrityError("failed"), FailureReasonEnum.TRANSFER_ERROR),
            (grzexc.ConfigurationError("failed"), FailureReasonEnum.CONFIGURATION_ERROR),
            (grzexc.PruefberichtGenerationError("failed"), FailureReasonEnum.PRUEFBERICHT_GENERATION_ERROR),
            (grzexc.PruefberichtRejectedError("failed"), FailureReasonEnum.PRUEFBERICHT_REJECTED),
            (KeyboardInterrupt(), FailureReasonEnum.INTERRUPTED),
            (TerminateInterrupt(), FailureReasonEnum.INTERRUPTED),
            (DataValidationError("failed"), FailureReasonEnum.VALIDATION_ERROR),
            (DuplicateTanGError(), FailureReasonEnum.DUPLICATE_TANG),
            (grzexc.IncompleteSubmissionError("failed"), FailureReasonEnum.INCOMPLETE_SUBMISSION),
            (grzexc.SubmissionCleanedError("failed"), FailureReasonEnum.SUBMISSION_CLEANED),
            (grzexc.DetailedQCError("failed"), FailureReasonEnum.DETAILED_QC_ERROR),
            (subprocess.CalledProcessError(returncode=3, cmd="some other command"), FailureReasonEnum.UNKNOWN),
            (RuntimeError("unexpected"), FailureReasonEnum.UNKNOWN),
            (Exception("generic"), FailureReasonEnum.UNKNOWN),
            (ValueError("some value error"), FailureReasonEnum.UNKNOWN),
        ],
    )
    def test_maps_correctly(self, db_context: DbContext, exception: BaseException, expected: FailureReasonEnum):
        result = db_context._map_exception_to_failure_reason(type(exception), exception)
        assert result == expected

    def test_an_unwrapped_pydantic_error_maps_to_unknown(self, db_context: DbContext):
        """Invalid metadata reaches DbContext as a SubmissionValidationError, so a bare pydantic error is a bug."""
        exc = ValidationError.from_exception_data("test", [])
        assert db_context._map_exception_to_failure_reason(type(exc), exc) == FailureReasonEnum.UNKNOWN

    def test_invalid_metadata_maps_to_a_validation_error(self, db_context: DbContext, tmp_path):
        metadata_file = tmp_path / "metadata.json"
        metadata_file.write_text('{"submission": {}}')

        with pytest.raises(grzexc.SubmissionValidationError) as excinfo:
            SubmissionMetadata(metadata_file)

        exc = excinfo.value
        assert db_context._map_exception_to_failure_reason(type(exc), exc) == FailureReasonEnum.VALIDATION_ERROR

    def test_maps_a_cause_of_an_unmapped_exception(self, db_context: DbContext):
        """An exception raised ``from`` a mapped one gets the failure reason of its cause."""
        with pytest.raises(RuntimeError) as exc_info:
            raise RuntimeError("processing failed") from grzexc.UploadError("upload failed")
        result = db_context._map_exception_to_failure_reason(exc_info.type, exc_info.value)
        assert result == FailureReasonEnum.TRANSFER_ERROR

    def test_none_exception_returns_unknown(self, db_context: DbContext):
        result = db_context._map_exception_to_failure_reason(type(None), None)
        assert result == FailureReasonEnum.UNKNOWN

    def test_all_enum_values_are_covered(self, db_context: DbContext):
        """Ensures every FailureReasonEnum value except UNKNOWN is reachable via a mapped exception."""
        mapped_results = {
            db_context._map_exception_to_failure_reason(type(exc), exc)
            for exc in [
                grzexc.MissingSubmissionFileError(),
                grzexc.DecryptionError(),
                grzexc.EncryptionError(),
                grzexc.TransferError(),
                grzexc.ConfigurationError(),
                grzexc.PruefberichtGenerationError(),
                grzexc.PruefberichtRejectedError(),
                KeyboardInterrupt(),
                DuplicateTanGError(),
                DuplicateInitialSubmissionError(1),
                grzexc.IncompleteSubmissionError(),
                grzexc.SubmissionCleanedError(),
                grzexc.DetailedQCError(),
                grzexc.SubmissionValidationError(),
            ]
        }
        recorded = {e for e in FailureReasonEnum if e != FailureReasonEnum.UNKNOWN} - RETIRED_FAILURE_REASONS
        unmapped = recorded - mapped_results
        assert not unmapped, f"These FailureReasonEnum values have no exception mapping: {unmapped}"

    def test_every_expected_failure_records_a_current_reason(self, db_context: DbContext):
        """A GrzError that maps to ``unknown`` would be recorded as a bug, and a retired reason not at all.

        This walks every subclass, not only the leaves: a class with subclasses of its own,
        such as ``TransferError``, can still be raised directly and needs its own reason.
        """

        def all_subclasses(cls: type[grzexc.GrzError]) -> list[type[grzexc.GrzError]]:
            direct = cls.__subclasses__()
            return direct + [sub for subclass in direct for sub in all_subclasses(subclass)]

        # grouping classes that no code raises directly, so they carry no failure reason of their own
        grouping_classes = {grzexc.GrzError, grzexc.SubmissionRejectedError}
        # takes file errors, and its reason is that of its decisive file error (see ``TestFilesFailed``)
        grouping_classes.add(FilesFailedError)

        reasons = {
            cls.__name__: db_context._map_exception_to_failure_reason(cls, cls("failed"))
            for cls in all_subclasses(grzexc.GrzError)
            if cls not in grouping_classes
        }

        wrong = {
            name: reason
            for name, reason in reasons.items()
            if reason in {FailureReasonEnum.UNKNOWN, *RETIRED_FAILURE_REASONS}
        }
        assert not wrong, f"These errors record no current failure reason: {wrong}"


def _files_failed(*errors: BaseException) -> FilesFailedError:
    return FilesFailedError("Processing failed", [FileError(f"file{i}", error) for i, error in enumerate(errors)])


class TestFilesFailed:
    """Several files fail in one run, and the most decisive of their errors sets the failure reason."""

    @pytest.mark.parametrize(
        "errors,expected",
        [
            (
                (grzexc.TransferError("slow"), grzexc.SubmissionValidationError("bad")),
                FailureReasonEnum.VALIDATION_ERROR,
            ),
            (
                (grzexc.TransferError("slow"), grzexc.ConfigurationError("refused")),
                FailureReasonEnum.CONFIGURATION_ERROR,
            ),
            ((grzexc.EncryptionError("failed"), grzexc.TransferError("slow")), FailureReasonEnum.TRANSFER_ERROR),
            (
                (grzexc.DecryptionError("failed"), grzexc.SubmissionValidationError("bad")),
                FailureReasonEnum.DECRYPTION_ERROR,
            ),
            ((RuntimeError("bug"), grzexc.EncryptionError("failed")), FailureReasonEnum.UNKNOWN),
        ],
        ids=[
            "submitter-before-transfer",
            "configuration-before-transfer",
            "transfer-before-other",
            "first-of-equal-rank",
            "first-of-other",
        ],
    )
    def test_the_most_decisive_error_sets_the_reason(
        self, db_context: DbContext, errors: tuple[BaseException, ...], expected: FailureReasonEnum
    ):
        exc = _files_failed(*errors)
        assert db_context._map_exception_to_failure_reason(type(exc), exc) == expected

    def test_is_an_expected_failure(self):
        """``grzctl.cli`` reports a ``GrzError`` and shows a traceback for any other exception."""
        assert issubclass(FilesFailedError, grzexc.GrzError)

    def test_every_file_error_is_recorded(self, ctx, mock_db):
        """The decisive error is the recorded one, and every file error is listed with its own reason."""
        exc = _files_failed(grzexc.TransferError("slow"), grzexc.SubmissionValidationError("bad"))

        ctx.__exit__(type(exc), exc, None)

        mock_db.update_submission_state.assert_called_once_with(
            ctx.submission_id,
            SubmissionStateEnum.ERROR,
            data={
                "error": "bad",
                "errors": [
                    {"file": "file0", "reason": "transfer_error", "message": "slow"},
                    {"file": "file1", "reason": "validation_error", "message": "bad"},
                ],
            },
            failure_reason=FailureReasonEnum.VALIDATION_ERROR,
            grzctl_versions=mock.ANY,
        )

    def test_an_unmapped_decisive_error_is_the_recorded_one(self, ctx, mock_db):
        """With no mapped reason, ``data.error`` is still the message of the decisive file error."""
        exc = _files_failed(RuntimeError("unexpected"))

        ctx.__exit__(type(exc), exc, None)

        mock_db.update_submission_state.assert_called_once_with(
            ctx.submission_id,
            SubmissionStateEnum.ERROR,
            data={
                "error": "unexpected",
                "errors": [{"file": "file0", "reason": "unknown", "message": "unexpected"}],
            },
            failure_reason=FailureReasonEnum.UNKNOWN,
            grzctl_versions=mock.ANY,
        )

    def test_a_file_error_without_a_message_is_listed_by_its_type(self, ctx, mock_db):
        """A bare ``MemoryError()`` has an empty ``str()``, so its entry names the type like ``data.error`` does."""
        exc = _files_failed(MemoryError())

        ctx.__exit__(type(exc), exc, None)

        data = mock_db.update_submission_state.call_args.kwargs["data"]
        assert data["errors"] == [{"file": "file0", "reason": "unknown", "message": "MemoryError"}]


class TestDbContextFailureReason:
    def test_file_not_found_maps_correctly(self, ctx, mock_db):
        exc = grzexc.MissingSubmissionFileError("missing file")
        ctx.__exit__(type(exc), exc, None)
        mock_db.update_submission_state.assert_called_once_with(
            ctx.submission_id,
            SubmissionStateEnum.ERROR,
            data={"error": str(exc)},
            failure_reason=FailureReasonEnum.FILE_NOT_FOUND,
            grzctl_versions=mock.ANY,
        )

    def test_an_interruption_is_recorded_by_its_type(self, ctx, mock_db):
        """A KeyboardInterrupt carries no message, so the recorded error names its type."""
        exc = KeyboardInterrupt()
        ctx.__exit__(type(exc), exc, None)
        mock_db.update_submission_state.assert_called_once_with(
            ctx.submission_id,
            SubmissionStateEnum.ERROR,
            data={"error": "KeyboardInterrupt"},
            failure_reason=FailureReasonEnum.INTERRUPTED,
            grzctl_versions=mock.ANY,
        )

    def test_unknown_exception_maps_to_unknown(self, ctx, mock_db):
        exc = RuntimeError("something unexpected")
        ctx.__exit__(type(exc), exc, None)
        mock_db.update_submission_state.assert_called_once_with(
            ctx.submission_id,
            SubmissionStateEnum.ERROR,
            data={"error": str(exc)},
            failure_reason=FailureReasonEnum.UNKNOWN,
            grzctl_versions=mock.ANY,
        )

    def test_no_exception_does_not_set_failure_reason(self, ctx, mock_db):
        ctx.__exit__(None, None, None)
        mock_db.update_submission_state.assert_called_once_with(
            ctx.submission_id,
            SubmissionStateEnum.ENCRYPTED,
            grzctl_versions=mock.ANY,
        )


class TestCheckPrerequisites:
    @staticmethod
    def _context(
        mock_db,
        start_state: SubmissionStateEnum,
        end_state: SubmissionStateEnum,
        transitions: Transitions = STEP_TRANSITIONS,
    ) -> DbContext:
        context = DbContext(
            configuration={},
            submission_id="123_2025-01-01_00000000",
            start_state=start_state,
            end_state=end_state,
            transitions=transitions,
            enabled=True,
        )
        context.db = mock_db  # bypass __enter__
        return context

    def test_first_state_does_not_warn_about_the_history(self, mock_db, caplog):
        """A brand-new submission has no prior state to find, so nothing is logged."""
        mock_db.get_submission.return_value = None
        mock_db.add_submission.return_value.get_latest_state.return_value = None
        mock_db.add_submission.return_value.states = []
        context = self._context(
            mock_db, SubmissionStateEnum.PROCESSING, SubmissionStateEnum.UPLOADING, PROCESS_TRANSITIONS
        )

        with caplog.at_level(logging.WARNING):
            context._check_prerequisites()

        assert caplog.records == []
        mock_db.add_submission.assert_called_once()

    def test_manual_upload_entry_starts_a_brand_new_submission(self, mock_db, caplog):
        """The step-by-step flow starts at ``UPLOADING``, which is an entry state too and needs no prior state."""
        mock_db.get_submission.return_value = None
        mock_db.add_submission.return_value.get_latest_state.return_value = None
        mock_db.add_submission.return_value.states = []
        context = self._context(mock_db, SubmissionStateEnum.UPLOADING, SubmissionStateEnum.UPLOADED)

        with caplog.at_level(logging.WARNING):
            context._check_prerequisites()

        assert caplog.records == []
        mock_db.add_submission.assert_called_once()

    def test_middle_state_raises_when_the_submission_does_not_exist(self, mock_db):
        """A manual step needs the submission to exist; a missing one is a hard error."""
        mock_db.get_submission.return_value = None
        context = self._context(mock_db, SubmissionStateEnum.DOWNLOADING, SubmissionStateEnum.DOWNLOADED)

        with pytest.raises(SubmissionNotFoundError):
            context._check_prerequisites()

        mock_db.add_submission.assert_not_called()

    @pytest.mark.parametrize(
        ("start_state", "end_state", "transitions"),
        [
            (SubmissionStateEnum.PROCESSING, SubmissionStateEnum.PROCESSED, PROCESS_TRANSITIONS),
            (SubmissionStateEnum.UPLOADING, SubmissionStateEnum.UPLOADED, STEP_TRANSITIONS),
        ],
    )
    def test_entry_state_adds_a_new_submission(self, db, start_state, end_state, transitions):
        """Both entry states start a submission that the DB does not know yet.

        A real DB, unlike a mock, fails if the new submission's states are not loaded.
        """
        context = DbContext(
            configuration={},
            submission_id="123456789_2025-01-01_a1b2c3d4",
            start_state=start_state,
            end_state=end_state,
            transitions=transitions,
            enabled=True,
        )
        context.db = db  # bypass __enter__

        context._check_prerequisites()

        assert db.get_submission(context.submission_id) is not None

    def test_later_state_warns_when_the_history_lacks_the_prior_state(self, mock_db, caplog):
        """The history check still runs for every state that has a prior state."""
        mock_db.get_submission.return_value.get_latest_state.return_value = None
        mock_db.get_submission.return_value.states = []
        context = self._context(mock_db, SubmissionStateEnum.ENCRYPTING, SubmissionStateEnum.ENCRYPTED)

        with caplog.at_level(logging.WARNING):
            context._check_prerequisites()

        assert "state history does not contain" in caplog.text


S = SubmissionStateEnum

PROCESS_WITHOUT_QC = [S.PROCESSING, S.PROCESSED, S.REPORTING, S.REPORTED, S.CLEANING, S.CLEANED]
STEP_BY_STEP = [
    *(S.UPLOADING, S.UPLOADED, S.DOWNLOADING, S.DOWNLOADED, S.DECRYPTING, S.DECRYPTED),
    *(S.VALIDATING, S.VALIDATED, S.ENCRYPTING, S.ENCRYPTED, S.ARCHIVING, S.ARCHIVED),
    *(S.REPORTING, S.REPORTED, S.CLEANING, S.CLEANED),
]


class TestTransitions:
    """Each flow expects the states of its own steps, so a successful run of either logs no warning."""

    SUBMISSION_ID = "123456789_2025-01-01_a1b2c3d4"

    def _record(self, db, transitions: Transitions, states: list[SubmissionStateEnum]) -> None:
        """Record ``states`` in order, and check the prerequisites of each start state of the flow first.

        A state that is no start state of the flow comes from elsewhere, such as ``grzctl db sync``.
        """
        db.add_submission(self.SUBMISSION_ID)
        for state in states:
            if state in transitions:
                context = DbContext(
                    configuration={},
                    submission_id=self.SUBMISSION_ID,
                    start_state=state,
                    end_state=state,
                    transitions=transitions,
                )
                context.db = db  # bypass __enter__
                context._check_prerequisites()
            db.update_submission_state(self.SUBMISSION_ID, state)

    @pytest.mark.parametrize(
        ("transitions", "states"),
        [
            (PROCESS_TRANSITIONS, PROCESS_WITHOUT_QC),
            (PROCESS_TRANSITIONS, [*PROCESS_WITHOUT_QC[:4], S.QCING, S.QCED, S.CLEANING, S.CLEANED]),
            (PROCESS_TRANSITIONS, [S.PROCESSING, S.PROCESSED, S.QCING, S.QCED]),
            (PROCESS_TRANSITIONS, [S.UPLOADING, S.UPLOADED, *PROCESS_WITHOUT_QC]),
            (PROCESS_TRANSITIONS, [*PROCESS_WITHOUT_QC[:4], *PROCESS_WITHOUT_QC[:2], S.CLEANING, S.CLEANED]),
            (PROCESS_TRANSITIONS, [*PROCESS_WITHOUT_QC[:3], S.ERROR, *PROCESS_WITHOUT_QC]),
            (STEP_TRANSITIONS, STEP_BY_STEP),
            (STEP_TRANSITIONS, [*STEP_BY_STEP[:-2], S.QCING, S.QCED, S.CLEANING, S.CLEANED]),
            (STEP_TRANSITIONS, [*PROCESS_WITHOUT_QC[:2], S.REPORTING, S.REPORTED, S.CLEANING, S.CLEANED]),
            (STEP_TRANSITIONS, [*STEP_BY_STEP[:13], S.ERROR, *STEP_BY_STEP[12:]]),
        ],
        ids=[
            "process",
            "process-with-detailed-qc",
            "process-without-pruefbericht",
            "process-after-db-sync",
            "process-rerun-after-no-clean-inbox",
            "process-rerun-after-a-failed-pruefbericht",
            "step-by-step",
            "step-by-step-with-detailed-qc",
            "pruefbericht-submit-and-clean-after-process",
            "step-by-step-retry-of-a-failed-pruefbericht",
        ],
    )
    def test_a_successful_run_logs_no_warning(self, db, caplog, transitions, states):
        with caplog.at_level(logging.WARNING, logger="grzctl.dbcontext"):
            self._record(db, transitions, states)

        assert caplog.records == []

    @pytest.mark.parametrize(
        ("transitions", "states"),
        [
            (STEP_TRANSITIONS, [*STEP_BY_STEP[:6], S.ENCRYPTING]),
            (STEP_TRANSITIONS, [*STEP_BY_STEP[:12], S.CLEANING]),
            (STEP_TRANSITIONS, [*STEP_BY_STEP[:13], S.ERROR, S.CLEANING]),
            (PROCESS_TRANSITIONS, [S.PROCESSING, S.REPORTING]),
        ],
        ids=[
            "encrypt-before-validate",
            "clean-before-pruefbericht",
            "pruefbericht-before-processed",
            "another-step-after-an-error",
        ],
    )
    def test_an_out_of_order_step_warns(self, db, caplog, transitions, states):
        with caplog.at_level(logging.WARNING, logger="grzctl.dbcontext"):
            self._record(db, transitions, states)

        assert "currently in state" in caplog.text
        assert "state history does not contain" in caplog.text
