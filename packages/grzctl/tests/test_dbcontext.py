import logging
import subprocess
from unittest import mock
from unittest.mock import MagicMock

import grz_common.exceptions as grzexc
import pytest
from grz_common.workers.submission import SubmissionMetadata
from grz_db.errors import DuplicateInitialSubmissionError, DuplicateTanGError, SubmissionNotFoundError
from grz_db.models.submission import RETIRED_FAILURE_REASONS, FailureReasonEnum, SubmissionStateEnum
from grzctl.dbcontext import DbContext
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
            (grzexc.ConfigurationError("failed"), FailureReasonEnum.CONFIGURATION_ERROR),
            (grzexc.ReportingError("failed"), FailureReasonEnum.REPORTING_ERROR),
            (KeyboardInterrupt(), FailureReasonEnum.INTERRUPTED),
            (DuplicateTanGError(), FailureReasonEnum.DUPLICATE_TANG),
            (grzexc.IncompleteSubmissionError("failed"), FailureReasonEnum.INCOMPLETE_SUBMISSION),
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
                grzexc.ReportingError(),
                KeyboardInterrupt(),
                DuplicateTanGError(),
                DuplicateInitialSubmissionError(1),
                grzexc.IncompleteSubmissionError(),
                grzexc.DetailedQCError(),
                grzexc.SubmissionValidationError(),
            ]
        }
        recorded = {e for e in FailureReasonEnum if e != FailureReasonEnum.UNKNOWN} - RETIRED_FAILURE_REASONS
        unmapped = recorded - mapped_results
        assert not unmapped, f"These FailureReasonEnum values have no exception mapping: {unmapped}"

    def test_every_expected_failure_records_a_current_reason(self, db_context: DbContext):
        """A GrzError that maps to ``unknown`` would be recorded as a bug, and a retired reason not at all."""

        def leaves(cls: type[grzexc.GrzError]) -> list[type[grzexc.GrzError]]:
            subclasses = cls.__subclasses__()
            return [leaf for subclass in subclasses for leaf in leaves(subclass)] if subclasses else [cls]

        reasons = {
            cls.__name__: db_context._map_exception_to_failure_reason(cls, cls("failed"))
            for cls in leaves(grzexc.GrzError)
        }

        wrong = {
            name: reason
            for name, reason in reasons.items()
            if reason in {FailureReasonEnum.UNKNOWN, *RETIRED_FAILURE_REASONS}
        }
        assert not wrong, f"These errors record no current failure reason: {wrong}"


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
    def _context(mock_db, start_state: SubmissionStateEnum, end_state: SubmissionStateEnum) -> DbContext:
        context = DbContext(
            configuration={},
            submission_id="123_2025-01-01_00000000",
            start_state=start_state,
            end_state=end_state,
            enabled=True,
        )
        context.db = mock_db  # bypass __enter__
        return context

    def test_first_state_does_not_warn_about_the_history(self, mock_db, caplog):
        """A brand-new submission has no prior state to find, so nothing is logged."""
        mock_db.get_submission.return_value = None
        mock_db.add_submission.return_value.get_latest_state.return_value = None
        mock_db.add_submission.return_value.states = []
        context = self._context(mock_db, SubmissionStateEnum.PROCESSING, SubmissionStateEnum.UPLOADING)

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
        ("start_state", "end_state"),
        [
            (SubmissionStateEnum.PROCESSING, SubmissionStateEnum.PROCESSED),
            (SubmissionStateEnum.UPLOADING, SubmissionStateEnum.UPLOADED),
        ],
    )
    def test_entry_state_adds_a_new_submission(self, db, start_state, end_state):
        """``grzctl process`` and the step-by-step flow both start a submission that the DB does not know yet.

        A real DB, unlike a mock, fails if the new submission's states are not loaded.
        """
        context = DbContext(
            configuration={},
            submission_id="123456789_2025-01-01_a1b2c3d4",
            start_state=start_state,
            end_state=end_state,
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
