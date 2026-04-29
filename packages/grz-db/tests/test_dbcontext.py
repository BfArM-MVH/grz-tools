from unittest.mock import MagicMock

import pytest
from grz_db.models.submission import FailureReasonEnum, SubmissionStateEnum
from grzctl.dbcontext import DbContext
from pydantic import ValidationError


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


class TestDbContextFailureReason:
    def test_file_not_found_maps_correctly(self, ctx, mock_db):
        exc = FileNotFoundError("missing file")
        ctx.__exit__(type(exc), exc, None)
        mock_db.update_submission_state.assert_called_once_with(
            ctx.submission_id,
            SubmissionStateEnum.ERROR,
            data={"error": str(exc)},
            failure_reason=FailureReasonEnum.FILE_NOT_FOUND,
        )

    def test_validation_error_maps_correctly(self, ctx, mock_db):
        exc = ValidationError.from_exception_data("test", [])
        ctx.__exit__(type(exc), exc, None)
        mock_db.update_submission_state.assert_called_once_with(
            ctx.submission_id,
            SubmissionStateEnum.ERROR,
            data={"error": str(exc)},
            failure_reason=FailureReasonEnum.VALIDATION_ERROR,
        )

    def test_unknown_exception_maps_to_unknown(self, ctx, mock_db):
        exc = RuntimeError("something unexpected")
        ctx.__exit__(type(exc), exc, None)
        mock_db.update_submission_state.assert_called_once_with(
            ctx.submission_id,
            SubmissionStateEnum.ERROR,
            data={"error": str(exc)},
            failure_reason=FailureReasonEnum.UNKNOWN,
        )

    def test_no_exception_does_not_set_failure_reason(self, ctx, mock_db):
        ctx.__exit__(None, None, None)
        mock_db.update_submission_state.assert_called_once_with(
            ctx.submission_id,
            SubmissionStateEnum.ENCRYPTED,
        )
