from unittest import mock
from unittest.mock import MagicMock

import pytest
from grz_common.exceptions import (
    DecryptionError,
    EncryptionError,
    IncompleteSubmissionError,
    NetworkError,
    UploadError,
)
from grz_db.errors import DuplicateTanGError
from grz_db.models.submission import FailureReasonEnum, SubmissionStateEnum
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
            (FileNotFoundError("missing file"), FailureReasonEnum.FILE_NOT_FOUND),
            (DecryptionError("failed"), FailureReasonEnum.DECRYPTION_ERROR),
            (EncryptionError("failed"), FailureReasonEnum.ENCRYPTION_ERROR),
            (NetworkError("failed"), FailureReasonEnum.NETWORK_ERROR),
            (UploadError("failed"), FailureReasonEnum.UPLOAD_ERROR),
            (DuplicateTanGError(), FailureReasonEnum.DUPLICATE_TANG),
            (IncompleteSubmissionError("failed"), FailureReasonEnum.INCOMPLETE_SUBMISSION),
            (RuntimeError("unexpected"), FailureReasonEnum.UNKNOWN),
            (Exception("generic"), FailureReasonEnum.UNKNOWN),
            (ValueError("some value error"), FailureReasonEnum.UNKNOWN),
        ],
    )
    def test_maps_correctly(self, db_context: DbContext, exception: BaseException, expected: FailureReasonEnum):
        result = db_context._map_exception_to_failure_reason(type(exception), exception)
        assert result == expected

    def test_validation_error_maps_correctly(self, db_context: DbContext):
        """ValidationError requires special construction so tested separately."""
        from pydantic import BaseModel

        class DummyModel(BaseModel):
            x: int

        try:
            DummyModel(x="not_an_int")  # type: ignore
        except ValidationError as e:
            result = db_context._map_exception_to_failure_reason(type(e), e)
            assert result == FailureReasonEnum.VALIDATION_ERROR

    def test_none_exception_returns_unknown(self, db_context: DbContext):
        result = db_context._map_exception_to_failure_reason(type(None), None)
        assert result == FailureReasonEnum.UNKNOWN

    def test_file_not_found_takes_priority_over_os_error(self, db_context: DbContext):
        """FileNotFoundError is a subclass of OSError — must map to FILE_NOT_FOUND not UNKNOWN."""
        exc = FileNotFoundError("file missing")
        result = db_context._map_exception_to_failure_reason(type(exc), exc)
        assert result == FailureReasonEnum.FILE_NOT_FOUND

    def test_all_enum_values_are_covered(self, db_context: DbContext):
        """Ensures every FailureReasonEnum value except UNKNOWN is reachable via a mapped exception."""
        from pydantic import BaseModel
        from pydantic import ValidationError as PydanticValidationError

        class _Dummy(BaseModel):
            x: int

        validation_exc = None
        try:
            _Dummy(x="not_an_int")  # type: ignore
        except PydanticValidationError as e:
            validation_exc = e

        assert validation_exc is not None, "Failed to construct a ValidationError for testing"

        mapped_results = {
            db_context._map_exception_to_failure_reason(type(exc), exc)
            for exc in [
                FileNotFoundError(),
                DecryptionError(),
                EncryptionError(),
                NetworkError(),
                UploadError(),
                DuplicateTanGError(),
                IncompleteSubmissionError(),
                validation_exc,
            ]
        }
        unmapped = {e for e in FailureReasonEnum if e != FailureReasonEnum.UNKNOWN} - mapped_results
        assert not unmapped, f"These FailureReasonEnum values have no exception mapping: {unmapped}"


class TestDbContextFailureReason:
    def test_file_not_found_maps_correctly(self, ctx, mock_db):
        exc = FileNotFoundError("missing file")
        ctx.__exit__(type(exc), exc, None)
        mock_db.update_submission_state.assert_called_once_with(
            ctx.submission_id,
            SubmissionStateEnum.ERROR,
            data={"error": str(exc)},
            failure_reason=FailureReasonEnum.FILE_NOT_FOUND,
            grzctl_versions=mock.ANY,
        )

    def test_validation_error_maps_correctly(self, ctx, mock_db):
        exc = ValidationError.from_exception_data("test", [])
        ctx.__exit__(type(exc), exc, None)
        mock_db.update_submission_state.assert_called_once_with(
            ctx.submission_id,
            SubmissionStateEnum.ERROR,
            data={"error": str(exc)},
            failure_reason=FailureReasonEnum.VALIDATION_ERROR,
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
