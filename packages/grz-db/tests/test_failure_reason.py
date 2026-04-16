"""Tests for FailureReasonEnum."""

from grz_db.models.submission import FailureReasonEnum


def test_failure_reason_enum_values():
    """Test that all enum values are correctly defined."""
    assert FailureReasonEnum.DUPLICATE_TANG == "Duplicate TanG"
    assert FailureReasonEnum.MISSING_DATA == "Missing Data"
    assert FailureReasonEnum.DECRYPTION_ERROR == "Decryption Error"
    assert FailureReasonEnum.NETWORK_ERROR == "Network Error"


def test_failure_reason_enum_case_insensitive():
    """Test case insensitive enum behavior."""
    assert FailureReasonEnum("duplicate tang") == FailureReasonEnum.DUPLICATE_TANG
    assert FailureReasonEnum("MISSING DATA") == FailureReasonEnum.MISSING_DATA
    assert FailureReasonEnum("Decryption Error") == FailureReasonEnum.DECRYPTION_ERROR


def test_failure_reason_enum_listable():
    """Test that the enum can be listed for CLI choices."""
    reasons = FailureReasonEnum.list()
    assert len(reasons) == 4
    assert "Duplicate TanG" in reasons
    assert "Missing Data" in reasons
    assert "Decryption Error" in reasons
    assert "Network Error" in reasons
