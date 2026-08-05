import pytest
from grz_pydantic_models.submission.metadata.versioning import Version


def test_version_comparisons():
    assert Version("1") == Version("1.0")
    assert Version("1.0") == Version("1.0.0.0")
    assert Version("3.2") == Version("3.2.0")

    assert Version("0.1") < Version("1.0.0")

    assert Version("3.2.0") > Version("3.1.99")
    assert Version("3.1.99") < Version("3.2")
    assert Version("3.2") > Version("3.1.99")

    assert Version("1") >= Version("1.0")
    assert Version("2.1.0") >= Version("2")


def test_version_trailing_zeros_do_not_fake_equality():
    """
    A shorter version padded with zeros compares equal, so '<=' must not accept a longer, larger one.
    Comparison used to consume the components, leaving '<=' to compare two empty remainders and call
    them equal, which made every upper bound pass.
    """
    assert Version("1.4.0") > Version("1.3")
    assert not Version("1.4.0") <= Version("1.3")
    assert not Version("1.3.1") <= Version("1.3")
    assert Version("1.3.0") <= Version("1.3")


def test_version_rejects_a_non_numeric_version():
    """Parsing was lazy, so a malformed version only failed later, during an unrelated comparison."""
    with pytest.raises(ValueError, match="Failed to parse"):
        Version("1.2.x")
