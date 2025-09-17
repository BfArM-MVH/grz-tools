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
