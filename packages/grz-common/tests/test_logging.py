import logging
from unittest.mock import patch

from grz_common import logging as grz_logging


def _format_test_record(logging_format: str) -> tuple[str, str, str]:
    record_name = "test.logger"
    record_level = logging.INFO
    record_message = "hello"
    record = logging.LogRecord(record_name, record_level, __file__, 123, record_message, (), None)
    formatted = logging.Formatter(fmt=logging_format, datefmt=grz_logging.LOGGING_DATEFMT).format(record)
    return record_name, record_message, formatted


def test_build_logging_format_includes_hostname_when_available():
    logging_format = grz_logging.build_logging_format("mock-host")

    record_name, record_message, formatted = _format_test_record(logging_format)

    assert "mock-host" in logging_format
    assert f"mock-host {record_name}: {record_message}" in formatted


def test_build_logging_format_omits_hostname_when_none():
    logging_format = grz_logging.build_logging_format(None)

    record_name, record_message, formatted = _format_test_record(logging_format)

    assert f"{record_name}: {record_message}" in formatted
    assert "None" not in logging_format


def test_get_hostname_returns_name_when_available():
    with patch("socket.gethostname", return_value="mock-host"):
        assert grz_logging.get_hostname() == "mock-host"


def test_get_hostname_returns_none_when_unavailable():
    with patch("socket.gethostname", side_effect=OSError):
        assert grz_logging.get_hostname() is None
