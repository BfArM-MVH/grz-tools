import logging
from importlib import reload
from unittest.mock import patch

import pytest
from grz_common import logging as grz_logging


@pytest.fixture(autouse=True)
def _restore_logging_module():
    yield
    reload(grz_logging)


def _create_and_format_test_record() -> tuple[str, str, str]:
    record_name = "test.logger"
    record_level = logging.INFO
    record_message = "hello"
    record = logging.LogRecord(record_name, record_level, __file__, 123, record_message, (), None)
    formatted = logging.Formatter(fmt=grz_logging.LOGGING_FORMAT, datefmt=grz_logging.LOGGING_DATEFMT).format(record)
    return record_name, record_message, formatted


def test_setup_cli_logging_adds_hostname_to_log_records_when_hostname_is_available():
    with patch("socket.gethostname", return_value="mock-host"):
        reload(grz_logging)
        grz_logging.setup_cli_logging(log_file=None, log_level="INFO")

        record_name, record_message, formatted = _create_and_format_test_record()

        assert grz_logging.HOSTNAME == "mock-host"
        assert grz_logging.HOSTNAME in grz_logging.LOGGING_FORMAT
        assert f"{grz_logging.HOSTNAME} {record_name}: {record_message}" in formatted


def test_setup_cli_logging_omits_hostname_when_unavailable():
    with patch("socket.gethostname", side_effect=OSError):
        reload(grz_logging)
        grz_logging.setup_cli_logging(log_file=None, log_level="INFO")

        record_name, record_message, formatted = _create_and_format_test_record()

        assert grz_logging.HOSTNAME is None
        assert f"{record_name}: {record_message}" in formatted
