"""``>>`` reports the error that stopped the stream, and logs one that closing the stages adds."""

import io
import logging

import pytest
from grz_common.pipeline.components import DevNullSink, Observer, ReadStream, Tee


class _StreamError(Exception):
    """The error that stops the stream."""


class _Panic(BaseException):
    """Stands in for pyo3's PanicException, which a validator can raise when it is closed."""


class _FailingSource(io.RawIOBase):
    def readable(self) -> bool:
        return True

    def readinto(self, buffer) -> int:
        raise _StreamError("simulated read failure")


class _PanicOnClose(Observer):
    def observe(self, chunk: bytes) -> None:
        pass

    def close(self) -> None:
        if not self.closed:
            try:
                raise _Panic("simulated panic on close")
            finally:
                super().close()


def test_a_close_error_after_a_stream_error_is_logged_not_raised(caplog):
    chain = ReadStream(_FailingSource()) | Tee(_PanicOnClose())

    with caplog.at_level(logging.WARNING), pytest.raises(_StreamError):
        chain >> DevNullSink()

    assert "Closing the pipeline after a failed stream failed as well" in caplog.text
