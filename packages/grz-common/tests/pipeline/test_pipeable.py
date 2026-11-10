"""``>>`` reports the error that stopped the stream, and logs one that closing the stages adds."""

import io
import logging

import pytest
from grz_common.pipeline.components import DevNullSink, Observer, ReadStream, Tee
from grz_common.pipeline.components.perf import StreamMetricsRegistry


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


class _RecordingSink(io.RawIOBase):
    """A sink that records whether the write was completed or aborted."""

    def __init__(self) -> None:
        super().__init__()
        self.events: list[str] = []

    def writable(self) -> bool:
        return True

    def write(self, data) -> int:
        return len(data)

    def close(self) -> None:
        if not self.closed:
            self.events.append("complete")
        super().close()

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        if exc_type is None:
            self.close()
            return
        self.events.append("abort")
        super().close()


def _wrapped_with_measure(sink):
    return StreamMetricsRegistry().measure("upload")(sink)


def _wrapped_with_observer(sink):
    return _NoopObserver() | sink


class _NoopObserver(Observer):
    def observe(self, chunk: bytes) -> None:
        pass


@pytest.mark.parametrize("wrap", [_wrapped_with_measure, _wrapped_with_observer])
def test_a_wrapped_sink_is_aborted_on_error(wrap):
    sink = _RecordingSink()
    with pytest.raises(_StreamError):
        ReadStream(_FailingSource()) >> wrap(sink)
    assert sink.events == ["abort"]


@pytest.mark.parametrize("wrap", [_wrapped_with_measure, _wrapped_with_observer])
def test_a_wrapped_sink_is_completed_on_success(wrap):
    sink = _RecordingSink()
    ReadStream(io.BytesIO(b"data")) >> wrap(sink)
    assert sink.events == ["complete"]
