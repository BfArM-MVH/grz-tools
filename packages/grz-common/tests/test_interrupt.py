"""wait_for_files stops the files of a pool one level further on each interrupt."""

import logging
import threading
from concurrent.futures import Future

import pytest
from grz_common.interrupt import TerminateInterrupt, wait_for_files

WAIT = {"wait": True, "cancel_futures": True}


class _Pool:
    """A pool that records its ``shutdown`` calls, and gets a Ctrl-C during the calls numbered in ``interrupted``."""

    def __init__(self, stop: threading.Event, interrupted: tuple[int, ...] = ()):
        self._stop = stop
        self._interrupted = interrupted
        self.calls: list[tuple[dict, bool]] = []

    def shutdown(self, **kwargs):
        # also record whether the running files were told to stop before the call
        self.calls.append((kwargs, self._stop.is_set()))
        if len(self.calls) in self._interrupted:
            raise KeyboardInterrupt


def _raising(exception: BaseException) -> Future:
    """A future whose ``result()`` raises ``exception``, as the wait does that an interrupt reaches."""
    future: Future = Future()
    future.set_exception(exception)
    return future


def test_the_first_ctrl_c_lets_the_running_files_finish(caplog):
    stop = threading.Event()
    pool = _Pool(stop)
    interrupt = KeyboardInterrupt()

    with caplog.at_level(logging.WARNING), pytest.raises(KeyboardInterrupt) as raised:
        wait_for_files(pool, [_raising(interrupt)], stop)

    assert raised.value is interrupt
    assert pool.calls == [(WAIT, False)]
    assert not stop.is_set()
    assert "Ctrl-C again" in caplog.text


def test_a_second_ctrl_c_stops_the_running_files(caplog):
    stop = threading.Event()
    pool = _Pool(stop, interrupted=(1,))
    interrupt = KeyboardInterrupt()

    with caplog.at_level(logging.WARNING), pytest.raises(KeyboardInterrupt) as raised:
        wait_for_files(pool, [_raising(interrupt)], stop)

    assert raised.value is interrupt, "the first interrupt is raised again"
    assert pool.calls == [(WAIT, False), (WAIT, True)]
    assert len(caplog.records) == 2, "each level says what another Ctrl-C does"


def test_a_third_ctrl_c_stops_the_waiting():
    stop = threading.Event()
    pool = _Pool(stop, interrupted=(1, 2))
    interrupt = KeyboardInterrupt()

    with pytest.raises(KeyboardInterrupt) as raised:
        wait_for_files(pool, [_raising(interrupt)], stop)

    assert raised.value is not interrupt, "the third interrupt propagates"
    assert pool.calls == [(WAIT, False), (WAIT, True)]


def test_sigterm_starts_at_the_second_level():
    """A supervisor follows SIGTERM with SIGKILL, so the running files would not get to finish."""
    stop = threading.Event()
    pool = _Pool(stop)
    interrupt = TerminateInterrupt()

    with pytest.raises(TerminateInterrupt) as raised:
        wait_for_files(pool, [_raising(interrupt)], stop)

    assert raised.value is interrupt
    assert pool.calls == [(WAIT, True)]


def test_a_failed_file_lets_the_running_files_finish():
    stop = threading.Event()
    pool = _Pool(stop)
    error = ValueError("the file failed")

    with pytest.raises(ValueError) as raised:
        wait_for_files(pool, [_raising(error)], stop)

    assert raised.value is error
    assert pool.calls == [({"cancel_futures": True}, False)]
