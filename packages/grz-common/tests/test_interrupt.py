"""wait_for_files stops the files of a pool one level further on each interrupt."""

import logging
import signal
import threading
import time
from concurrent.futures import Future, ThreadPoolExecutor

import pytest
from grz_common.interrupt import TerminateInterrupt, wait_for_files


class _Pool:
    """A pool that records its ``shutdown`` calls."""

    def __init__(self):
        self.calls: list[dict] = []

    def shutdown(self, **kwargs):
        self.calls.append(kwargs)


@pytest.fixture
def waits(mocker):
    """Replace the wait for the running files. ``waits.calls`` records whether they were told to stop before each wait.

    Set ``waits.interrupted`` to the numbers of the waits that get a Ctrl-C.
    """

    class _Waits:
        def __init__(self):
            self.stop = threading.Event()
            self.interrupted: tuple[int, ...] = ()
            self.calls: list[bool] = []

        def __call__(self, futures):
            self.calls.append(self.stop.is_set())
            if len(self.calls) in self.interrupted:
                raise KeyboardInterrupt

    waits = _Waits()
    mocker.patch("grz_common.interrupt.wait", waits)
    return waits


def _raising(exception: BaseException) -> Future:
    """A future whose ``result()`` raises ``exception``, as the wait does that an interrupt reaches."""
    future: Future = Future()
    future.set_exception(exception)
    return future


def test_the_first_ctrl_c_lets_the_running_files_finish(caplog, waits):
    pool = _Pool()
    interrupt = KeyboardInterrupt()

    with caplog.at_level(logging.WARNING), pytest.raises(KeyboardInterrupt) as raised:
        wait_for_files(pool, [_raising(interrupt)], waits.stop)

    assert raised.value is interrupt
    assert pool.calls == [{"wait": False, "cancel_futures": True}], "no queued file starts"
    assert waits.calls == [False]
    assert not waits.stop.is_set()
    assert "Ctrl-C again" in caplog.text


def test_a_second_ctrl_c_stops_the_running_files(caplog, waits):
    waits.interrupted = (1,)
    interrupt = KeyboardInterrupt()

    with caplog.at_level(logging.WARNING), pytest.raises(KeyboardInterrupt) as raised:
        wait_for_files(_Pool(), [_raising(interrupt)], waits.stop)

    assert raised.value is interrupt, "the first interrupt is raised again"
    assert waits.calls == [False, True]
    assert len(caplog.records) == 2, "each level says what another Ctrl-C does"


def test_a_third_ctrl_c_stops_the_waiting(waits):
    waits.interrupted = (1, 2)
    interrupt = KeyboardInterrupt()

    with pytest.raises(KeyboardInterrupt) as raised:
        wait_for_files(_Pool(), [_raising(interrupt)], waits.stop)

    assert raised.value is not interrupt, "the third interrupt propagates"
    assert waits.calls == [False, True]


def test_sigterm_starts_at_the_second_level(waits):
    """A supervisor follows SIGTERM with SIGKILL, so the running files would not get to finish."""
    interrupt = TerminateInterrupt()

    with pytest.raises(TerminateInterrupt) as raised:
        wait_for_files(_Pool(), [_raising(interrupt)], waits.stop)

    assert raised.value is interrupt
    assert waits.calls == [True]


def test_a_failed_file_lets_the_running_files_finish(waits):
    pool = _Pool()
    error = ValueError("the file failed")

    with pytest.raises(ValueError) as raised:
        wait_for_files(pool, [_raising(error)], waits.stop)

    assert raised.value is error
    assert pool.calls == [{"wait": False, "cancel_futures": True}], "no queued file starts"
    assert waits.calls == [False]


def test_a_ctrl_c_after_a_failed_file_stops_the_running_files(waits):
    """The run fails anyway, so the interrupt only cuts the wait short. The error of the file is still raised."""
    waits.interrupted = (1,)
    error = ValueError("the file failed")

    with pytest.raises(ValueError) as raised:
        wait_for_files(_Pool(), [_raising(error)], waits.stop)

    assert raised.value is error
    assert waits.calls == [False, True]


def test_the_second_level_waits_for_a_running_file_to_stop():
    """A real Ctrl-C interrupts the wait itself, which in Python 3.12 must not be ``Thread.join``.

    An interrupted join there marks the running thread as finished, and every later join returns at once.
    """
    stop = threading.Event()
    started = threading.Event()
    cleaned_up = threading.Event()

    def file():
        started.set()
        while not stop.is_set():
            time.sleep(0.01)
        time.sleep(0.2)  # a stopped file aborts its upload and removes its partial file
        cleaned_up.set()

    def ctrl_c_twice():
        started.wait()
        for _ in range(2):
            time.sleep(0.2)
            signal.pthread_kill(threading.main_thread().ident, signal.SIGINT)

    previous_handler = signal.signal(signal.SIGINT, signal.default_int_handler)
    pool = ThreadPoolExecutor(max_workers=1)
    try:
        futures = [pool.submit(file)]
        threading.Thread(target=ctrl_c_twice, daemon=True).start()
        with pytest.raises(KeyboardInterrupt):
            wait_for_files(pool, futures, stop)
        assert cleaned_up.is_set()
    finally:
        stop.set()
        pool.shutdown()
        signal.signal(signal.SIGINT, previous_handler)
