import io
import os

import pytest
from grz_common.pipeline.components.perf import MeasuringObserver, MeasuringStream, MetricsRegistry


@pytest.fixture
def registry():
    return MetricsRegistry()


@pytest.fixture
def random_data():
    return os.urandom(1024 * 1024)  # 1MB of random data


def test_measuring_stream_integrity(registry, random_data):
    """
    Verifies that MeasuringStream passes data through exactly 1:1
    and records the read metric.
    """
    source = io.BytesIO(random_data)
    stream = MeasuringStream(source, "test_source", registry)
    output = stream.read()

    assert output == random_data

    stats = registry.metrics.get("test_source")
    assert stats is not None
    assert stats["bytes"] == len(random_data)
    assert stats["time"] > 0


def test_measuring_observer_chaining(registry, random_data):
    """
    Verifies that MeasuringObserver allows data to flow into a sink
    and does not swallow it.
    """
    sink = io.BytesIO()

    observer = MeasuringObserver("test_observer", registry)
    observer.set_sink(sink)

    bytes_written = observer.write(random_data)
    assert bytes_written == len(random_data)
    assert sink.getvalue() == random_data

    stats = registry.metrics.get("test_observer")
    assert stats is not None
    assert stats["bytes"] == len(random_data)


def test_measurement_small_reads(registry):
    """Ensure measuring very small chunks doesn't cause errors."""
    stream = MeasuringStream(io.BytesIO(b"abc"), "small", registry)
    assert stream.read(1) == b"a"
    assert stream.read(1) == b"b"
    assert stream.read(1) == b"c"
    assert stream.read(1) == b""

    assert registry.metrics["small"]["bytes"] == 3
