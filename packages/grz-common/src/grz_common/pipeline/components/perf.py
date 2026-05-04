import threading
import time
from collections.abc import Buffer
from typing import Any

from . import Readable, ReadStream, Writable, WriteStream


class StreamMetricsRegistry:
    """Thread-safe registry to aggregate metrics."""

    def __init__(self, enabled: bool = True):
        self.enabled = enabled
        self._lock = threading.Lock()
        self.metrics: dict[str, dict[str, float]] = {}

    def update(self, name: str, size: int, duration: float):
        with self._lock:
            if name not in self.metrics:
                self.metrics[name] = {"bytes": 0, "time": 0.0}
            self.metrics[name]["bytes"] += size
            self.metrics[name]["time"] += duration

    def report(self):
        """Log the throughput per stage."""
        with self._lock:
            stats = []
            for name, data in self.metrics.items():
                mb = data["bytes"] / (1024 * 1024)
                if data["time"] > 0:
                    mb_s = mb / data["time"]
                    stats.append(f"{name}: {mb:.2f}MB in {data['time']:.2f}s ({mb_s:.2f} MB/s)")
            return " | ".join(stats)

    def measure(self, name: str, stream: Any = None) -> type | Any | None:
        """Return a measuring wrapper for pipeline integration.

        Without a stream, returns a class for use with ``|``::

            pipeline |= metrics.measure("1_Source")
            pipeline = pipeline | Encrypt(...) | metrics.measure("4_Encrypt")

        With a stream, wraps it directly (for non-Pipeable objects like file handles)::

            writer = metrics.measure("2b_Write", writer)

        Returns None (or the original stream) when metrics are disabled,
        which ``Pipeable.__or__`` treats as a no-op.
        """
        if not self.enabled:
            return stream  # None for pipe usage, original stream for direct wrap

        registry = self

        class _MeasuringStage:
            """Factory dispatching to MeasuringReadStream or MeasuringWriteStream."""

            def __new__(cls, s: Any) -> Any:
                if isinstance(s, Readable) and s.readable():
                    return MeasuringReadStream(s, name, registry)
                if isinstance(s, Writable) and s.writable():
                    return MeasuringWriteStream(s, name, registry)
                raise TypeError(f"Cannot measure stream of type {type(s).__name__}")

        if stream is not None:
            return _MeasuringStage(stream)
        return _MeasuringStage


class MeasuringReadStream(ReadStream):
    """Wraps a ReadStream to measure read latency and throughput (Pull)."""

    def __init__(self, source: Readable | None, name: str, registry: StreamMetricsRegistry):
        super().__init__(source)
        self.name = name
        self.registry = registry

    def read(self, size: int | None = -1) -> bytes:
        start = time.perf_counter()
        if self.source is None:
            raise RuntimeError("Stream source not set")
        chunk = self.source.read(size)
        duration = time.perf_counter() - start

        self.registry.update(self.name, len(chunk), duration)
        return chunk


class MeasuringWriteStream(WriteStream):
    """Wraps a WriteStream to measure processing time (Push)."""

    def __init__(self, sink: Writable | None, name: str, registry: StreamMetricsRegistry):
        super().__init__(sink)
        self.name = name
        self.registry = registry

    def write(self, data: Buffer) -> int:
        start = time.perf_counter()
        # Measure total time (own processing + downstream push)
        bytes_written = super().write(data)
        duration = time.perf_counter() - start

        self.registry.update(self.name, bytes_written, duration)
        return bytes_written
