import threading
import time
from collections.abc import Buffer

from . import ReadStream, WriteStream


class MetricsRegistry:
    """Thread-safe registry to aggregate metrics."""

    def __init__(self):
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


class MeasuringReadStream(ReadStream):
    """Wraps a ReadStream to measure read latency and throughput (Pull)."""

    def __init__(self, source, name: str, registry: MetricsRegistry):
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

    def __init__(self, name: str, registry: MetricsRegistry):
        super().__init__()
        self.name = name
        self.registry = registry

    def write(self, data: Buffer) -> int:
        start = time.perf_counter()
        # Measure total time (own processing + downstream push)
        bytes_written = super().write(data)
        duration = time.perf_counter() - start

        self.registry.update(self.name, bytes_written, duration)
        return bytes_written
