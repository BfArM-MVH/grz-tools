import threading
import time

from . import StreamWrapper


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


class MeasuringStream(StreamWrapper):
    """Wraps a stream to measure read latency and throughput."""

    def __init__(self, source, name: str, registry: MetricsRegistry):
        super().__init__(source)
        self.name = name
        self.registry = registry

    def read(self, size: int | None = -1) -> bytes:
        start = time.perf_counter()
        chunk = self.source.read(size)
        duration = time.perf_counter() - start

        self.registry.update(self.name, len(chunk), duration)
        return chunk
