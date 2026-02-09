"""Checksum computation pipeline stages."""

from __future__ import annotations

import hashlib

from .base import StreamObserver


class Sha256Checksummer(StreamObserver):
    """
    Computes SHA256 checksum of data passing through.

    Usage:
        checksummer = Sha256Checksummer()
        checksummer.initialize(context)
        checksummer.observe(data1)
        checksummer.observe(data2)
        checksummer.finalize()
        checksum = checksummer.hexdigest()
    """

    def __init__(self, name: str | None = None):
        super().__init__(name or "Sha256Checksummer")
        self._hasher = hashlib.sha256()
        self._bytes_processed = 0

    def observe(self, data: bytes) -> None:
        """Update the checksum with observed data."""
        self._hasher.update(data)
        self._bytes_processed += len(data)

    def finalize(self) -> None:
        """Store the final checksum in the context."""
        self.context.checksums["sha256"] = self._hasher.hexdigest()

    def hexdigest(self) -> str:
        """Return the current checksum as a hex string."""
        return self._hasher.hexdigest()

    def digest(self) -> bytes:
        """Return the current checksum as bytes."""
        return self._hasher.digest()

    @property
    def bytes_processed(self) -> int:
        """Return the number of bytes processed."""
        return self._bytes_processed

    def get_result(self) -> str:
        """Return the checksum hex digest."""
        return self.hexdigest()
