"""Compression/decompression pipeline stages."""

from __future__ import annotations

import types
import zlib
from typing import Any

from .base import PipelineError, StreamTransformer

# Try to use isal for faster gzip operations (2-4x faster than zlib)
try:
    from isal import isal_zlib

    _ZLIB: types.ModuleType = isal_zlib
    _USING_ISAL: bool = True
except ImportError:
    _ZLIB = zlib
    _USING_ISAL = False


class GzipDecompressor(StreamTransformer):
    """
    Decompresses gzip-compressed data.

    Uses Intel ISA-L (isal) for 2-4x faster decompression if available,
    otherwise falls back to standard zlib.
    """

    # Gzip magic number
    GZIP_MAGIC = b"\x1f\x8b"

    def __init__(self, auto_detect: bool = True, name: str | None = None):
        """
        Initialize the gzip decompressor.

        :param auto_detect: If True, detect gzip format from first bytes.
                           If False, assume all input is gzipped.
        :param name: Stage name for logging
        """
        super().__init__(name or "GzipDecompressor")
        self._auto_detect = auto_detect
        self._is_gzipped: bool | None = None  # None = not yet determined
        self._decompressor: Any = None
        self._header_buffer = b""
        self._bytes_in = 0
        self._bytes_out = 0

    def process(self, data: bytes) -> bytes:
        """
        Process a chunk of potentially compressed data.

        :param data: Input bytes (may be gzip-compressed)
        :returns: Decompressed bytes
        """
        if not data:
            return b""

        self._bytes_in += len(data)

        # Auto-detect gzip on first chunk
        if self._is_gzipped is None and self._auto_detect:
            self._header_buffer += data
            if len(self._header_buffer) < 2:
                # Need more data to detect
                return b""

            if self._header_buffer[:2] == self.GZIP_MAGIC:
                self._is_gzipped = True
                # wbits = 16 + MAX_WBITS for gzip format
                self._decompressor = _ZLIB.decompressobj(16 + zlib.MAX_WBITS)
                backend = "isal" if _USING_ISAL else "zlib"
                self._log.debug(f"Detected gzip stream, using {backend} for decompression")
            else:
                self._is_gzipped = False
                self._log.debug("Stream is not gzip-compressed, passing through")

            # Process buffered data
            data = self._header_buffer
            self._header_buffer = b""
        elif self._is_gzipped is None and not self._auto_detect:
            self._is_gzipped = True
            self._decompressor = _ZLIB.decompressobj(16 + zlib.MAX_WBITS)

        # Pass through if not gzipped
        if not self._is_gzipped:
            self._bytes_out += len(data)
            return data

        # Decompress
        try:
            if self._decompressor is None:
                raise RuntimeError("Decompressor not initialized")
            decompressed = self._decompressor.decompress(data)
            self._bytes_out += len(decompressed)
            return decompressed
        except zlib.error as e:
            raise PipelineError(f"Gzip decompression failed: {e}", self.name, e) from e

    def flush(self) -> bytes:
        """
        Flush any remaining buffered data.

        :returns: Any remaining decompressed data
        """
        # Return any remaining header buffer (non-gzipped case with tiny input)
        if self._header_buffer:
            data = self._header_buffer
            self._header_buffer = b""
            self._bytes_out += len(data)
            return data

        # Flush decompressor
        if self._decompressor is not None:
            try:
                remaining = self._decompressor.flush()
                self._bytes_out += len(remaining)
                return remaining
            except zlib.error as e:
                raise PipelineError(f"Gzip flush failed: {e}", self.name, e) from e

        return b""

    def finalize(self) -> None:
        """Record decompression stats in context."""
        self.context.metadata["gzip_bytes_in"] = self._bytes_in
        self.context.metadata["gzip_bytes_out"] = self._bytes_out
        if self._is_gzipped:
            ratio = self._bytes_out / self._bytes_in if self._bytes_in > 0 else 0
            self.context.metadata["gzip_ratio"] = ratio
            self._log.debug(f"Decompressed {self._bytes_in} -> {self._bytes_out} bytes (ratio: {ratio:.2f})")

    @property
    def is_gzipped(self) -> bool | None:
        """Return whether the stream was detected as gzipped."""
        return self._is_gzipped

    @property
    def bytes_in(self) -> int:
        """Return bytes received."""
        return self._bytes_in

    @property
    def bytes_out(self) -> int:
        """Return bytes produced."""
        return self._bytes_out
