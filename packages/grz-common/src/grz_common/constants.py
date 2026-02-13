"""Constants for progress bars, JSON schema validation, S3 multipart uploads, and other settings."""

import logging

from tqdm.auto import tqdm

TQDM_BAR_FORMAT = "{desc} ▕{bar:50}▏ {n_fmt:>10}/{total_fmt:<10} ({rate_fmt:>12}, ETA: {remaining:>6}) {postfix}"
TQDM_DEFAULTS = {
    "bar_format": TQDM_BAR_FORMAT,
    "unit": "iB",
    "unit_scale": True,
    "miniters": 1,
    "smoothing": 0.00001,
    "colour": "cyan",
    "ascii": "░▒█",
}


class TqdmLoggingHandler(logging.Handler):
    """Logging handler that writes through tqdm to avoid progress bar corruption."""

    def emit(self, record: logging.LogRecord) -> None:
        """Emit a log record using tqdm.write() to avoid progress bar issues."""
        try:
            msg = self.format(record)
            tqdm.write(msg)
        except Exception:
            self.handleError(record)


# Maximum number of parts for multipart upload
MULTIPART_MAX_PARTS = 1000

# Minimum part size
MULTIPART_MIN_PART_SIZE = 5 * 1024 * 1024  # 5 MiB

# Threshold for when to use multipart upload (boto3 default)
MULTIPART_THRESHOLD = 8 * 1024 * 1024  # 8 MiB

# Default chunk/part size for multipart uploads
# Using 256MB as default to reduce number of parts for large files
MULTIPART_DEFAULT_PART_SIZE = 256 * 1024 * 1024  # 256 MiB

# Default chunk size for streaming operations (matches crypt4gh segment size)
STREAMING_CHUNK_SIZE = 64 * 1024  # 64 KiB
