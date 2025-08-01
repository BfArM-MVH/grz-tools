"""Constants for logging configuration, JSON schema validation, and other settings."""

LOGGING_FORMAT = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
LOGGING_DATEFMT = "%Y-%m-%d %I:%M %p"

TQDM_SMOOTHING: float = 0.00001
TQDM_BAR_FORMAT = "{desc} ▕{bar:50}▏ {n_fmt:>10}/{total_fmt:<10} ({rate_fmt:>12}, ETA: {remaining:>6}) {postfix}"
TQDM_DEFAULTS = {
    "bar_format": TQDM_BAR_FORMAT,
    "unit": "iB",
    "unit_scale": True,
    "miniters": 1,
    "smoothing": TQDM_SMOOTHING,
    "colour": "cyan",
}

REDACTED_TAN = "0" * 64
