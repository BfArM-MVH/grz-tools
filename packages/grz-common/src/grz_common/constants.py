"""Constants for progress bars, JSON schema validation, and other settings."""

TQDM_BAR_FORMAT = "{desc} ▕{bar:50}▏ {n_fmt:>10}/{total_fmt:<10} ({rate_fmt:>12}, ETA: {remaining:>6}) {postfix}"
TQDM_DEFAULTS = {
    "bar_format": TQDM_BAR_FORMAT,
    "unit": "iB",
    "unit_scale": True,
    "miniters": 1,
    "smoothing": 0.00001,
    "colour": "cyan",
}
