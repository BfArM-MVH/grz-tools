"""Constants for logging configuration, JSON schema validation, and other settings."""

# This will be overridden by the package that imports it
PACKAGE_ROOT = "grz_common"

LOGGING_FORMAT = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
LOGGING_DATEFMT = "%Y-%m-%d %I:%M %p"
LOGGING_CONFIG = {
    "version": 1,
    "disable_existing_loggers": True,
    "formatters": {
        "standard": {
            "format": LOGGING_FORMAT,
            "datefmt": LOGGING_DATEFMT,
        },
    },
    "handlers": {
        "default": {
            "level": "INFO",
            "formatter": "standard",
            "class": "logging.StreamHandler",
        },
    },
    "loggers": {
        "": {  # root logger
            "handlers": ["default"],
            "level": "WARNING",
            "propagate": False,
        },
        PACKAGE_ROOT: {"handlers": ["default"], "level": "INFO", "propagate": False},
        "__main__": {  # if __name__ == '__main__'
            "handlers": ["default"],
            "level": "DEBUG",
            "propagate": False,
        },
    },
}

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
