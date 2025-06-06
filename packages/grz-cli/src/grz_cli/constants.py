"""Constants for logging configuration, JSON schema validation, and other settings."""

PACKAGE_ROOT = "grz_cli"

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

TQDM_SMOOTHING: float = 0.0001
