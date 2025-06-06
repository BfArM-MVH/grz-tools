import functools
from warnings import warn


def deprecated(msg=None):
    def decorator(f):
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            if msg is None:
                warn(f"{f.__name__} is deprecated.")
            else:
                warn(msg)
            return f(*args, **kwargs)

        return wrapper

    return decorator
