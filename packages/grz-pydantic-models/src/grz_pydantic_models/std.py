import functools
from warnings import warn


def deprecated(func=None, /, msg=None):
    def decorator(f):
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            if msg is None:
                warn(f"{f.__name__} is deprecated.")
            else:
                warn(msg)
            return f(*args, **kwargs)

        return wrapper

    # stolen from dataclass implementation
    # determine if called as @deprecated or @deprecated()
    if func is None:
        # @deprecated()
        return decorator

    # @deprecated
    return decorator(func)
