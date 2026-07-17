"""
Command modules for the grzctl package.
"""

import functools

import click

limit = click.option("--limit", type=click.IntRange(min=0), default=10)


def grzctl_configuration(f):
    """Decorator that reads the unified GrzctlConfig from the click context.

    The root CLI loads the single config file and stores the resulting
    :class:`~grzctl.models.config.GrzctlConfig` in ``ctx.obj["configuration"]``.
    This decorator retrieves it and passes it as the ``configuration`` keyword argument.
    """

    @click.pass_context
    def wrapper(ctx, *args, **kwargs):
        ctx.ensure_object(dict)
        config = ctx.obj.get("configuration")
        if config is None:
            raise click.UsageError("Missing required option '--config'.")
        kwargs["configuration"] = config
        return ctx.invoke(f, *args, **kwargs)

    return functools.update_wrapper(wrapper, f)
