"""
Command modules for the grzctl package.
"""

import click

limit = click.option("--limit", type=click.IntRange(min=0), default=10)
