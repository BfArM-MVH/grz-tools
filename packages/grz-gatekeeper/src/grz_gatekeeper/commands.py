import logging

import click
import uvicorn

from . import config
from .config import get_gatekeeper_config
from .main import app
from .session_store import init_gatekeeper_sessions_db

log = logging.getLogger(__name__)


@click.command()
@click.option("--config-file", help="Path to the configuration file.", required=True)
@click.option("--host", default="127.0.0.1")
@click.option("--port", default=54321)
def run(config_file, host, port):
    """Starts the GRZ Gatekeeper API server."""
    # Set the global config path before the app starts
    config.CONFIG_FILE_PATH = config_file

    # Initialize session DB table on startup
    init_gatekeeper_sessions_db()

    # Validate config on startup to fail early
    get_gatekeeper_config()
    log.info(f"Starting server with config from: {config.CONFIG_FILE_PATH}")

    uvicorn.run(app, host=host, port=port)
