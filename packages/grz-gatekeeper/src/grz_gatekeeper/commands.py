import logging

import click
import uvicorn
from grz_common.transfer import init_s3_client
from sqlmodel import Session

from . import config
from .config import get_gatekeeper_config
from .dependencies import get_s3_options
from .main import app
from .services import session_cleanup
from .session_store import gatekeeper_db_engine, init_gatekeeper_sessions_db

log = logging.getLogger(__name__)


@click.command()
@click.option("--config-file", help="Path to the configuration file.", required=True)
@click.option("--host", default="127.0.0.1")
@click.option("--port", default=54321)
def run(config_file, host, port):
    """Start GRZ Gatekeeper API server."""
    # Set the global config path before the app starts
    config.CONFIG_FILE_PATH = config_file

    # Initialize session DB table on startup
    init_gatekeeper_sessions_db()

    # Validate config on startup to fail early
    get_gatekeeper_config()
    log.info(f"Starting server with config from: {config.CONFIG_FILE_PATH}")

    uvicorn.run(app, host=host, port=port)


@click.command("clean-stale-sessions")
@click.option("--config-file", help="Path to the configuration file.", required=True)
@click.option("--older-than-days", default=10, show_default=True, help="Clean up sessions older than this many days.")
def clean_stale_sessions(config_file, older_than_days):
    """
    Scan for and remove abandoned/stale upload sessions from grz-gatekeeper database and S3.
    """
    config.CONFIG_FILE_PATH = config_file
    gatekeeper_config = get_gatekeeper_config()
    s3_options = get_s3_options(gatekeeper_config)
    s3_client = init_s3_client(s3_options)

    log.info("Starting stale session cleanup job.")
    with Session(gatekeeper_db_engine) as session:
        session_cleanup(
            session=session, s3_client=s3_client, bucket_name=s3_options.bucket, older_than_days=older_than_days
        )
    log.info("Cleanup job finished.")
