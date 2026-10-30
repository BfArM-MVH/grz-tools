"""add grz_cli_version to submissions

NULL = never recorded (submissions run before this change or with --no-update-db)
VERSION FILE NOT FOUND = Signifies that the grz-cli version object was not found. Constant defined in transfer.py
VERSION READ FAILED = Signifies an error occurred during the grz-cli version check. Constant defined in transfer.py

Revision ID: f49beda8afac
Revises: fe9309db9fca
Create Date: 2026-10-23 13:38:45.251444+00:00

"""

from collections.abc import Sequence

import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision: str = "f49beda8afac"
down_revision: str | Sequence[str] | None = "fe9309db9fca"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    """Upgrade schema."""
    op.add_column("submissions", sa.Column("grz_cli_version", sa.String(), nullable=True))


def downgrade() -> None:
    """Downgrade schema."""
    raise RuntimeError("Downgrades not supported.")
