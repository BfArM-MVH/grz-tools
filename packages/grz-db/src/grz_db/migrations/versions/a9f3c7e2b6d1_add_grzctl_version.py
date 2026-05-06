"""add grzctl_version to submission_states

Revision ID: a9f3c7e2b6d1
Revises: 834bd50b8734
Create Date: 2026-05-05 00:00:00.000000
"""

from collections.abc import Sequence
from alembic import op
import sqlalchemy as sa

revision: str = "a9f3c7e2b6d1"
down_revision: str | Sequence[str] | None = "834bd50b8734"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    op.add_column(
        "submission_states",
        sa.Column("grzctl_version", sa.String(length=64), nullable=True),
    )
    op.create_index("ix_submission_states_grzctl_version", "submission_states", ["grzctl_version"])


def downgrade() -> None:
    op.drop_index("ix_submission_states_grzctl_version", table_name="submission_states")
    op.drop_column("submission_states", "grzctl_version")
