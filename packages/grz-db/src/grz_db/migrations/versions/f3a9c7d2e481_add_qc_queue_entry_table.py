"""add qc_queue table

Revision ID: f3a9c7d2e481
Revises: 9023ab55c239
Create Date: 2026-03-31 12:00:00.000000+00:00

"""

from collections.abc import Sequence

import sqlalchemy as sa
from alembic import op
from sqlmodel.sql.sqltypes import AutoString

revision: str = "f3a9c7d2e481"
down_revision: str | Sequence[str] | None = "9023ab55c239"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    """Upgrade schema."""
    op.create_table(
        "qc_queue",
        sa.Column("submission_id", AutoString(), primary_key=True),
        sa.Column(
            "basic_qc_passed_at",
            sa.DateTime(timezone=True),
            nullable=False,
            server_default=sa.func.now(),
        ),
        sa.ForeignKeyConstraint(["submission_id"], ["submissions.id"], ondelete="CASCADE"),
    )
    op.create_index(
        index_name="ix_qc_queue_basic_qc_passed_at",
        table_name="qc_queue",
        columns=["basic_qc_passed_at"],
    )


def downgrade() -> None:
    """Downgrade schema."""
    raise RuntimeError("Downgrades not supported.")
