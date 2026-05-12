"""add failure_reason to submission_states

Revision ID: 09602efd9105
Revises: 834bd50b8734
Create Date: 2026-04-29 09:35:54.981907+00:00
"""

from collections.abc import Sequence

import sqlalchemy as sa
from alembic import op

revision: str = "09602efd9105"
down_revision: str | Sequence[str] | None = "c5f6a8b9d2e1"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None

ENUM_NAME = "failurereasonenum"
FAILURE_REASONS = [
    "duplicate_tang",
    "incomplete_submission",
    "decryption_error",
    "network_error",
    "validation_error",
    "file_not_found",
    "encryption_error",
    "upload_error",
    "unknown",
]


def upgrade() -> None:
    """Upgrade schema."""
    context = op.get_context()
    dialect_name = context.bind.dialect.name if context.bind else "sqlite"

    if dialect_name == "postgresql":
        failure_reason_enum = sa.Enum(*FAILURE_REASONS, name=ENUM_NAME)
        failure_reason_enum.create(op.get_bind(), checkfirst=True)
        op.add_column(
            "submission_states",
            sa.Column("failure_reason", sa.Enum(*FAILURE_REASONS, name=ENUM_NAME), nullable=True),
        )
    else:
        with op.batch_alter_table("submission_states", schema=None) as batch_op:
            batch_op.add_column(sa.Column("failure_reason", sa.Enum(*FAILURE_REASONS, name=ENUM_NAME), nullable=True))


def downgrade() -> None:
    """Downgrade schema."""
    raise RuntimeError("Downgrades not supported.")
