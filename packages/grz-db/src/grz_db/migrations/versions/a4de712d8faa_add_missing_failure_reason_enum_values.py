"""add_missing_failure_reason_enum_values

Revision ID: a4de712d8faa
Revises: 524b2052f143
Create Date: 2026-04-27 12:14:37.441830+00:00

"""

from collections.abc import Sequence

from alembic import op

# revision identifiers, used by Alembic.
revision: str = "a4de712d8faa"
down_revision: str | Sequence[str] | None = "524b2052f143"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    context = op.get_context()
    if context.bind and context.bind.dialect.name == "postgresql":
        for value in ["validation_error", "file_not_found", "encryption_error", "upload_error", "unknown"]:
            op.execute(f"ALTER TYPE failurereasonenum ADD VALUE IF NOT EXISTS '{value}'")


def downgrade() -> None:
    """Downgrade schema."""
    raise RuntimeError("Downgrades not supported.")
