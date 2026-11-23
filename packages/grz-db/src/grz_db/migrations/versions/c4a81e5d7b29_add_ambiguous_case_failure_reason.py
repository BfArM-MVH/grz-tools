"""add ambiguous_case failure reason

Revision ID: c4a81e5d7b29
Revises: fe9309db9fca
Create Date: 2026-11-09 00:00:00.000000+00:00

Adds ``ambiguous_case`` to the PostgreSQL ``failurereasonenum`` type. It is recorded when a
submission's case key matches more than one case, which the resolver refuses to guess.
"""

from collections.abc import Sequence

from alembic import op

revision: str = "c4a81e5d7b29"
down_revision: str | Sequence[str] | None = "f49beda8afac"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    bind = op.get_bind()
    # SQLite stores the enum as plain VARCHAR, so only PostgreSQL needs the type extended.
    # ADD VALUE inside a transaction requires PostgreSQL 12+, and even then the new value
    # cannot be used in the same transaction; this migration does not use it.
    if bind.dialect.name == "postgresql":
        op.execute("ALTER TYPE failurereasonenum ADD VALUE IF NOT EXISTS 'ambiguous_case'")


def downgrade() -> None:
    """Downgrade schema."""
    raise RuntimeError("Downgrades not supported.")
