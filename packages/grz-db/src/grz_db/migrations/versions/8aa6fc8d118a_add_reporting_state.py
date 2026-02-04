"""add 'reporting' state

Revision ID: 8aa6fc8d118a
Revises: fb3df229a77b
Create Date: 2026-02-04 11:15:07.740306+00:00

"""

from collections.abc import Sequence

import sqlalchemy as sa
from alembic import op


# revision identifiers, used by Alembic.
revision: str = "8aa6fc8d118a"
down_revision: str | Sequence[str] | None = "fb3df229a77b"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


ENUM_NAME = "submissionstateenum"

OLD_STATES = [
    "UPLOADING",
    "UPLOADED",
    "DOWNLOADING",
    "DOWNLOADED",
    "DECRYPTING",
    "DECRYPTED",
    "VALIDATING",
    "VALIDATED",
    "ENCRYPTING",
    "ENCRYPTED",
    "ARCHIVING",
    "ARCHIVED",
    "REPORTED",
    "QCING",
    "QCED",
    "CLEANING",
    "CLEANED",
    "FINISHED",
    "ERROR",
]

NEW_STATES = [
    "UPLOADING",
    "UPLOADED",
    "DOWNLOADING",
    "DOWNLOADED",
    "DECRYPTING",
    "DECRYPTED",
    "VALIDATING",
    "VALIDATED",
    "ENCRYPTING",
    "ENCRYPTED",
    "ARCHIVING",
    "ARCHIVED",
    "REPORTING",
    "REPORTED",
    "QCING",
    "QCED",
    "CLEANING",
    "CLEANED",
    "FINISHED",
    "ERROR",
]


def upgrade() -> None:
    context = op.get_context()
    dialect_name = context.bind.dialect.name if context.bind else "sqlite"

    # see https://stackoverflow.com/questions/1771543/adding-a-new-value-to-an-existing-enum-type
    # TODO: also check https://bakkenbaeck.com/tech/enums-views-alembic-migrations
    if dialect_name == "postgresql":
        op.execute("COMMIT")
        op.execute(f"ALTER TYPE {ENUM_NAME} ADD VALUE 'REPORTING' BEFORE 'REPORTED'")
    else:
        with op.batch_alter_table("submission_states", schema=None) as batch_op:
            batch_op.alter_column(
                "state",
                type_=sa.Enum(*NEW_STATES, name=ENUM_NAME),
                existing_type=sa.Enum(*OLD_STATES, name=ENUM_NAME),
                existing_nullable=False,
            )


def downgrade() -> None:
    """Downgrade schema."""
    raise RuntimeError("Downgrades not supported.")
