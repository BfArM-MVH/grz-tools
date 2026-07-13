"""add change-request audit columns

Revision ID: c8d4a7f2b1e6
Revises: 66f36abbea34
Create Date: 2026-05-08 00:00:00.000000+00:00

"""

from collections.abc import Sequence

import sqlalchemy as sa
from alembic import op
from sqlmodel.sql.sqltypes import AutoString

# revision identifiers, used by Alembic.
revision: str = "c8d4a7f2b1e6"
down_revision: str | Sequence[str] | None = "66f36abbea34"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    """Upgrade schema."""
    raw_content_type_enum = sa.Enum("PDF", "PNG", name="requestrawcontenttype")
    raw_content_type_enum.create(op.get_bind(), checkfirst=True)

    # All columns are nullable so existing rows (created before these columns existed)
    # remain valid. Required-ness for *new* entries is enforced at the application layer
    # (CLI / Pydantic input validation), not by the schema.
    with op.batch_alter_table("submission_change_requests") as batch_op:
        batch_op.add_column(sa.Column("requester_name", AutoString(), nullable=True))
        batch_op.add_column(sa.Column("requester_email", AutoString(), nullable=True))
        batch_op.add_column(sa.Column("requested_at", sa.Date(), nullable=True))
        batch_op.add_column(sa.Column("request_email_content", sa.Text(), nullable=True))
        batch_op.add_column(sa.Column("request_raw_content", sa.LargeBinary(), nullable=True))
        batch_op.add_column(
            sa.Column(
                "request_raw_content_type",
                sa.Enum("PDF", "PNG", name="requestrawcontenttype", create_type=False),
                nullable=True,
            )
        )
        # The pair-check passes for old rows too (both NULL → equality holds).
        batch_op.create_check_constraint(
            "chk_change_request_raw_content_type_paired",
            "(request_raw_content IS NULL) = (request_raw_content_type IS NULL)",
        )


def downgrade() -> None:
    """Downgrade schema."""
    raise RuntimeError("Downgrades not supported.")
