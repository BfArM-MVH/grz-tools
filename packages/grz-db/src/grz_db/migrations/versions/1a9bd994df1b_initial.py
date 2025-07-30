"""initial

Revision ID: 1a9bd994df1b
Revises: 9bf36a91c87e
Create Date: 2025-07-30 09:43:06.671741+00:00

"""

from collections.abc import Sequence

import sqlalchemy as sa
from alembic import op
from sqlmodel.sql.sqltypes import AutoString

# revision identifiers, used by Alembic.
revision: str = "1a9bd994df1b"
down_revision: str | Sequence[str] | None = None
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    """Upgrade schema."""
    op.create_table(
        "submissions",
        sa.Column("id", AutoString(), primary_key=True),
        sa.Column("tan_g", AutoString()),
        sa.Column("pseudonym", AutoString()),
        if_not_exists=True,
    )
    op.create_index(
        index_name="ix_submissions_tan_g", table_name="submissions", columns=["tan_g"], unique=True, if_not_exists=True
    )
    op.create_index(index_name="ix_submissions_id", table_name="submissions", columns=["id"], if_not_exists=True)
    op.create_index(
        index_name="ix_submissions_pseudonym", table_name="submissions", columns=["pseudonym"], if_not_exists=True
    )

    op.create_table(
        "submission_states",
        sa.Column(
            "state",
            sa.Enum(
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
                name="submissionstateenum",
            ),
            nullable=False,
        ),
        sa.Column("data", sa.JSON()),
        sa.Column("timestamp", sa.DateTime(timezone=True), nullable=False),
        sa.Column("id", sa.Integer(), primary_key=True),
        sa.Column("submission_id", AutoString(), sa.ForeignKey("submissions.id"), nullable=False),
        sa.Column("author_name", AutoString(), nullable=False),
        sa.Column("signature", AutoString(), nullable=False),
        if_not_exists=True,
    )
    op.create_index(
        index_name="ix_submission_states_author_name",
        table_name="submission_states",
        columns=["author_name"],
        if_not_exists=True,
    )
    op.create_index(
        index_name="ix_submission_states_submission_id",
        table_name="submission_states",
        columns=["submission_id"],
        if_not_exists=True,
    )
    op.create_index(
        index_name="ix_submission_states_id", table_name="submission_states", columns=["id"], if_not_exists=True
    )

    op.create_table(
        "submission_change_requests",
        sa.Column("change", sa.Enum("MODIFY", "DELETE", "TRANSFER", name="changerequestenum"), nullable=False),
        sa.Column("data", sa.JSON()),
        sa.Column("timestamp", sa.DateTime(timezone=True), nullable=False),
        sa.Column("id", sa.Integer(), primary_key=True),
        sa.Column("submission_id", AutoString(), sa.ForeignKey("submissions.id"), nullable=False),
        sa.Column("author_name", AutoString(), nullable=False),
        sa.Column("signature", AutoString(), nullable=False),
        if_not_exists=True,
    )
    op.create_index(
        index_name="ix_submission_change_requests_author_name",
        table_name="submission_change_requests",
        columns=["author_name"],
        if_not_exists=True,
    )
    op.create_index(
        index_name="ix_submission_change_requests_id",
        table_name="submission_change_requests",
        columns=["id"],
        if_not_exists=True,
    )
    op.create_index(
        index_name="ix_submission_change_requests_submission_id",
        table_name="submission_change_requests",
        columns=["submission_id"],
        if_not_exists=True,
    )


def downgrade() -> None:
    """Downgrade schema."""
    raise RuntimeError("Downgrades not supported.")
