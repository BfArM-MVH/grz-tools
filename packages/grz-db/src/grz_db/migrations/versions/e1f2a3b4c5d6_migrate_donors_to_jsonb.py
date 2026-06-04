"""migrate donors to jsonb

Revision ID: e1f2a3b4c5d6
Revises: d4e5f6a7b8c9
Create Date: 2026-05-30 14:30:00.000000

"""

from collections.abc import Sequence

import sqlalchemy as sa
from alembic import op
import sqlalchemy.dialects.postgresql as sa_psql

# revision identifiers, used by Alembic.
revision: str = "e1f2a3b4c5d6"
down_revision: str | Sequence[str] | None = "d4e5f6a7b8c9"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    # 1. Add new JSONB columns as nullable
    op.add_column("donors", sa.Column("library_types_jsonb", sa.JSON().with_variant(sa_psql.JSONB, "postgresql"), nullable=True))
    op.add_column("donors", sa.Column("sequence_types_jsonb", sa.JSON().with_variant(sa_psql.JSONB, "postgresql"), nullable=True))
    op.add_column("donors", sa.Column("sequence_subtypes_jsonb", sa.JSON().with_variant(sa_psql.JSONB, "postgresql"), nullable=True))
    op.add_column("donors", sa.Column("research_consent_missing_justifications_jsonb", sa.JSON().with_variant(sa_psql.JSONB, "postgresql"), nullable=True))

    # 2. Migrate data
    # SemicolonSeparatedStringSet stores as "val1;val2;val3"
    # We need to convert this to JSON array ["val1", "val2", "val3"]
    # PostgreSQL string_to_array and array_to_json are useful.
    op.execute(
        """
        UPDATE donors
        SET 
            library_types_jsonb = (SELECT json_agg(x) FROM unnest(string_to_array(library_types, ';')) x WHERE x <> ''),
            sequence_types_jsonb = (SELECT json_agg(x) FROM unnest(string_to_array(sequence_types, ';')) x WHERE x <> ''),
            sequence_subtypes_jsonb = (SELECT json_agg(x) FROM unnest(string_to_array(sequence_subtypes, ';')) x WHERE x <> ''),
            research_consent_missing_justifications_jsonb = (SELECT json_agg(x) FROM unnest(string_to_array(research_consent_missing_justifications, ';')) x WHERE x <> '')
        """
    )
    
    # Handle empty sets (null in DB) - ensured by the query above if they are NULL they stay NULL or empty.
    # But wait, if they are NULL, string_to_array(NULL, ';') is NULL. json_agg of NULL is NULL.
    
    # 3. Drop old columns and rename new ones
    op.drop_column("donors", "library_types")
    op.drop_column("donors", "sequence_types")
    op.drop_column("donors", "sequence_subtypes")
    op.drop_column("donors", "research_consent_missing_justifications")
    
    op.alter_column("donors", "library_types_jsonb", new_column_name="library_types", nullable=False)
    op.alter_column("donors", "sequence_types_jsonb", new_column_name="sequence_types", nullable=False)
    op.alter_column("donors", "sequence_subtypes_jsonb", new_column_name="sequence_subtypes", nullable=False)
    op.alter_column("donors", "research_consent_missing_justifications_jsonb", new_column_name="research_consent_missing_justifications", nullable=True)


def downgrade() -> None:
    # 1. Add old columns back as nullable
    op.add_column("donors", sa.Column("library_types_str", sa.String(), nullable=True))
    op.add_column("donors", sa.Column("sequence_types_str", sa.String(), nullable=True))
    op.add_column("donors", sa.Column("sequence_subtypes_str", sa.String(), nullable=True))
    op.add_column("donors", sa.Column("research_consent_missing_justifications_str", sa.String(), nullable=True))

    # 2. Migrate data back
    # Convert JSON array ["val1", "val2"] back to "val1;val2"
    op.execute(
        """
        UPDATE donors
        SET 
            library_types_str = (SELECT array_to_string(array(SELECT json_array_elements_text(library_types)), ';')),
            sequence_types_str = (SELECT array_to_string(array(SELECT json_array_elements_text(sequence_types)), ';')),
            sequence_subtypes_str = (SELECT array_to_string(array(SELECT json_array_elements_text(sequence_subtypes)), ';')),
            research_consent_missing_justifications_str = (SELECT array_to_string(array(SELECT json_array_elements_text(research_consent_missing_justifications)), ';'))
        """
    )

    # 3. Drop JSONB columns and rename new ones
    op.drop_column("donors", "library_types")
    op.drop_column("donors", "sequence_types")
    op.drop_column("donors", "sequence_subtypes")
    op.drop_column("donors", "research_consent_missing_justifications")

    op.alter_column("donors", "library_types_str", new_column_name="library_types", nullable=False)
    op.alter_column("donors", "sequence_types_str", new_column_name="sequence_types", nullable=False)
    op.alter_column("donors", "sequence_subtypes_str", new_column_name="sequence_subtypes", nullable=False)
    op.alter_column("donors", "research_consent_missing_justifications_str", new_column_name="research_consent_missing_justifications", nullable=True)
