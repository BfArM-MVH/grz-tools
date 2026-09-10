"""add cases table and link submissions

Revision ID: f8c1a4b7e2d9
Revises: c8d4a7f2b1e6
Create Date: 2026-07-08 00:00:00.000000+00:00

Renames ``submissions.pseudonym`` to ``submissions.local_case_id``, the name of what it holds
(the submitter's ``localCaseId``), introduces the ``cases`` table and links each submission to its
case via ``submissions.case_id``.

A case's authoritative identity is ``psn`` (the RKI pseudonym), unique once assigned.
``submitter_id`` and ``local_case_id`` are resolution keys used to locate a case before a ``psn``
exists, not the authoritative identity themselves. The partial unique index
``ux_cases_submitter_local_case`` keeps the pair unique wherever both halves are present; neither
is required, since a future flow may resolve a case by ``psn`` alone. Existing submissions are
grouped by ``(submitter_id, local_case_id)`` and backfilled into cases. ``submitter_id`` stays on
the submission row as well as on the case.

Some rows are left unlinked: rows missing either half of the key, ``test`` submissions, and rows
whose key is shared by more than one QC-passed ``initial`` submission. A ``test`` submission is
never case-tracked. A shared key identifies no single patient, so no case is created for it and
every submission carrying it keeps ``case_id`` NULL. Those keys are logged.

A shared key leaves NULL meaning not yet resolved, and a later flow can still resolve those
submissions by ``psn``. A ``test`` submission's NULL is permanent instead. Merging distinct
patients into one case has to be unpicked one submission at a time with
``grzctl db case unlink``, and the key that caused the merge is no help in deciding which
submission belongs where.

Also extends ``failurereasonenum`` with ``duplicate_initial``, recorded on a duplicate initial
submission: one being validated for a case whose initial submission already passed basic QC.
"""

import logging
from collections.abc import Sequence

import sqlalchemy as sa
from alembic import op
from sqlmodel.sql.sqltypes import AutoString

# revision identifiers, used by Alembic.
revision: str = "f8c1a4b7e2d9"
down_revision: str | Sequence[str] | None = "c8d4a7f2b1e6"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None

# Mirrors grz_pydantic_models.submission.metadata.LOCAL_CASE_ID_PLACEHOLDERS; hardcoded so the
# migration stays frozen. The archive contains both spellings as redaction placeholders; neither
# identifies a patient, so neither is a case key. The same check runs at resolution time, in
# SubmitterLocalCaseResolver.find_case, via is_redacted_local_case_id.
LOCAL_CASE_ID_NON_KEYS = ["", "REDACTED_LOCAL_CASE_ID"]

# Alembic's own progress logger, so these lines appear in the same stream as "Running upgrade".
logger = logging.getLogger("alembic.runtime.migration")


def _can_key_a_case(rows: sa.FromClause) -> sa.ColumnElement[bool]:
    """Rows carrying both halves of a usable ``(submitter_id, local_case_id)`` key."""
    return sa.and_(
        rows.c.submitter_id.is_not(None),
        rows.c.local_case_id.is_not(None),
        rows.c.local_case_id.not_in(LOCAL_CASE_ID_NON_KEYS),
    )


def _qc_passed_initial(rows: sa.FromClause) -> sa.ColumnElement[bool]:
    """An ``initial`` submission that passed basic QC.

    One patient has exactly one. Earlier uploads that failed basic QC are retries of the same
    patient rather than further patients, which is why this counts QC-passed rows and not every
    ``initial``: counting every one would read a submitter who fixed a failed upload as a
    submitter who reused the key.
    """
    return sa.and_(
        rows.c.submission_type == "initial",
        rows.c.basic_qc_passed.is_(True),
    )


def upgrade() -> None:
    """Upgrade schema."""
    bind = op.get_bind()

    # --- Name the local case ID column for what it holds ---
    # Done first so that the reflection below and the backfill see the new name. The index is
    # dropped and recreated because SQLite cannot rename one.
    op.alter_column("submissions", "pseudonym", new_column_name="local_case_id")
    op.drop_index("ix_submissions_pseudonym", table_name="submissions")
    op.create_index("ix_submissions_local_case_id", "submissions", ["local_case_id"])

    submissions = sa.Table("submissions", sa.MetaData(), autoload_with=bind)

    # Report the keys the backfill below will refuse to group. Logged rather than left for the
    # operator to notice: an upgrade that links only part of the data looks exactly like one that
    # linked all of it.
    untrusted_keys = bind.execute(
        sa.select(
            submissions.c.submitter_id,
            submissions.c.local_case_id,
            sa.func.count().label("n"),
        )
        .where(_can_key_a_case(submissions), _qc_passed_initial(submissions))
        .group_by(submissions.c.submitter_id, submissions.c.local_case_id)
        .having(sa.func.count() > 1)
    ).fetchall()

    if untrusted_keys:
        logger.warning(
            "cases backfill: %d submitter-local case key(s) are shared by more than one QC-passed "
            "'initial' submission, so they identify no single patient. No case is created for "
            "them and every submission carrying one keeps case_id NULL:\n%s",
            len(untrusted_keys),
            "\n".join(f"  ({row[0]}, {row[1]}): {row[2]} QC-passed initial submissions" for row in untrusted_keys),
        )

    # --- Create the cases table ---
    op.create_table(
        "cases",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column("psn", AutoString(), nullable=True),
        sa.Column("submitter_id", AutoString(), nullable=True),
        sa.Column("local_case_id", AutoString(), nullable=True),
    )
    # psn is the authoritative identity: unique when present.
    op.create_index(
        "ux_cases_psn",
        table_name="cases",
        columns=["psn"],
        unique=True,
        postgresql_where=sa.text("psn IS NOT NULL"),
        sqlite_where=sa.text("psn IS NOT NULL"),
    )
    # A submitter's local case ID identifies one patient, so it identifies one case. Unique only
    # where both halves are present: a case resolved by psn alone carries neither, and any number
    # of those may coexist. Without this two processes populating one patient's submissions both
    # find no case and both create one, which leaves the key ambiguous for good.
    op.create_index(
        "ux_cases_submitter_local_case",
        table_name="cases",
        columns=["submitter_id", "local_case_id"],
        unique=True,
        postgresql_where=sa.text("submitter_id IS NOT NULL AND local_case_id IS NOT NULL"),
        sqlite_where=sa.text("submitter_id IS NOT NULL AND local_case_id IS NOT NULL"),
    )

    # --- Link column on submissions ---
    op.add_column("submissions", sa.Column("case_id", sa.Integer(), nullable=True))
    op.create_index("ix_submissions_case_id", "submissions", ["case_id"])
    if bind.dialect.name == "postgresql":
        # SQLite cannot ALTER TABLE ADD CONSTRAINT; the ORM declares the FK regardless.
        op.create_foreign_key(
            "fk_submissions_case_id_cases",
            source_table="submissions",
            referent_table="cases",
            local_cols=["case_id"],
            remote_cols=["id"],
        )

    # --- Backfill: one case per distinct (submitter_id, local_case_id), storing the keys, then link ---
    # Lightweight, standalone Table objects (not reflected) describing just the columns we
    # need. A fresh MetaData per table avoids re-defining the ``submissions`` table reflected above.
    cases = sa.Table(
        "cases",
        sa.MetaData(),
        sa.Column("id", sa.Integer(), primary_key=True),
        sa.Column("submitter_id", AutoString()),
        sa.Column("local_case_id", AutoString()),
    )
    submissions_link = sa.Table(
        "submissions",
        sa.MetaData(),
        sa.Column("case_id", sa.Integer()),
        sa.Column("submitter_id", AutoString()),
        sa.Column("local_case_id", AutoString()),
        # declared as the existing native enum so comparisons render with the right type
        sa.Column(
            "submission_type",
            sa.Enum("initial", "followup", "addition", "correction", "test", name="submissiontype"),
        ),
        sa.Column("basic_qc_passed", sa.Boolean()),
    )

    # The same rule as the report above, asked once per row: how many QC-passed initial
    # submissions share this row's key? More than one and the key denotes more than one patient,
    # so it is not a case key at all. Written as a correlated count rather than an IN list of the
    # keys already fetched, so the statement stays one statement however many keys are bad.
    peers = submissions_link.alias("peers")
    key_denotes_one_patient = (
        sa.select(sa.func.count())
        .select_from(peers)
        .where(
            peers.c.submitter_id == submissions_link.c.submitter_id,
            peers.c.local_case_id == submissions_link.c.local_case_id,
            _qc_passed_initial(peers),
        )
        .correlate(submissions_link)
        .scalar_subquery()
        <= 1
    )

    # Only rows that can form a (submitter_id, local_case_id) key participate in the backfill.
    # Test submissions are never case-tracked; the explicit NULL check keeps rows with an unknown
    # type in the backfill, since SQL's three-valued logic would otherwise exclude them too.
    # A row whose key fails ``key_denotes_one_patient`` is left unlinked whatever its own type or
    # QC state: the key is what is untrustworthy, so the whole group it names stays out, not just
    # the initial submissions that revealed the problem.
    has_keys = sa.and_(
        _can_key_a_case(submissions_link),
        sa.or_(
            submissions_link.c.submission_type.is_(None),
            submissions_link.c.submission_type != "test",
        ),
        key_denotes_one_patient,
    )
    # INSERT INTO cases (submitter_id, local_case_id)
    # SELECT DISTINCT submitter_id, local_case_id FROM submissions WHERE has_keys
    bind.execute(
        sa.insert(cases).from_select(
            ["submitter_id", "local_case_id"],
            sa.select(submissions_link.c.submitter_id, submissions_link.c.local_case_id).where(has_keys).distinct(),
        )
    )
    # Scalar subquery, one per submissions row, that looks up the matching case's id:
    #   SELECT c.id FROM cases c
    #   WHERE c.submitter_id = submissions.submitter_id AND c.local_case_id = submissions.local_case_id
    matching_case_id = (
        sa.select(cases.c.id)
        .where(
            cases.c.submitter_id == submissions_link.c.submitter_id,
            cases.c.local_case_id == submissions_link.c.local_case_id,
        )
        # mark `submissions_link` as coming from the enclosing UPDATE
        .correlate(submissions_link)
        # mark the SELECT as a single-value expression usable in `SET case_id = ...`.
        .scalar_subquery()
    )
    bind.execute(sa.update(submissions_link).values(case_id=matching_case_id).where(has_keys))

    # --- Enforce at most one QC-passed initial submission per case ---
    # Several initial submissions may be linked while pending basic QC (the data alone cannot
    # tell competing uploads apart); whichever passes basic QC first becomes the case's
    # QC-passed initial submission.
    # CREATE UNIQUE INDEX fails if the rows already break the rule, so the backfill above has to
    # have left it intact. Skipping the untrusted keys is what does that: their rows are unlinked,
    # and a unique index does not constrain a NULL case_id.
    op.create_index(
        "ux_submissions_one_initial_per_case",
        table_name="submissions",
        columns=["case_id"],
        unique=True,
        postgresql_where=sa.text("submission_type = 'initial' AND basic_qc_passed IS TRUE"),
        sqlite_where=sa.text("submission_type = 'initial' AND basic_qc_passed IS TRUE"),
    )

    # --- Failure reason for the duplicate initial submission ---
    # SQLite stores the enum as plain VARCHAR, so only PostgreSQL needs the type extended.
    # ADD VALUE inside a transaction requires PostgreSQL 12+, and even then the new value
    # cannot be used in the same transaction; this migration does not use it.
    if bind.dialect.name == "postgresql":
        op.execute("ALTER TYPE failurereasonenum ADD VALUE IF NOT EXISTS 'duplicate_initial'")


def downgrade() -> None:
    """Downgrade schema."""
    raise RuntimeError("Downgrades not supported.")
