import itertools
import logging
from operator import itemgetter

import rich.pretty
import rich.table
import rich.text
import sqlalchemy.orm
import textual
from cryptography.hazmat.primitives.serialization import SSHPublicKeyTypes
from grz_db.models.submission import Submission, SubmissionDb, SubmissionStateLog
from grz_pydantic_models.submission.metadata import SubmissionType
from sqlalchemy import func as sqlfn
from sqlmodel import select
from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import HorizontalScroll
from textual.validation import Number, Regex
from textual.widgets import DataTable, Footer, Input, Static, TabbedContent, TabPane

from . import _verify_signature

logger = logging.getLogger(__name__)
_DEFAULT_SEARCH_LIMIT = 20


class SubmissionCountByStateTable(Static):
    def on_mount(self) -> None:
        self.loading = True
        self.border_title = "Submission count by state"
        self.styles.border = ("round", self.app.theme_variables["foreground"])
        self.styles.min_width = 31

    @textual.work
    async def load(self, database: SubmissionDb) -> None:
        self.loading = True
        with database._get_session() as session:
            # credit: https://stackoverflow.com/a/28090544
            log1 = sqlalchemy.orm.aliased(SubmissionStateLog, name="log1")
            log2 = sqlalchemy.orm.aliased(SubmissionStateLog, name="log2")

            statement = (
                select(log1.state, sqlfn.count(log1.state))
                .join(
                    log2, (log1.submission_id == log2.submission_id) & (log1.timestamp < log2.timestamp), isouter=True
                )
                .where(log2.timestamp.is_(None))
                .group_by(log1.state)
            )

            counts_by_state = session.exec(statement).all()

        table = rich.table.Table(
            rich.table.Column(header="State"),
            rich.table.Column(header="Submission Count", justify="center"),
            show_header=False,
        )
        total = 0
        for state, count in counts_by_state:
            table.add_row(str(state), rich.pretty.Pretty(count))
            total += count
        table.add_section()
        table.add_row("[bold]Total[/bold]", rich.pretty.Pretty(total))
        self.update(table)
        self.loading = False


class SubmissionCountByConsentTable(Static):
    def on_mount(self) -> None:
        self.loading = True
        self.border_title = "Submission count by consent"
        self.styles.border = ("round", self.app.theme_variables["foreground"])
        self.styles.min_width = 33

    @textual.work
    async def load(self, database: SubmissionDb) -> None:
        self.loading = True
        with database._get_session() as session:
            statement = select(Submission.consented, sqlfn.count(Submission.consented)).group_by(Submission.consented)
            counts_by_consent = session.exec(statement).all()

        table = rich.table.Table(
            rich.table.Column(header="Consented"),
            rich.table.Column(header="Submission Count", justify="center"),
            show_header=False,
        )
        for consented, count in counts_by_consent:
            table.add_row(rich.pretty.Pretty(consented), rich.pretty.Pretty(count))
        self.update(table)
        self.loading = False


class SubmissionCountByDetailedQCByLETable(Static):
    def on_mount(self) -> None:
        self.loading = True
        self.border_title = "Detailed QC by LE (Non-Test)"
        self.styles.border = ("round", self.app.theme_variables["foreground"])

    @textual.work
    async def load(self, database: SubmissionDb) -> None:
        self.loading = True
        with database._get_session() as session:
            statement = (
                select(
                    Submission.submitter_id,
                    Submission.detailed_qc_passed,
                    sqlfn.count(Submission.submitter_id),  # can't count detailed_qc_pass because NULLs are not counted
                )
                .where(Submission.submission_type != SubmissionType.test)
                .group_by(Submission.submitter_id, Submission.detailed_qc_passed)
            )
            counts_by_pass_by_le = session.exec(statement).all()

        rows = []
        for submitter_id, group in itertools.groupby(
            sorted(counts_by_pass_by_le, key=itemgetter(0)), key=itemgetter(0)
        ):
            qced = 0
            total = 0
            for _, detailed_qc_passed, count in group:
                if detailed_qc_passed is not None:
                    qced += count
                total += count
            rows.append([submitter_id, qced, total])

        table = rich.table.Table(
            rich.table.Column(header="Submitter ID", justify="center"),
            rich.table.Column(header="QCed", justify="center"),
        )
        for submitter_id, qced, total in rows:
            qced_prop = qced / total
            qced_prop_text = (
                rich.text.Text(f"{qced}/{total} (")
                + rich.text.Text(f"{qced_prop:.1%}", style="green" if qced_prop >= 0.02 else "red")
                + rich.text.Text(")")
            )
            table.add_row(str(submitter_id), qced_prop_text)
        self.update(table)
        self.loading = False


class SearchResultsDataTable(DataTable):
    def on_mount(self) -> None:
        self.loading = True
        self.add_columns("Submission ID", "Pseudonym")

    @textual.work
    async def search(
        self,
        database: SubmissionDb,
        submission_id: str | None = None,
        pseudonym: str | None = None,
        submitter_id: str | None = None,
        limit: int | None = None,
    ) -> None:
        self.loading = True
        self.clear()
        with database._get_session() as session:
            statement = select(Submission.id, Submission.pseudonym)
            if submission_id:
                statement = statement.where(Submission.id == submission_id)
            if pseudonym:
                statement = statement.where(Submission.pseudonym == pseudonym)
            if submitter_id:
                statement = statement.where(Submission.submitter_id == submitter_id)
            if isinstance(limit, int) and limit >= 0:
                statement = statement.limit(limit)
            else:
                raise ValueError("limit must be >=0")
            results = session.exec(statement).all()
        logger.debug("Populating search table from the following statement: '%s'", str(statement))
        self.add_rows(results)

        self.loading = False


class DatabaseBrowser(App):
    """A Textual app to browse a GRZ submission database."""

    BINDINGS = [
        Binding(key="r", action="refresh", description="Refresh data"),
        Binding(key="q", action="quit", description="Quit the app"),
    ]
    CSS_PATH = "tui.css"

    def __init__(self, database: SubmissionDb, public_keys: dict[str, SSHPublicKeyTypes], **kwargs) -> None:
        super().__init__(**kwargs)
        self._database = database
        self._public_keys = public_keys

    def compose(self) -> ComposeResult:
        with TabbedContent():
            with TabPane("Overview", id="pane-overview"):
                with HorizontalScroll(id="container-overview-charts"):
                    yield SubmissionCountByStateTable()
                    yield SubmissionCountByConsentTable()
                    yield SubmissionCountByDetailedQCByLETable()
                yield DataTable(id="table-states-latest")
            with TabPane("Search", id="pane-search"):
                with HorizontalScroll(id="container-search-inputs"):
                    input_id = Input(
                        id="search-input-id",
                        placeholder="123456789_2025-07-01_a2b2c3d4",
                        validators=[Regex(r"^[0-9]{9}_\d{4}-\d{2}-\d{2}_[a-f0-9]{8}$")],
                    )
                    input_id.border_title = "Submission ID"
                    # explicitly validate here because is_valid initializes to True and validation only automatically runs on change events
                    input_id.validate(input_id.value)
                    yield input_id
                    input_pseudonym = Input(id="search-input-pseudonym", placeholder="CASE12345")
                    input_pseudonym.border_title = "Pseudonym"
                    input_pseudonym.validate(input_pseudonym)
                    yield input_pseudonym
                    input_submitter = Input(
                        id="search-input-submitter", placeholder="123456789", validators=[Regex(r"^[0-9]{9}$")]
                    )
                    input_submitter.border_title = "Submitter"
                    input_submitter.validate(input_submitter.value)
                    yield input_submitter
                    input_limit = Input(
                        id="search-input-limit", placeholder=str(_DEFAULT_SEARCH_LIMIT), validators=[Number(minimum=0)]
                    )
                    input_limit.border_title = "Limit"
                    input_limit.validate(input_limit.value)
                    yield input_limit
                yield SearchResultsDataTable(id="table-search-results")
        yield Footer()

    def on_mount(self) -> None:
        for table in self.query(DataTable):
            table.loading = True
        table_states_latest = self.query_exactly_one("#table-states-latest")
        table_states_latest.add_columns("Timestamp", "State", "Submission ID", "Steward", "Signature")
        table_states_latest.cursor_type = "row"
        self.action_refresh()

    @textual.work
    async def _refresh_overview(self) -> None:
        self.query_exactly_one(SubmissionCountByStateTable).load(self._database)
        self.query_exactly_one(SubmissionCountByConsentTable).load(self._database)
        self.query_exactly_one(SubmissionCountByDetailedQCByLETable).load(self._database)
        table_states_latest = self.query_exactly_one("#table-states-latest")
        table_states_latest.loading = True
        with self._database._get_session() as session:
            statement = (
                select(SubmissionStateLog).order_by(SubmissionStateLog.timestamp.desc()).limit(_DEFAULT_SEARCH_LIMIT)
            )
            latest_states = session.exec(statement).all()
        table_states_latest.clear()
        for state in latest_states:
            signature_status, verifying_key_comment = _verify_signature(self._public_keys, state.author_name, state)
            table_states_latest.add_row(
                state.timestamp,
                state.state,
                state.submission_id,
                state.author_name,
                signature_status.rich_display(verifying_key_comment),
            )
        table_states_latest.loading = False

    @textual.on(Input.Submitted)
    def _on_input_submitted(self, event: Input.Submitted) -> None:
        self._refresh_search()

    @textual.work
    async def _refresh_search(self) -> None:
        input_id = self.query_exactly_one("#search-input-id")
        submission_id = input_id.value if input_id.is_valid else None

        input_pseudonym = self.query_exactly_one("#search-input-pseudonym")
        pseudonym = input_pseudonym.value if input_pseudonym.is_valid else None

        input_submitter = self.query_exactly_one("#search-input-submitter")
        submitter = input_submitter.value if input_submitter.is_valid else None

        input_limit = self.query_exactly_one("#search-input-limit")
        # default some limit in case of invalid input instead of listing whole database
        limit = int(input_limit.value) if input_limit.is_valid else _DEFAULT_SEARCH_LIMIT

        logger.debug("search(submission_id = '%s')", submission_id)
        self.query_exactly_one(SearchResultsDataTable).search(
            database=self._database,
            submission_id=submission_id,
            pseudonym=pseudonym,
            submitter_id=submitter,
            limit=limit,
        )

    def action_refresh(self) -> None:
        self._refresh_overview()
        self._refresh_search()
