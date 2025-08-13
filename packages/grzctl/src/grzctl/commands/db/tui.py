import json
import logging
from datetime import datetime

import rich.panel
import rich.table
import rich.text
from grz_db.models.submission import Submission, SubmissionStateEnum
from textual import events
from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal, Vertical
from textual.widgets import DataTable, Footer, Header, Input, Static
from textual.widgets._data_table import ColumnKey

from .common import _verify_signature

log = logging.getLogger(__name__)


class SubmissionBrowser(App):
    CSS_PATH = "tui.css"
    BINDINGS = [
        Binding("q", "quit", "Quit"),
        Binding("f", "focus_filter", "Filter"),
        Binding("/", "focus_filter", "Filter", show=False),
        Binding("s", "toggle_sort_column", "Sort"),
    ]

    ID_COLUMNS = {
        "id": "ID",
    }
    STATE_COLUMNS = {
        "latest_state": "Latest State",
        "timestamp": "Timestamp",
        "data_steward": "Data Steward",
        "signature": "Signature",
    }
    DETAIL_COLUMNS = {
        "tan_g": "tanG",
        "pseudonym": "Pseudonym",
        "submission_date": "Subm. Date",
        "submission_type": "Subm. Type",
        "submitter_id": "Submitter ID",
        "data_node_id": "Data Node",
        "disease_type": "Disease",
        "library_type": "Library",
        "consented": "Consented",
        "basic_qc_passed": "Basic QC",
        "detailed_qc_passed": "Detailed QC",
    }

    def __init__(self, submissions: list[Submission], public_keys: dict, **kwargs) -> None:
        super().__init__(**kwargs)
        self.submission_map = {s.id: s for s in submissions}
        self.public_keys = public_keys
        self.title = "GRZ DB Submission Browser"
        self._sort_by: ColumnKey | None = None
        self._sort_reverse = False
        self._labels: dict[str, str] = {}
        self._table: DataTable | None = None
        self._current_submission_id: str | None = None

    def compose(self) -> ComposeResult:
        yield Header()
        yield Input(placeholder="Filter submissions...", id="filter-input")
        with Horizontal():
            table = DataTable(id="list-table", cursor_type="cell")
            # start the app with focus on the table
            table.focus()
            yield table
            with Vertical(id="details-view"):
                yield Static(id="details-panel", expand=True)
                yield DataTable(id="history-table", cursor_type="row")
        yield Footer()

    def on_mount(self) -> None:
        list_table = self.query_one("#list-table", DataTable)
        self._table = list_table
        # submission_id is always first!
        self._labels = self.ID_COLUMNS | self.STATE_COLUMNS | self.DETAIL_COLUMNS
        for key, label in self._labels.items():
            list_table.add_column(label, key=key)

        history_table = self.query_one("#history-table", DataTable)
        history_table.add_columns("Log ID", "Timestamp (UTC)", "State", "Data", "Data Steward", "Signature")

        self.populate_table()

        if self._table.row_count > 0:
            first_row_key = self._table.get_row_at(0)[0]
            if first_row_key is not None and (submission := self.submission_map.get(str(first_row_key))):
                self._current_submission_id = str(first_row_key)
                self._update_details_view(submission)
        else:
            self.query_one("#details-panel", Static).update(
                rich.panel.Panel("No submissions to display.", title="Details", border_style="dim")
            )

        self.query_one("#filter-input").display = False

    def on_data_table_cell_highlighted(self, event: DataTable.CellHighlighted) -> None:
        if event.control.id != "list-table" or event.cell_key is None:
            return

        try:
            submission_id = self._table.get_row(event.cell_key.row_key)[0]
            if submission_id == self._current_submission_id:
                return
            self._current_submission_id = submission_id
            if submission := self.submission_map.get(submission_id):
                self._update_details_view(submission)
        except KeyError:
            log.debug("Could not find row key %s during highlight.", event.cell_key.row_key)

    def _update_details_view(self, submission: Submission) -> None:
        details_panel_static = self.query_one("#details-panel", Static)
        grid = rich.table.Table.grid(expand=True, padding=(0, 1))
        grid.add_column(style="bold blue")
        grid.add_column()
        details_map = self.ID_COLUMNS | self.DETAIL_COLUMNS
        for attr, label in details_map.items():
            value = getattr(submission, attr, None)
            value_str = value.value if hasattr(value, "value") else value
            grid.add_row(
                f"{label}:", str(value_str) if value_str is not None else rich.text.Text("missing", "italic yellow")
            )
        details_panel_static.update(
            rich.panel.Panel(grid, title=f"Details for {submission.id}", border_style="none", padding=(0, 0))
        )

        history_table = self.query_one("#history-table", DataTable)
        history_table.clear()
        if not submission.states:
            return
        sorted_states = sorted(submission.states, key=lambda s: s.timestamp, reverse=False)
        for state_log in sorted_states:
            data_str = json.dumps(state_log.data) if state_log.data else ""
            state = state_log.state
            state_text = rich.text.Text(state.value, style="red" if state == SubmissionStateEnum.ERROR else "default")
            signature_status, comment = _verify_signature(self.public_keys, state_log.author_name, state_log)
            history_table.add_row(
                str(state_log.id),
                state_log.timestamp.isoformat(),
                state_text,
                data_str,
                state_log.author_name,
                signature_status.rich_display(comment),
            )

    def populate_table(self, submissions: list[Submission] | None = None) -> None:  # noqa: C901
        self._table.clear()
        submissions_to_process = list(self.submission_map.values()) if submissions is None else submissions
        sorted_submissions = sorted(
            submissions_to_process,
            key=lambda s: max(s.states, key=lambda st: st.timestamp).timestamp if s.states else datetime.min,
            reverse=True,
        )
        for submission in sorted_submissions:
            row_data = []
            latest_state = max(submission.states, key=lambda s: s.timestamp) if submission.states else None
            for key in self._labels:
                cell = "Error"
                if key in self.ID_COLUMNS or key in self.DETAIL_COLUMNS:
                    value = getattr(submission, key, None)
                    if isinstance(value, bool):
                        cell = rich.text.Text("Yes", style="green") if value else rich.text.Text("No", style="red")
                    else:
                        value_str = value.value if hasattr(value, "value") else value
                        cell = str(value_str) if value_str is not None else "N/A"
                elif key in self.STATE_COLUMNS:
                    if not latest_state:
                        cell = "N/A"
                    elif key == "latest_state":
                        cell = rich.text.Text(
                            latest_state.state.value,
                            style="red" if latest_state.state == SubmissionStateEnum.ERROR else "default",
                        )
                    elif key == "timestamp":
                        cell = latest_state.timestamp.isoformat()
                    elif key == "data_steward":
                        cell = latest_state.author_name
                    elif key == "signature":
                        sig_status, comment = _verify_signature(
                            self.public_keys, latest_state.author_name, latest_state
                        )
                        cell = sig_status.rich_display(comment)
                row_data.append(cell)
            self._table.add_row(*row_data, key=submission.id)

    def _do_sort(self, column_key: ColumnKey | None) -> None:
        if column_key is None:
            return
        if self._sort_by == column_key:
            self._sort_reverse = not self._sort_reverse
        else:
            self._sort_by = column_key
            self._sort_reverse = False

        column_index = self._table.get_column_index(column_key)

        def sort_key(item: tuple) -> str | int | float:
            value = item[column_index]
            if isinstance(value, rich.text.Text):
                return value.plain.casefold()
            if isinstance(value, rich.text.Text) and value.plain in ("Yes", "No"):
                return 1 if value.plain == "Yes" else 0
            return value if value is not None else ""

        self._table.sort(key=sort_key, reverse=self._sort_reverse)
        self._update_header_sorting_indicator()

        if self._table.row_count > 0:
            try:
                submission_id = self._table.get_row_at(self._table.cursor_row)[0]
                if submission_id != self._current_submission_id:
                    self._current_submission_id = submission_id
                    if submission := self.submission_map.get(submission_id):
                        self._update_details_view(submission)
            except IndexError:
                pass

    def _update_header_sorting_indicator(self) -> None:
        indicator = " ▼" if self._sort_reverse else " ▲"
        for key, column in self._table.columns.items():
            original_label = self._labels.get(key.value, "")
            new_label = f"{original_label}{indicator}" if key == self._sort_by else original_label
            column.label = rich.text.Text(new_label)

    def on_data_table_header_selected(self, event: DataTable.HeaderSelected) -> None:
        if event.control.id == "list-table":
            self._do_sort(event.column_key)

    def on_blur(self, event: events.Blur) -> None:
        if event.widget.id == "filter-input" and not event.widget.value:
            event.widget.display = False

    def on_input_changed(self, event: Input.Changed) -> None:
        query = event.value.lower()
        if not query:
            self.populate_table()
            return

        def filterable_values(submission: Submission):
            values = set(submission.model_dump().values())
            if submission.states:
                values |= set(submission.states[0].model_dump().values()) | {submission.states[0].author_name}
            return values

        filtered = [
            submission
            for submission in self.submission_map.values()
            if any(query in str(v).casefold() for v in filterable_values(submission))
        ]
        self.populate_table(filtered)

    def action_toggle_sort_column(self) -> None:
        try:
            column_key = list(self._table.columns.keys())[self._table.cursor_column]
            self._do_sort(column_key)
        except IndexError:
            pass

    def action_focus_filter(self) -> None:
        filter_input = self.query_one("#filter-input")
        filter_input.display = True
        filter_input.focus()
