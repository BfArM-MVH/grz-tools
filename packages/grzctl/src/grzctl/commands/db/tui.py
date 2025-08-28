from cryptography.hazmat.primitives.serialization import SSHPublicKeyTypes
from grz_db.models.submission import SubmissionDb
from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import HorizontalScroll
from textual.widgets import DataTable, Footer, TabbedContent, TabPane


class DatabaseBrowser(App):
    """A Textual app to browse a GRZ submission database."""

    BINDINGS = [
        Binding(key="q", action="quit", description="Quit the app"),
    ]

    def __init__(self, database: SubmissionDb, public_keys: dict[str, SSHPublicKeyTypes], **kwargs) -> None:
        super().__init__(**kwargs)
        self._database = database
        self._public_keys = public_keys

    def compose(self) -> ComposeResult:
        with TabbedContent():
            with TabPane("Overview", id="pane-overview"):
                yield HorizontalScroll(id="container-overview-charts")
                yield DataTable(id="table-states-latest")
            with TabPane("Issues", id="pane-issues"):
                yield DataTable(id="table-issues")
            with TabPane("Search", id="pane-search"):
                yield HorizontalScroll(id="container-search-inputs")
                yield DataTable(id="table-search-results")
        yield Footer()
