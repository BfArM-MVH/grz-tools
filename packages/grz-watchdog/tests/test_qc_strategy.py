import datetime
import math
import random
import subprocess
import sys
from collections.abc import Generator
from pathlib import Path
from unittest.mock import patch

import pytest
import yaml

# can't directly import from the workflow scripts, so temporarily add to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = PROJECT_ROOT / "workflow/scripts"
sys.path.insert(0, str(SCRIPTS_DIR))
import qc_strategy

sys.path.pop(0)

SUBMITTER_ID = "123456789"


@pytest.fixture(scope="function")
def setup_test_db(tmp_path: Path) -> Generator[tuple[Path, Path], None, None]:
    """
    Create a clean database and config for each test function.
    """
    db_path = tmp_path / "submissions.sqlite"
    config_path = tmp_path / "db_config.yaml"
    known_keys_path = tmp_path / "known_keys"

    db_config_content = {
        "db": {
            "database_url": f"sqlite:///{db_path.resolve()}",
            "author": {
                "name": "test-author",
                "private_key": "-----BEGIN OPENSSH PRIVATE KEY-----\nfoo\n-----END OPENSSH PRIVATE KEY-----",
                "private_key_passphrase": "bar",
            },
            "known_public_keys": str(known_keys_path.resolve()),
        }
    }
    with open(config_path, "w") as f:
        yaml.dump(db_config_content, f)

    with open(known_keys_path, "w") as f:
        f.write("")

    try:
        subprocess.run(
            ["grzctl", "db", "--config-file", str(config_path), "init"],
            check=True,
            capture_output=True,
            text=True,
        )
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        pytest.fail(
            f"Failed to initialize database with `grzctl db init`. Ensure `grzctl` is in your PATH.\nError: {e.stderr or e}"
        )

    yield config_path, db_path


def _add_db_history(
    db_path: Path, submission_id: str, submitter_id: str, states: list[str], initial_timestamp: datetime.datetime
):
    """
    Helper to manually insert a submission and its state history into the submission DB.
    Necessary because we need to control the timestamps of the submission additions / state updates for these tests.
    """
    date_str = initial_timestamp.date().isoformat()

    subprocess.run(
        [
            "sqlite3",
            str(db_path),
            f"INSERT OR IGNORE INTO submissions (id, submitter_id, submission_date) VALUES ('{submission_id}', '{submitter_id}', '{date_str}');",
        ],
        check=True,
    )

    current_timestamp = initial_timestamp
    for state in states:
        current_ts_str = current_timestamp.strftime("%Y-%m-%d %H:%M:%S")
        state_value = state.capitalize()
        subprocess.run(
            [
                "sqlite3",
                str(db_path),
                f"INSERT INTO submission_states (submission_id, state, timestamp, author_name, signature) VALUES ('{submission_id}', '{state_value}', '{current_ts_str}', 'baz', 'dummy-signature');",
            ],
            check=True,
        )
        current_timestamp += datetime.timedelta(minutes=1)


class TestQcStrategy:
    """Tests the QC selection strategy function should_run_qc."""

    @pytest.mark.parametrize(
        "test_date",
        [datetime.date(2025, 10, 15), datetime.date(2025, 1, 1)],
    )
    def test_first_of_month(self, setup_test_db, test_date):
        """Test that the first validated submission of a month is always QCed."""
        config_path, _ = setup_test_db

        with patch("qc_strategy.datetime.date") as mock_date:
            mock_date.today.return_value = test_date
            result = qc_strategy.should_run_qc(
                db_config_path=str(config_path),
                submitter_id=SUBMITTER_ID,
                target_percentage=2.0,
                salt="any_salt",
            )
            assert result is True, "The first submission of the month was not selected for QC."

    def test_deterministic_selection(self, setup_test_db):
        """Test that the random selection is deterministic based on the seed."""
        config_path, db_path = setup_test_db
        test_date = datetime.date(2025, 10, 20)
        salt = "a_very_specific_salt"
        target_percentage = 20.0
        block = math.floor(1 / (target_percentage / 100.0))

        # add one qced submission to db
        _add_db_history(
            db_path,
            "submission-qc-1",
            SUBMITTER_ID,
            ["uploaded", "validated", "qcing"],
            datetime.datetime(2025, 10, 10, 12, 0),
        )
        block_index = 1

        # pre-calculate the target index for the current block
        seed = f"{SUBMITTER_ID}-{test_date.year}-{test_date.month}-{block_index}-{salt}"
        rng = random.Random(seed)
        target_index_for_hit = rng.randint(0, int(block) - 1)

        # add validated submissions to match the target index
        for i in range(target_index_for_hit):
            _add_db_history(
                db_path,
                f"submission-val-hit-{i}",
                SUBMITTER_ID,
                ["uploaded", "validated"],
                datetime.datetime(2025, 10, 11 + i, 10, 0),
            )

        with patch("qc_strategy.datetime.date") as mock_date:
            mock_date.today.return_value = test_date
            selected_for_qc = qc_strategy.should_run_qc(
                db_config_path=str(config_path),
                submitter_id=SUBMITTER_ID,
                target_percentage=target_percentage,
                salt=salt,
            )
            assert selected_for_qc is True, "Should have been selected for qc."
