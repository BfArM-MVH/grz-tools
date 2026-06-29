import datetime

import pytest
from grz_db.common import serialize_datetime_to_iso_z


@pytest.mark.parametrize(
    "dt,description",
    [
        # SQLite: reads back as naive (no tzinfo)
        (datetime.datetime(2024, 1, 1, 12, 0, 0), "SQLite naive datetime"),
        # PostgreSQL: reads back as timezone-aware UTC
        (datetime.datetime(2024, 1, 1, 12, 0, 0, tzinfo=datetime.UTC), "PostgreSQL UTC datetime"),
        # Non-UTC timezone (should be converted to UTC)
        (
            datetime.datetime(2024, 1, 1, 13, 0, 0, tzinfo=datetime.timezone(datetime.timedelta(hours=1))),
            "Non-UTC timezone",
        ),
    ],
)
def test_serialize_datetime_to_iso_z_database_agnostic(dt, description):
    """
    serialize_datetime_to_iso_z must produce identical 'Z'-suffixed output
    regardless of whether the datetime came from SQLite (naive) or PostgreSQL (timezone-aware).
    This ensures signature verification works consistently across database backends.
    """
    result = serialize_datetime_to_iso_z(dt)
    assert result.endswith("Z"), f"{description}: Expected 'Z' suffix, got: {result}"
    assert "+00:00" not in result, f"{description}: Should not contain '+00:00'"
    # All should normalize to the same UTC representation
    assert result == "2024-01-01T12:00:00Z", f"{description}: Expected normalized UTC, got: {result}"
