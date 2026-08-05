from datetime import UTC, datetime, tzinfo

from pydantic import BaseModel, ConfigDict
from pydantic.alias_generators import to_camel


class StrictBaseModel(BaseModel):
    model_config = ConfigDict(
        extra="forbid",
        validate_assignment=True,
        use_enum_values=True,
        alias_generator=to_camel,
    )


def as_aware_datetime(value: datetime, zone_default: tzinfo = UTC) -> datetime:
    """
    Attach ``zone_default`` to a naive datetime, leaving an already aware one untouched.

    Naive datetimes are treated as UTC rather than local time so that comparisons stay consistent
    across environments, no matter whether the value was parsed from a string without an offset or
    constructed programmatically.

    :param value: datetime to normalize.
    :param zone_default: tzinfo to attach when ``value`` carries none.
    :returns: timezone-aware datetime.
    """
    return value.replace(tzinfo=zone_default) if value.tzinfo is None else value
