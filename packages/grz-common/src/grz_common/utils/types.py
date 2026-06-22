import datetime
import enum
import logging

log = logging.getLogger(__name__)


class CaseInsensitiveStrEnum(enum.StrEnum):
    """
    A StrEnum that is case-insensitive for member lookup and comparison with strings.
    """

    @classmethod
    def _missing_(cls, value):
        """
        Override to allow case-insensitive lookup of enum members by value.
        e.g., MyEnum('value') will match MyEnum.VALUE.
        """
        if isinstance(value, str):
            for member in cls:
                if member.value.casefold() == value.casefold():
                    return member
        return None

    def __eq__(self, other):
        """
        Override to allow case-insensitive comparison of enum members by value.
        """
        if isinstance(other, enum.Enum):
            return self is other
        if isinstance(other, str):
            return self.value.casefold() == other.casefold()
        return NotImplemented

    def __hash__(self):
        """
        Override to make hash consistent with eq.
        """
        return hash(self.value.casefold())


class ListableEnum(enum.StrEnum):
    """Mixin for enum classes whose members can be listed."""

    @classmethod
    def list(cls) -> list[str]:
        """Returns a list of enum members."""
        return list(map(lambda c: c.value, cls))


def serialize_datetime_to_iso_z(dt: datetime.datetime) -> str:
    """
    Serializes a datetime object to a canonical ISO 8601 string format with 'Z' for UTC.
    """
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=datetime.UTC)

    if dt.tzinfo != datetime.UTC and dt.utcoffset() != datetime.timedelta(0):
        dt = dt.astimezone(datetime.UTC)

    return dt.isoformat()
