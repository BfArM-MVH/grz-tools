import calendar
import datetime


def quarter_date_bounds(year: int, quarter: int) -> tuple[datetime.date, datetime.date]:
    """Return the first and last day of a quarter (1-based)."""
    first_month = (quarter - 1) * 3 + 1
    return month_date_bounds(year, first_month)[0], month_date_bounds(year, first_month + 2)[1]


def month_date_bounds(year: int, month: int) -> tuple[datetime.date, datetime.date]:
    """Return the first and last day of a month."""
    _first_weekday, days_in_month = calendar.monthrange(year, month)
    return datetime.date(year, month, 1), datetime.date(year, month, days_in_month)


def date_to_quarter_year(date: datetime.date) -> tuple[int, int]:
    """Return (1-based quarter, year) given a date."""
    return (date.month - 1) // 3 + 1, date.year
