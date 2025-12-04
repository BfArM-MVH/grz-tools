import calendar
import datetime


def quarter_date_bounds(year: int, quarter: int) -> tuple[datetime.date, datetime.date]:
    quarter_start_date = datetime.date(year=year, month=((quarter - 1) * 3) + 1, day=1)
    quarter_end_month = quarter_start_date.month + 2
    _quarter_end_month_first_weekday, days_in_quarter_end_month = calendar.monthrange(year, quarter_end_month)
    quarter_end_date = quarter_start_date.replace(month=quarter_end_month, day=days_in_quarter_end_month)

    return quarter_start_date, quarter_end_date


def date_to_quarter_year(date: datetime.date) -> tuple[int, int]:
    """Return (1-based quarter, year) given a date."""
    return (date.month - 1) // 3 + 1, date.year
