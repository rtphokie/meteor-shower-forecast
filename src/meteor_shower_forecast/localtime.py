"""Coordinate -> IANA timezone, for displaying results in local time.

Internally every datetime in this package is UTC; this is purely a
presentation concern for the CLI.
"""

from __future__ import annotations

from datetime import timezone as dt_timezone
from functools import lru_cache
from zoneinfo import ZoneInfo


@lru_cache(maxsize=1)
def _finder():
    from timezonefinder import TimezoneFinder

    return TimezoneFinder()


def zone_at(lat: float, lon: float) -> tuple[ZoneInfo | dt_timezone, str]:
    """Return (tzinfo, name) for the given coordinates.

    Falls back to UTC for the rare coordinates outside every timezone
    polygon (e.g. open ocean far from any boundary).
    """
    name = _finder().timezone_at(lat=lat, lng=lon)
    if name is None:
        return dt_timezone.utc, "UTC"
    return ZoneInfo(name), name
