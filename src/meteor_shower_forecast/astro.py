"""Small self-contained spherical-astronomy helpers.

No third-party dependencies: Julian date / GMST / topocentric alt-az are
standard closed-form formulas (Meeus, *Astronomical Algorithms*), accurate to
well under a degree over centuries -- more than sufficient for judging a
meteor radiant's altitude above the horizon.
"""

from __future__ import annotations

import math
import re
from datetime import date, datetime, timedelta, timezone


def parse_ra(s: str) -> float:
    """Parse a right ascension string like '03h 04m 00s' into decimal degrees.

    Catalog data is stored as decimal degrees; this (and parse_dec) exist
    for scripts/build_shower_catalog.py's curated-source-data table, not
    package runtime code.
    """
    parts = [p for p in re.split(r"[hms\s]+", s.strip()) if p]
    h = float(parts[0])
    m = float(parts[1]) if len(parts) > 1 else 0.0
    sec = float(parts[2]) if len(parts) > 2 else 0.0
    hours = h + m / 60 + sec / 3600
    return hours * 15.0


def parse_dec(s: str) -> float:
    """Parse a declination string like '+58° 00\\' 00\"' into signed decimal degrees."""
    s = s.strip()
    sign = -1 if s.startswith("-") else 1
    parts = [p for p in re.split(r"[°'\"\s]+", s.lstrip("+-")) if p]
    d = float(parts[0])
    m = float(parts[1]) if len(parts) > 1 else 0.0
    sec = float(parts[2]) if len(parts) > 2 else 0.0
    return sign * (d + m / 60 + sec / 3600)


def julian_date(dt_utc: datetime) -> float:
    """Julian date for a timezone-aware (or naive-UTC) datetime."""
    if dt_utc.tzinfo is not None:
        dt_utc = dt_utc.astimezone(timezone.utc).replace(tzinfo=None)
    y, m = dt_utc.year, dt_utc.month
    d = (
        dt_utc.day
        + dt_utc.hour / 24
        + dt_utc.minute / 1440
        + (dt_utc.second + dt_utc.microsecond / 1e6) / 86400
    )
    if m <= 2:
        y -= 1
        m += 12
    a = y // 100
    b = 2 - a + a // 4
    return (
        math.floor(365.25 * (y + 4716))
        + math.floor(30.6001 * (m + 1))
        + d
        + b
        - 1524.5
    )


def gmst_hours(dt_utc: datetime) -> float:
    """Greenwich Mean Sidereal Time, in hours (Meeus ch. 12)."""
    jd = julian_date(dt_utc)
    t = (jd - 2451545.0) / 36525.0
    gmst_deg = (
        280.46061837
        + 360.98564736629 * (jd - 2451545.0)
        + 0.000387933 * t * t
        - t * t * t / 38710000.0
    )
    return (gmst_deg % 360.0) / 15.0


def alt_az(ra_deg: float, dec_deg: float, lat_deg: float, lon_deg: float, dt_utc: datetime) -> tuple[float, float]:
    """Topocentric altitude/azimuth (degrees) of a fixed-RA/Dec object.

    lon_deg is east-positive. Azimuth is measured from north through east.
    """
    lst_hours = (gmst_hours(dt_utc) + lon_deg / 15.0) % 24.0
    hour_angle_deg = (lst_hours * 15.0 - ra_deg) % 360.0

    ha = math.radians(hour_angle_deg)
    dec = math.radians(dec_deg)
    lat = math.radians(lat_deg)

    sin_alt = math.sin(dec) * math.sin(lat) + math.cos(dec) * math.cos(lat) * math.cos(ha)
    alt = math.asin(max(-1.0, min(1.0, sin_alt)))

    cos_az = (math.sin(dec) - math.sin(alt) * math.sin(lat)) / (math.cos(alt) * math.cos(lat) + 1e-12)
    az = math.acos(max(-1.0, min(1.0, cos_az)))
    if math.sin(ha) > 0:
        az = 2 * math.pi - az

    return math.degrees(alt), math.degrees(az)


def moon_illumination_pct(dt_utc: datetime) -> float:
    """Approximate Moon illuminated fraction (0-100%) at dt_utc.

    Derived from the Moon's mean elongation from the Sun (Meeus,
    *Astronomical Algorithms*, ch. 49's D term) via k = (1 - cos(D)) / 2.
    This is a low-precision approximation -- it omits the periodic
    correction terms (lunar/solar anomaly, evection, etc.) a full lunar
    ephemeris would include, so it can be off by several percent -- but
    needs no ephemeris data and is accurate enough to judge how much
    moonlight will degrade meteor visibility on a given night.
    """
    jd = julian_date(dt_utc)
    t = (jd - 2451545.0) / 36525.0
    d = (
        297.8501921
        + 445267.1114034 * t
        - 0.0018819 * t * t
        + t**3 / 545868.0
        - t**4 / 113065000.0
    ) % 360.0
    return (1 - math.cos(math.radians(d))) / 2 * 100.0


def solar_longitude(dt_utc: datetime) -> float:
    """Apparent geocentric ecliptic longitude of the Sun (degrees, 0-360),
    equinox of date -- the "solar longitude" (lambda-sun) convention used
    throughout meteor-shower literature (IAU MDC, IMO) to specify activity
    periods independent of calendar leap-year drift.

    Low-precision formula (Meeus, *Astronomical Algorithms*, ch. 25),
    accurate to about 0.01 degrees (a few minutes of time). Used both by
    scripts/build_shower_catalog.py (solar-longitude -> calendar-date) and
    at runtime by sun_ra_dec (-> sun_altitude, for excluding daylight from
    the best-observing-time search in forecast.py).
    """
    jd = julian_date(dt_utc)
    t = (jd - 2451545.0) / 36525.0
    l0 = 280.46646 + 36000.76983 * t + 0.0003032 * t * t
    m = math.radians(357.52911 + 35999.05029 * t - 0.0001537 * t * t)
    c = (
        (1.914602 - 0.004817 * t - 0.000014 * t * t) * math.sin(m)
        + (0.019993 - 0.000101 * t) * math.sin(2 * m)
        + 0.000289 * math.sin(3 * m)
    )
    true_longitude = l0 + c
    omega = math.radians(125.04 - 1934.136 * t)
    apparent = true_longitude - 0.00569 - 0.00478 * math.sin(omega)
    return apparent % 360.0


def sun_ra_dec(dt_utc: datetime) -> tuple[float, float]:
    """Apparent geocentric equatorial coordinates of the Sun (RA, Dec;
    decimal degrees, J2000-ish equinox of date), converted from
    solar_longitude via the standard ecliptic->equatorial rotation by
    Earth's mean obliquity (Meeus ch. 25).
    """
    jd = julian_date(dt_utc)
    t = (jd - 2451545.0) / 36525.0
    epsilon = math.radians(23.439291 - 0.0130042 * t)
    lam = math.radians(solar_longitude(dt_utc))
    ra = math.degrees(math.atan2(math.cos(epsilon) * math.sin(lam), math.cos(lam))) % 360.0
    dec = math.degrees(math.asin(math.sin(epsilon) * math.sin(lam)))
    return ra, dec


def sun_altitude(lat_deg: float, lon_deg: float, dt_utc: datetime) -> float:
    """Topocentric altitude of the Sun (degrees) -- used to keep the
    best-observing-time search (forecast.py) from ever landing in
    daylight.
    """
    ra_deg, dec_deg = sun_ra_dec(dt_utc)
    alt, _ = alt_az(ra_deg, dec_deg, lat_deg, lon_deg, dt_utc)
    return alt


def date_for_solar_longitude(target_deg: float, year: int) -> datetime:
    """Inverse of solar_longitude: the UTC datetime nearest `year`-01-01
    at which the Sun reaches `target_deg` of apparent ecliptic longitude.

    Newton's method using the (near-constant, ~0.9856 deg/day) rate of
    change of solar longitude; converges to sub-second precision in a
    handful of iterations since the rate varies by only ~3% over a year.
    """
    target_deg = target_deg % 360.0
    t = datetime(year, 1, 1, tzinfo=timezone.utc) + timedelta(days=(target_deg / 360.0) * 365.2422)
    for _ in range(8):
        lam = solar_longitude(t)
        diff = (target_deg - lam + 180) % 360 - 180
        t += timedelta(days=diff / 0.9856)
    return t


def resolve_nearest_peak(peak_month: int, peak_day: int, peak_hour_utc: float, reference: datetime) -> tuple[datetime, float]:
    """Find the calendar instance of an annual peak (month/day/hour-UTC)
    nearest to `reference`, handling year-boundary showers (e.g. a shower
    peaking Jan 3 considered from a late-December reference date).

    Returns (peak_datetime_utc, signed_days_from_peak) where a negative
    delta means the peak is still ahead of `reference`.
    """
    if reference.tzinfo is None:
        reference = reference.replace(tzinfo=timezone.utc)
    ref_date = reference.astimezone(timezone.utc).date()

    best: tuple[datetime, float] | None = None
    for year_offset in (-1, 0, 1):
        try:
            peak_date = date(ref_date.year + year_offset, peak_month, peak_day)
        except ValueError:
            continue
        peak_dt = datetime(
            peak_date.year, peak_date.month, peak_date.day, tzinfo=timezone.utc
        ) + timedelta(hours=peak_hour_utc)
        delta_days = (reference - peak_dt).total_seconds() / 86400.0
        if best is None or abs(delta_days) < abs(best[1]):
            best = (peak_dt, delta_days)
    assert best is not None
    return best


def next_peak(peak_month: int, peak_day: int, peak_hour_utc: float, reference: datetime) -> tuple[datetime, float]:
    """Soonest peak instance at or after `reference` (unlike
    resolve_nearest_peak, never returns one that has already passed).

    Returns (peak_datetime_utc, days_until_peak), days_until_peak >= 0.
    """
    if reference.tzinfo is None:
        reference = reference.replace(tzinfo=timezone.utc)
    ref_date = reference.astimezone(timezone.utc).date()

    best: tuple[datetime, float] | None = None
    for year_offset in (0, 1):
        try:
            peak_date = date(ref_date.year + year_offset, peak_month, peak_day)
        except ValueError:
            continue
        peak_dt = datetime(
            peak_date.year, peak_date.month, peak_date.day, tzinfo=timezone.utc
        ) + timedelta(hours=peak_hour_utc)
        if peak_dt < reference:
            continue
        days_until = (peak_dt - reference).total_seconds() / 86400.0
        if best is None or days_until < best[1]:
            best = (peak_dt, days_until)
    assert best is not None
    return best
