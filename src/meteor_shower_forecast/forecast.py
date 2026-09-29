"""Realistic per-location meteor rate estimates.

Combines, for a given shower, coordinates, and time:

  1. the shower's day-decayed effective ZHR (IMO/NASA-MEO-style asymmetric
     log-linear decay around peak date),
  2. the radiant's actual altitude above the observer's local horizon at
     that specific place and moment -- topocentric, so it changes
     continuously through the night as Earth rotates the fixed celestial
     radiant point under the observer (a radiant at the zenith contributes
     its full ZHR; one low on the horizon is heavily foreshortened -- ZHR
     is defined *at* zenith). The radiant itself can be up to ~10 degrees
     below that local horizon and meteors are still geometrically visible,
     because the meteors themselves ablate ~100km up, well above the
     radiant's own point on the horizon -- see HORIZON_DIP_DEG. (This is
     independent of "earthgrazers", the unrelated term for meteors that
     enter the atmosphere at a shallow trajectory angle and skip back out.), and
  3. how much the local sky brightness (from the site's Bortle class, and
     optionally lunar illumination) degrades the naked-eye limiting
     magnitude, which disproportionately hides a shower's fainter meteors.

The combination `zhr_effective * sin(radiant_alt) * lm_factor` and its
constituent formulas are adapted from the DarkHours project
(mbeher2200/DarkHours, darkhours/predictor.py and darkhours/targets.py),
MIT licensed. DarkHours has no PyPI package; the formulas are reimplemented
here in a small, dependency-free form rather than vendoring the (much
larger, service-oriented) upstream codebase.

If neither `bortle` nor `sqm` is given, sky brightness is instead derived
from (lat, lon) via the `bortlefinder` package (extracted from this
project's own earlier coordinate-based Bortle estimation and published
standalone to PyPI) -- see its README for the data source and accuracy
caveat.

The catalog (catalog.py) covers the ~38 showers on the IMO Working List,
but real decay-rate data -- required, alongside ZHR/population index, for
an actual rate estimate -- only exists for a small, well-observed subset
(see ShowerDef.has_rate_data). For the rest,
`zhr_effective`/`limiting_magnitude`/`lm_factor`/`estimated_rate_per_hour`
are all None; radiant geometry, peak timing, and parent body are still
reported.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone

import bortlefinder

from . import astro
from .catalog import SHOWERS, ShowerDef, find_shower
from .skybrightness import moon_delta_mag

LOW_RADIANT_ALT_DEG = 25.0  # below this, local rate collapses below ~45% of zhr_effective
ZHR_DECAY_FLOOR = 2.0  # meteors/hr -- below this, indistinguishable from sporadic background

# A radiant below the observer's local horizon isn't necessarily invisible:
# meteors ablate ~100km up, so the luminous trail can still clear that
# (topocentric, moment-by-moment) horizon even when the radiant direction
# itself can't.
#
# HORIZON_DIP_DEG is an order-of-magnitude ESTIMATE of how far below the
# horizon that can still hold, not a value taken from a meteor-science
# source -- no published near-horizon visibility model was found. It
# combines two separately well-established facts: the standard geodesic
# "dip of the horizon" formula for an object at height h above a sphere of
# radius R (dip = arccos(R / (R + h)); the same formula used for how far
# off a ship's mast is visible over the sea horizon), applied with
# h = 100km, the middle of the commonly cited ~80-120km meteor ablation
# altitude range. Treat this as a plausible geometric ceiling, not a
# validated cutoff -- meteors could in practice stop being visible well
# before it, and the exact number moves with the (variable) ablation
# altitude assumed.
_EARTH_RADIUS_KM = 6371.0
_METEOR_LAYER_KM = 100.0  # middle of the commonly cited ~80-120km ablation range
HORIZON_DIP_DEG = math.degrees(math.acos(_EARTH_RADIUS_KM / (_EARTH_RADIUS_KM + _METEOR_LAYER_KM)))


@dataclass
class RateEstimate:
    shower: str
    iau_code: str
    parent_body: str | None
    when_utc: datetime
    when_is_best_time: bool  # True if when_utc was found by the observer-specific best-time search, not requested/derived
    rate_at_catalog_peak: float | None  # what the same site/sky would show AT peak_dt_utc; only set when when_is_best_time
    peak_dt_utc: datetime
    days_from_peak: float
    radiant_alt_deg: float
    radiant_az_deg: float
    site_sqm: float
    site_bortle_class: int
    site_bortle_desc: str
    site_bortle_estimated: bool
    moon_illumination_pct: float
    moon_illumination_estimated: bool
    # None throughout this group means "no published ZHR/decay-rate data
    # for this shower" (see catalog.py) -- not "computed to be zero".
    zhr_effective: float | None
    limiting_magnitude: float | None
    lm_factor: float | None
    estimated_rate_per_hour: float | None
    notes: list[str]


def effective_zhr(peak_zhr: float, days_from_peak: float, b_rise: float, b_decline: float) -> float:
    """IMO/NASA-MEO-style asymmetric log-linear decay: ZHR(t) = peak_zhr * 10^(-B*|t|).

    b_rise applies before peak, b_decline at/after peak -- real showers
    typically rise faster approaching peak than they decline afterward.
    """
    if peak_zhr <= 0:
        return 0.0
    b = b_rise if days_from_peak < 0 else b_decline
    return peak_zhr * (10 ** (-b * abs(days_from_peak)))


def meaningful_activity_half_window(peak_zhr: float, b_rise: float, b_decline: float, floor: float = ZHR_DECAY_FLOOR) -> float:
    """Day-offset beyond which decayed ZHR drops below `floor` on either
    side -- solves peak_zhr * 10^(-B*t) = floor for t, using the shallower
    (slower-decaying) of b_rise/b_decline so the window is generous in
    both directions. Bounds how far from the catalog peak it's worth
    searching for a better observing time.
    """
    b_min = min(b_rise, b_decline)
    if peak_zhr <= floor or b_min <= 0:
        return 0.0
    return math.log10(peak_zhr / floor) / b_min


_SUN_ALT_DARK_DEG = -18.0  # astronomical twilight -- matches the fully-dark-sky assumption baked into Bortle/SQM
_SEARCH_STEP = timedelta(minutes=15)
_SEARCH_MAX_HALF_WINDOW_DAYS = 30.0
_SEARCH_MIN_HALF_WINDOW_DAYS = 0.5  # always search at least +-12h, even for a fast-decaying shower


def _raw_rate(
    shower: ShowerDef, lat: float, lon: float, when_utc: datetime, site_sqm: float, moon_illumination_pct: float
) -> float:
    """The core rate formula, stripped of notes/rounding -- shared by the
    detailed estimate below and the best-observing-time search, so the two
    can never drift apart.
    """
    _, days_from_peak = astro.resolve_nearest_peak(
        shower.peak_month, shower.peak_day, shower.peak_hour_utc, when_utc
    )
    zhr_eff = effective_zhr(shower.peak_zhr, days_from_peak, shower.b_rise, shower.b_decline)
    radiant_alt, _ = astro.alt_az(shower.radiant_ra_deg, shower.radiant_dec_deg, lat, lon, when_utc)
    if radiant_alt <= -HORIZON_DIP_DEG or zhr_eff <= 0:
        return 0.0
    delta_mag = 0.0
    if moon_illumination_pct > 0:
        delta_mag = moon_delta_mag(moon_illumination_pct, site_sqm, target_alt_deg=radiant_alt)
    limiting_mag = bortlefinder.nelm_from_sqm(site_sqm - delta_mag)
    lm_factor = min(1.0, shower.population_index ** (limiting_mag - 6.5))
    return zhr_eff * math.sin(math.radians(abs(radiant_alt))) * lm_factor


def _search_best_observing_time(
    shower: ShowerDef,
    lat: float,
    lon: float,
    site_sqm: float,
    moon_illumination_pct: float | None,
    reference: datetime,
) -> tuple[datetime, bool]:
    """Grid-search the shower's active window for the datetime that
    maximizes the estimated rate at (lat, lon).

    The catalog peak is when Earth crosses the shower's dust trail most
    densely -- a global, radiant-independent property. It says nothing
    about when the radiant is highest in *this* observer's sky, which
    cycles roughly every 24h; a few hours off the catalog peak, with the
    radiant higher up, commonly beats the peak itself. Restricted to
    astronomical night (Sun below _SUN_ALT_DARK_DEG) so a daytime
    altitude maximum is never selected; never searches before `reference`.

    Returns (best_datetime, found_dark_time). When found_dark_time is
    False, the search window never reached astronomical darkness at this
    location (e.g. high-latitude summer) and best_datetime is instead the
    least-bright moment found, as a fallback.
    """
    peak_dt, _ = astro.next_peak(shower.peak_month, shower.peak_day, shower.peak_hour_utc, reference)
    half_window_days = min(
        _SEARCH_MAX_HALF_WINDOW_DAYS,
        max(_SEARCH_MIN_HALF_WINDOW_DAYS, meaningful_activity_half_window(shower.peak_zhr, shower.b_rise, shower.b_decline)),
    )
    start = max(reference, peak_dt - timedelta(days=half_window_days))
    end = peak_dt + timedelta(days=half_window_days)

    best_dt, best_rate = start, -1.0
    fallback_dt, fallback_sun_alt = start, float("inf")
    found_dark = False

    t = start
    while t <= end:
        sun_alt = astro.sun_altitude(lat, lon, t)
        if sun_alt < _SUN_ALT_DARK_DEG:
            found_dark = True
            moon_pct = moon_illumination_pct if moon_illumination_pct is not None else astro.moon_illumination_pct(t)
            rate = _raw_rate(shower, lat, lon, t, site_sqm, moon_pct)
            if rate > best_rate:
                best_rate, best_dt = rate, t
        elif sun_alt < fallback_sun_alt:
            fallback_sun_alt, fallback_dt = sun_alt, t
        t += _SEARCH_STEP

    return (best_dt, True) if found_dark else (fallback_dt, False)


def estimate_shower_rate(
    shower: ShowerDef | str,
    lat: float,
    lon: float,
    when_utc: datetime | None = None,
    *,
    bortle: int | None = None,
    sqm: float | None = None,
    moon_illumination_pct: float | None = None,
    reference: datetime | None = None,
) -> RateEstimate:
    """Estimate the realistic per-hour meteor rate a naked-eye observer at
    (lat, lon) would see from a specific shower at a specific time.

    If `when_utc` is omitted, it defaults to *this observer's* best
    upcoming observing opportunity -- not the shower's catalog peak, which
    is a global, radiant-independent instant (when Earth crosses the dust
    trail most densely) that says nothing about whether the radiant is
    even above this observer's horizon then. See
    _search_best_observing_time for how that's found. `reference` is the
    anchor that search starts from (soonest occurrence at or after it) --
    it defaults to now, but e.g. Jan 1 of a given year finds/searches that
    year's occurrence instead (past years work too: nothing here assumes
    `reference` is in the future). Ignored if `when_utc` is given.

    At most one of `bortle` (1-9) or `sqm` (mag/arcsec^2) may be given. If
    neither is given, sky brightness is derived from (lat, lon) via
    `bortlefinder.estimate` -- this requires having run
    `bortlefinder-fetch-grid` once first (see README).

    If `moon_illumination_pct` is omitted, it is computed automatically
    from `when_utc` (see astro.moon_illumination_pct) rather than ignored;
    pass 0 explicitly to disregard moonlight entirely.
    """
    if bortle is not None and sqm is not None:
        raise ValueError("Pass at most one of bortle= or sqm=")

    notes: list[str] = []
    site_bortle_estimated = False
    if bortle is not None:
        site_sqm = bortlefinder.bortle_to_sqm(bortle)
    elif sqm is not None:
        site_sqm = sqm
    else:
        site_sqm = bortlefinder.estimate(lat, lon).sqm
        site_bortle_estimated = True

    site_bortle_class, site_bortle_desc = bortlefinder.sqm_to_bortle(site_sqm)

    if isinstance(shower, str):
        shower = find_shower(shower)

    when_is_best_time = False
    if when_utc is None:
        reference = reference or datetime.now(timezone.utc)
        if reference.tzinfo is None:
            reference = reference.replace(tzinfo=timezone.utc)
        if not shower.has_rate_data:
            # No ZHR/decay data to optimize against -- nothing to search.
            when_utc, _ = astro.next_peak(shower.peak_month, shower.peak_day, shower.peak_hour_utc, reference)
        else:
            when_utc, found_dark = _search_best_observing_time(
                shower, lat, lon, site_sqm, moon_illumination_pct, reference
            )
            when_is_best_time = True
            if not found_dark:
                notes.append(
                    "search window never reached astronomical darkness at this "
                    "location -- showing the least-bright available time instead"
                )
    elif when_utc.tzinfo is None:
        when_utc = when_utc.replace(tzinfo=timezone.utc)

    moon_illumination_estimated = moon_illumination_pct is None
    if moon_illumination_estimated:
        moon_illumination_pct = astro.moon_illumination_pct(when_utc)

    peak_dt, days_from_peak = astro.resolve_nearest_peak(
        shower.peak_month, shower.peak_day, shower.peak_hour_utc, when_utc
    )
    radiant_alt, radiant_az = astro.alt_az(
        shower.radiant_ra_deg, shower.radiant_dec_deg, lat, lon, when_utc
    )

    # What this same site/sky would show AT the catalog peak, for
    # comparison against the searched best time -- computed with the
    # peak's own moon illumination (auto) or the user's fixed override.
    rate_at_peak = None
    if when_is_best_time:
        peak_moon_pct = moon_illumination_pct if not moon_illumination_estimated else astro.moon_illumination_pct(peak_dt)
        rate_at_peak = round(_raw_rate(shower, lat, lon, peak_dt, site_sqm, peak_moon_pct), 1)

    def _bail(**rate_fields) -> RateEstimate:
        return RateEstimate(
            shower=shower.name,
            iau_code=shower.iau_code,
            parent_body=shower.parent_body,
            when_utc=when_utc,
            when_is_best_time=when_is_best_time,
            rate_at_catalog_peak=rate_at_peak,
            peak_dt_utc=peak_dt,
            days_from_peak=days_from_peak,
            radiant_alt_deg=round(radiant_alt, 1),
            radiant_az_deg=round(radiant_az, 1),
            site_sqm=round(site_sqm, 2),
            site_bortle_class=site_bortle_class,
            site_bortle_desc=site_bortle_desc,
            site_bortle_estimated=site_bortle_estimated,
            moon_illumination_pct=round(moon_illumination_pct, 1),
            moon_illumination_estimated=moon_illumination_estimated,
            notes=notes,
            **rate_fields,
        )

    if not shower.has_rate_data:
        missing = "ZHR" if shower.peak_zhr is None else "population-index/decay-rate"
        notes.append(
            f"no published {missing} data for this shower -- radiant geometry "
            "only, no rate estimate"
        )
        return _bail(zhr_effective=None, limiting_magnitude=None, lm_factor=None, estimated_rate_per_hour=None)

    zhr_eff = round(effective_zhr(shower.peak_zhr, days_from_peak, shower.b_rise, shower.b_decline), 1)

    if radiant_alt <= -HORIZON_DIP_DEG:
        notes.append(
            f"radiant is {abs(radiant_alt):.0f}° below the local horizon -- beyond the "
            f"~{HORIZON_DIP_DEG:.0f}° estimated geometric limit (see HORIZON_DIP_DEG), "
            "no meteors expected"
        )
        return _bail(zhr_effective=zhr_eff, limiting_magnitude=None, lm_factor=None, estimated_rate_per_hour=0.0)

    if zhr_eff < ZHR_DECAY_FLOOR:
        notes.append("shower not meaningfully active on this date")
    if radiant_alt < 0:
        notes.append(
            f"radiant is {abs(radiant_alt):.0f}° below the local horizon -- rate estimate "
            "in this band is a rough geometric extrapolation, not a validated model"
        )
    elif radiant_alt < LOW_RADIANT_ALT_DEG:
        notes.append(f"low radiant altitude ({radiant_alt:.0f}°) sharply foreshortens the rate")

    delta_mag = 0.0
    if moon_illumination_pct > 0:
        delta_mag = moon_delta_mag(moon_illumination_pct, site_sqm, target_alt_deg=radiant_alt)
        if delta_mag > 0.1:
            notes.append(f"moonlight brightens the sky by {delta_mag:.2f} mag/arcsec²")

    limiting_mag = bortlefinder.nelm_from_sqm(site_sqm - delta_mag)
    r = shower.population_index
    lm_factor = min(1.0, r ** (limiting_mag - 6.5))

    # abs(): within HORIZON_DIP_DEG below the horizon, mirror the same small
    # foreshortening factor a radiant that far *above* the horizon would
    # give, rather than an implausible hard jump from "near zero" to
    # "exactly zero" at the moment the radiant dips under the horizon.
    base_rate = zhr_eff * math.sin(math.radians(abs(radiant_alt)))
    local_rate = base_rate * lm_factor
    est_rate = round(local_rate, 1)

    if rate_at_peak is not None and est_rate - rate_at_peak > 0.05:
        notes.append(
            f"{est_rate - rate_at_peak:+.1f} more meteors/hour than observing at the "
            "shower's catalog peak from this location"
        )

    return _bail(
        zhr_effective=zhr_eff,
        limiting_magnitude=round(limiting_mag, 2),
        lm_factor=round(lm_factor, 3),
        estimated_rate_per_hour=est_rate,
    )


def next_shower(when_utc: datetime | None = None) -> ShowerDef:
    """The catalog shower with the soonest upcoming peak at or after `when_utc`,
    restricted to showers with real rate data (ShowerDef.has_rate_data) --
    there's no point defaulting to a shower this package cannot produce a
    rate estimate for. Pass a specific ShowerDef/name to
    estimate_shower_rate directly to look up any of the other IMO-listed
    showers.
    """
    when_utc = when_utc or datetime.now(timezone.utc)
    if when_utc.tzinfo is None:
        when_utc = when_utc.replace(tzinfo=timezone.utc)
    rated = [s for s in SHOWERS if s.has_rate_data]
    return min(
        rated,
        key=lambda s: astro.next_peak(s.peak_month, s.peak_day, s.peak_hour_utc, when_utc)[1],
    )


def estimate_all_showers(
    lat: float,
    lon: float,
    when_utc: datetime | None = None,
    *,
    bortle: int | None = None,
    sqm: float | None = None,
    moon_illumination_pct: float | None = None,
    reference: datetime | None = None,
    active_only: bool = True,
) -> list[RateEstimate]:
    """Estimate every catalog shower's rate, sorted chronologically by peak.

    Showers with no published ZHR data (see catalog.py) are always
    excluded -- there's nothing to rank them by. With active_only=True
    (default), showers whose decayed ZHR has fallen below the sporadic-
    background floor are also omitted. See estimate_shower_rate for
    `reference` (ignored if `when_utc` is given).
    """
    estimates = [
        estimate_shower_rate(
            s, lat, lon, when_utc, bortle=bortle, sqm=sqm,
            moon_illumination_pct=moon_illumination_pct, reference=reference,
        )
        for s in SHOWERS
        if s.has_rate_data
    ]
    if active_only:
        estimates = [e for e in estimates if e.zhr_effective >= ZHR_DECAY_FLOOR]
    return sorted(estimates, key=lambda e: e.peak_dt_utc)
