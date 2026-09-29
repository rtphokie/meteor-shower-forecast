from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone

from . import localtime
from .catalog import SHOWERS
from .forecast import RateEstimate, estimate_all_showers, estimate_shower_rate, next_shower


def _parse_when(value: str) -> datetime:
    dt = datetime.fromisoformat(value)
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    return dt


def _bortle_str(e: RateEstimate) -> str:
    tag = " (estimated from coordinates)" if e.site_bortle_estimated else ""
    return f"Bortle {e.site_bortle_class} ({e.site_bortle_desc}, SQM {e.site_sqm:.2f}){tag}"


def _moon_str(e: RateEstimate) -> str:
    tag = " (estimated for this date/time)" if e.moon_illumination_estimated else ""
    return f"{e.moon_illumination_pct:.0f}% illuminated{tag}"


def _fmt_local(dt: datetime, tz, *, with_weekday: bool = True) -> str:
    """Human-friendly local time, e.g. 'Wed, Oct 21, 2026, 2:54 AM EDT'.

    Uses the tz abbreviation (EDT/PST/...) rather than the IANA zone name.
    """
    local = dt.astimezone(tz)
    hour12 = local.hour % 12 or 12
    weekday = f"{local.strftime('%a')}, " if with_weekday else ""
    return f"{weekday}{local.strftime('%b')} {local.day}, {local.year}, {hour12}:{local.strftime('%M %p')} {local.tzname()}"


def _print_estimate(e: RateEstimate, tz) -> None:
    print(f"{e.shower} ({e.iau_code})")
    print(f"  parent body:         {e.parent_body or 'unknown'}")
    when_label = "best:" if e.when_is_best_time else "time:"
    print(f"  {when_label:<21}{_fmt_local(e.when_utc, tz)}")
    print(f"  peak:                {_fmt_local(e.peak_dt_utc, tz)}")
    print(f"  days from peak:      {e.days_from_peak:+.1f}")
    if e.zhr_effective is not None:
        print(f"  effective ZHR:       {e.zhr_effective:.1f}")
    print(f"  radiant altitude:    {e.radiant_alt_deg:.1f}° (az {e.radiant_az_deg:.1f}°)")
    print(f"  sky brightness:      {_bortle_str(e)}")
    print(f"  moon:                {_moon_str(e)}")
    if e.limiting_magnitude is not None:
        print(f"  limiting magnitude:  {e.limiting_magnitude:.2f}")
        print(f"  faint-meteor factor: {e.lm_factor:.3f}")
    if e.estimated_rate_per_hour is not None:
        print(f"  ~ estimated rate:    {e.estimated_rate_per_hour:.1f} meteors/hour")
    else:
        print("  ~ estimated rate:    no data")
    for note in e.notes:
        print(f"  note: {note}")


_PEAK_COL_WIDTH = 30
_PARENT_COL_WIDTH = 26
_SHOWER_COL_WIDTH = 26


def _print_table(estimates: list[RateEstimate], tz) -> None:
    if not estimates:
        print("No showers are meaningfully active at this location and time.")
        return
    print(f"Sky brightness: {_bortle_str(estimates[0])}")
    print()
    header = (
        f"{'Shower':<{_SHOWER_COL_WIDTH}}{'Code':<6}{'Rate/hr':>9}{'ZHR eff':>9}{'Radiant alt':>13}{'Moon':>6}  "
        f"{'Peak':>{_PEAK_COL_WIDTH}}  {'Parent body':<{_PARENT_COL_WIDTH}}"
    )
    print(header)
    print("-" * len(header))
    for e in estimates:
        peak_str = _fmt_local(e.peak_dt_utc, tz, with_weekday=False)
        parent = (e.parent_body or "unknown")[:_PARENT_COL_WIDTH]
        print(
            f"{e.shower:<{_SHOWER_COL_WIDTH}}{e.iau_code:<6}{e.estimated_rate_per_hour:>9.1f}{e.zhr_effective:>9.1f}"
            f"{e.radiant_alt_deg:>12.1f}°{e.moon_illumination_pct:>5.0f}%  {peak_str:>{_PEAK_COL_WIDTH}}"
            f"  {parent:<{_PARENT_COL_WIDTH}}"
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="meteor-shower-forecast",
        description=(
            "Estimate the realistic naked-eye meteor rate (meteors/hour) for a "
            "location, accounting for the shower's day-decayed ZHR, the "
            "radiant's altitude above the horizon, and local sky brightness "
            "(Bortle class and, optionally, moonlight)."
        ),
    )
    parser.add_argument("--lat", type=float, help="Observer latitude, degrees (north positive)")
    parser.add_argument("--lon", type=float, help="Observer longitude, degrees (east positive)")
    when_group = parser.add_mutually_exclusive_group()
    when_group.add_argument(
        "--when",
        type=_parse_when,
        default=None,
        help=(
            "Observation time as ISO-8601 (default: each shower's own next best "
            "observing opportunity, not now). Naive times are treated as UTC."
        ),
    )
    when_group.add_argument(
        "--year",
        type=int,
        default=None,
        metavar="YYYY",
        help=(
            "Look at that year's occurrence instead of the next upcoming one (past "
            "years work too). Still searches for the best time within it, unlike "
            "--when, which pins one exact instant."
        ),
    )
    sky_group = parser.add_mutually_exclusive_group()
    sky_group.add_argument(
        "--bortle",
        type=int,
        choices=range(1, 10),
        metavar="1-9",
        help="Bortle dark-sky class of the site (default: estimated from --lat/--lon)",
    )
    sky_group.add_argument("--sqm", type=float, help="Sky brightness directly, in mag/arcsec^2, instead of --bortle")
    parser.add_argument(
        "--moon-illumination",
        type=float,
        default=None,
        metavar="PCT",
        help=(
            "Moon illumination percentage (0-100) to factor into sky brightness "
            "(default: computed automatically from the observation date/time; pass 0 to ignore moonlight)"
        ),
    )
    parser.add_argument(
        "--shower",
        type=str,
        default=None,
        help="Estimate a single named shower, or 'all' for every active shower (default: the next upcoming shower)",
    )
    parser.add_argument("--list-showers", action="store_true", help="List known showers and exit")
    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.list_showers:
        for s in SHOWERS:
            if s.has_rate_data:
                zhr_str = f"ZHR {s.peak_zhr:.0f}"
            elif s.peak_zhr is not None:
                zhr_str = f"ZHR {s.peak_zhr:.0f}, but no rate estimate (no decay-rate data)"
            else:
                zhr_str = "no rate data (variable/outburst shower)"
            parent = f", parent {s.parent_body}" if s.parent_body else ""
            print(f"{s.iau_code:<5} {s.name:<32} peak {s.peak_month:02d}-{s.peak_day:02d}  {zhr_str}{parent}")
        return

    missing = [
        flag
        for flag, val in (("--lat", args.lat), ("--lon", args.lon))
        if val is None
    ]
    if missing:
        parser.error(f"the following arguments are required: {', '.join(missing)}")

    year_start = datetime(args.year, 1, 1, tzinfo=timezone.utc) if args.year is not None else None
    kwargs = dict(
        bortle=args.bortle, sqm=args.sqm, moon_illumination_pct=args.moon_illumination, reference=year_start
    )
    tz, _ = localtime.zone_at(args.lat, args.lon)

    try:
        if args.shower and args.shower.lower() != "all":
            estimate = estimate_shower_rate(args.shower, args.lat, args.lon, args.when, **kwargs)
            _print_estimate(estimate, tz)
        elif args.shower:
            estimates = estimate_all_showers(args.lat, args.lon, args.when, **kwargs)
            _print_table(estimates, tz)
        else:
            shower = next_shower(args.when or year_start)
            estimate = estimate_shower_rate(shower, args.lat, args.lon, args.when, **kwargs)
            _print_estimate(estimate, tz)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
    except FileNotFoundError as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()
