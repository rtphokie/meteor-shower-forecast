#!/usr/bin/env python3
"""Build src/meteor_shower_forecast/data/shower_catalog.csv from the IMO
Working List of Visual Meteor Showers.

The catalog is intentionally pared down to just the ~38 showers on this
list (rather than the IAU Meteor Data Center's full ~113-shower orbital
catalog): those are the showers with a real, currently-maintained ZHR and
population index behind them. The IAU MDC (still used here only for
parent-body cross-referencing) is an orbital/radiant database and
publishes neither ZHR nor population index.

IMO_WORKING_LIST below is transcribed from Table 5, "Working List of
Visual Meteor Showers", in the IMO's own Meteor Shower Calendar --
https://www.imo.net/ShCal27s.pdf (IMO INFO(3-26), the 2027 edition,
"correct according to the best information available in June 2026, with
maximum dates accurate only for 2027"). Re-transcribe this table from a
current edition periodically; IMO revises it as new observations come in.
Per IMO's usual citation practice, work using this data should credit:
  Molau, S.; Barentsen, G.; Crivello, S.; et al., International Meteor
    Organization, Meteor Shower Calendar.

ZHR values ending in the source's own "+" (e.g. Quadrantids' "80+") are
recorded here as that base number -- IMO's own notation for "this much or
more", not a precise figure. Showers marked "Var" (var iable/outburst-only,
e.g. alpha-Monocerotids) have no baseline ZHR at all and are recorded with
zhr=None.

Real published activity decay rates (b_rise/b_decline -- IMO's table
doesn't include these) come from the Global Meteor Network's operational
flux-monitoring data -- GMN_DECAY_RATES below is transcribed from
share/flux_showers.csv in github.com/CroatianMeteorNetwork/RMS (the
camera-network software behind https://globalmeteornetwork.org/flux/),
its "NASA meteoroid Environment Office" section only -- the other
sections (showers outside that list, and single-year outbursts) use an
undifferentiated Bp=Bm=0.2 placeholder rather than a real per-shower fit,
so they're excluded here rather than presented as measured data. GMN uses
the identical asymmetric log-linear decay model this package does
(RMS/Formats/Showers.py Shower.computeZHRFloat); a `Bm` of exactly 0 in
their table means "symmetric, reuse Bp" (their own fallback, replicated
here at transcription time), not "no decline". Per GMN's usual citation
practice, work using this data should credit:
  Vida, D.; Segon, D.; Gural, P.S.; et al. 2021, "The Global Meteor
    Network -- Methodology and First Results", MNRAS 506(4), 5046-5074.
  Vida, D.; Blaauw Erskine, R.C.; Brown, P.G.; et al. 2022, "Computing
    optical meteor flux using Global Meteor Network data", MNRAS 515(2),
    2322-2339.

Only 21 of these 38 showers have a credible (non-placeholder) GMN decay
rate -- see GMN_DECAY_RATES. The rest get real ZHR/population-index/
radiant/activity data but no rate estimate (nothing to build a decay
curve from).

Per the IAU MDC's own citation request, work using its parent-body data
should also cite:
  Jenniskens, P.; Jopek, T.J.; Janches, D.; Hajdukova, M.; Kokhirova, G.I.;
    Rudawska, R.; 2020, Planetary and Space Science, Vol. 182, 104821.
  Jopek, T.J.; Kanuchova, Z.; 2017, Planetary and Space Science, 143, 3-6.

Usage:
    uv run python scripts/build_shower_catalog.py
"""

from __future__ import annotations

import csv
import importlib.util
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

# Load astro.py directly by path rather than `from meteor_shower_forecast
# import astro` -- the package's __init__ imports catalog.py, which reads
# the very CSV this script is about to (re)generate, so a normal package
# import would fail on a clean checkout / after deleting the CSV.
_spec = importlib.util.spec_from_file_location(
    "astro", REPO_ROOT / "src" / "meteor_shower_forecast" / "astro.py"
)
astro = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(astro)

MDC_URL = "https://www.ta3.sk/IAUC22DB/MDC2022/Etc/streamestablisheddata2026.txt"
CACHE_PATH = REPO_ROOT / ".cache" / "established_showers.txt"
OUT_PATH = REPO_ROOT / "src" / "meteor_shower_forecast" / "data" / "shower_catalog.csv"

MDC_COL = dict(iau_number=1, ad_no=2, code=3, origin=31, ra=11, de=12, los=10, n_orbits=28)

REFERENCE_YEAR = 2027  # matches the IMO edition transcribed below; only month/day/hour are kept

# name, iau_number, code, solar longitude of max (deg), radiant RA/Dec
# (J2000 deg), ZHR ("+"-suffixed values kept as their base number; None ==
# "Var", no baseline rate). Transcribed from IMO Table 5 -- see module
# docstring.
IMO_WORKING_LIST = [
    ("Quadrantids", 10, "QUA", 283.15, 230.0, 49.0, 80.0),
    ("gamma-Ursae Minorids", 404, "GUM", 298.0, 228.0, 67.0, 3.0),
    ("alpha-Centaurids", 102, "ACE", 319.4, 211.0, -58.0, 6.0),
    ("April Lyrids", 6, "LYR", 32.32, 271.0, 34.0, 18.0),
    ("pi-Puppids", 137, "PPU", 33.5, 110.0, -45.0, None),
    ("eta-Aquariids", 31, "ETA", 45.5, 338.0, -1.0, 50.0),
    ("eta-Lyrids", 145, "ELY", 50.0, 291.0, 43.0, 3.0),
    ("Daytime Arietids", 171, "ARI", 76.7, 43.0, 24.0, 30.0),
    ("June Bootids", 170, "JBO", 90.3, 221.0, 48.0, None),
    ("July Pegasids", 175, "JPE", 108.0, 347.0, 11.0, 3.0),
    ("July gamma-Draconids", 184, "GDR", 125.13, 280.0, 51.0, 5.0),
    ("alpha-Capricornids", 1, "CAP", 128.0, 307.0, -10.0, 5.0),
    ("Southern delta-Aquariids", 5, "SDA", 128.0, 340.0, -16.0, 25.0),
    ("eta-Eridanids", 191, "ERI", 135.0, 41.0, -11.0, 3.0),
    ("Perseids", 7, "PER", 140.0, 48.0, 58.0, 110.0),
    ("kappa-Cygnids", 12, "KCG", 144.0, 286.0, 59.0, 3.0),
    ("Aurigids", 206, "AUR", 158.6, 91.0, 39.0, 6.0),
    ("September epsilon-Perseids", 208, "SPE", 166.7, 48.0, 40.0, 8.0),
    ("September Lyncids", 81, "SLY", 170.0, 113.0, 56.0, 3.0),
    ("Daytime Sextantids", 221, "DSX", 188.0, 156.0, -2.0, 5.0),
    ("October Camelopardalids", 281, "OCT", 192.58, 164.0, 79.0, 5.0),
    ("Draconids", 9, "DRA", 195.4, 263.0, 56.0, 5.0),
    ("epsilon-Geminids", 23, "EGE", 205.0, 102.0, 27.0, 3.0),
    ("Orionids", 8, "ORI", 208.0, 95.0, 16.0, 20.0),
    ("Leonis Minorids", 22, "LMI", 211.0, 162.0, 37.0, 2.0),
    ("Southern Taurids", 2, "STA", 223.0, 52.0, 15.0, 7.0),
    ("Northern Taurids", 17, "NTA", 230.0, 58.0, 22.0, 5.0),
    ("Leonids", 13, "LEO", 235.27, 152.0, 22.0, 15.0),
    ("alpha-Monocerotids", 246, "AMO", 239.32, 117.0, 1.0, None),
    ("November Orionids", 250, "NOO", 246.0, 91.0, 16.0, 3.0),
    ("Phoenicids", 254, "PHO", 249.5, 8.0, -27.0, None),
    ("Andromedids", 18, "AND", 254.0, 25.0, 51.0, 5.0),
    ("Puppid-Velids", 301, "PUP", 255.0, 123.0, -45.0, 10.0),
    ("Monocerotids", 19, "MON", 257.0, 100.0, 8.0, 3.0),
    ("sigma-Hydrids", 16, "HYD", 257.0, 125.0, 2.0, 7.0),
    ("Geminids", 4, "GEM", 262.2, 112.0, 33.0, 150.0),
    ("Ursids", 15, "URS", 270.7, 217.0, 76.0, 10.0),
    ("Comae Berenicids", 20, "COM", 271.0, 164.0, 29.0, 3.0),
]

# Population index r per code, transcribed from the same IMO table row.
IMO_POPULATION_INDEX = {
    "QUA": 2.1, "GUM": 3.0, "ACE": 2.0, "LYR": 2.1, "PPU": 2.0, "ETA": 2.4,
    "ELY": 3.0, "ARI": 2.8, "JBO": 2.2, "JPE": 3.0, "GDR": 3.0, "CAP": 2.5,
    "SDA": 2.5, "ERI": 3.0, "PER": 2.2, "KCG": 3.0, "AUR": 2.5, "SPE": 2.5,
    "SLY": 3.0, "DSX": 2.5, "OCT": 2.5, "DRA": 2.6, "EGE": 3.0, "ORI": 2.5,
    "LMI": 3.0, "STA": 2.3, "NTA": 2.3, "LEO": 2.5, "AMO": 2.4, "NOO": 3.0,
    "PHO": 2.8, "AND": 3.0, "PUP": 2.9, "MON": 3.0, "HYD": 3.0, "GEM": 2.6,
    "URS": 2.8, "COM": 3.0,
}

# (Bp, Bm) per code, transcribed verbatim from GMN's flux_showers.csv
# "NASA meteoroid Environment Office" section (see module docstring) --
# these become b_rise/b_decline directly (GMN's Bp/Bm are per degree of
# solar longitude; this package treats that as per-day, ~1.5% error most
# of the year, per effective_zhr's own documented caveat). Two showers in
# GMN's list (Southern mu Sagittariids/SSG, Piscis Austrinids/PAU) aren't
# on the IMO Working List this catalog is otherwise built from, so they're
# omitted -- GMN_DECAY_RATES is looked up by code, so missing keys are
# harmless, but there's no point carrying entries that can never match.
#
# A GMN Bm of exactly 0 means "symmetric activity, reuse Bp" (their own
# fallback in Shower.computeZHRFloat) rather than "no decline" -- applied
# here at transcription time, not literally recorded as 0.
GMN_DECAY_RATES = {
    "CAP": (0.059, 0.091),
    "STA": (0.026, 0.026),
    "GEM": (0.150, 0.462),
    "SDA": (0.109, 0.071),
    "LYR": (1.229, 4.620),
    "PER": (0.350, 0.350),
    "ORI": (0.120, 0.119),
    "DRA": (7.225, 7.225),
    "QUA": (0.865, 0.827),
    "KCG": (0.069, 0.069),
    "LEO": (0.550, 0.550),
    "URS": (2.292, 1.258),
    "HYD": (0.100, 0.100),
    "NTA": (0.026, 0.026),
    "MON": (0.250, 0.250),
    "LMI": (0.093, 0.079),
    "ETA": (0.125, 0.079),
    "ARI": (0.065, 0.055),
    "AUR": (0.190, 0.190),
    "SPE": (0.193, 0.193),
    "DSX": (0.063, 0.167),
}


def _download(url: str, dest: Path) -> None:
    if dest.exists():
        print(f"using cached {dest}")
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"downloading {url}")
    with urllib.request.urlopen(url, timeout=60) as resp:
        dest.write_bytes(resp.read())


def _parse_mdc_rows(path: Path) -> list[list[str]]:
    rows = []
    with path.open(encoding="latin-1") as f:
        for line in f:
            if not line.startswith('"'):
                continue
            rows.append([c.strip().strip('"').strip() for c in line.rstrip("\n").split("|")])
    return rows


def _int_or(s: str, default: int = 0) -> int:
    try:
        return int(s)
    except ValueError:
        return default


def _parent_bodies_by_code(mdc_rows: list[list[str]]) -> dict[str, str]:
    """Best (largest-orbit-sample) non-empty Origin field per IAU code."""
    by_code: dict[str, list[list[str]]] = {}
    for r in mdc_rows:
        by_code.setdefault(r[MDC_COL["code"]], []).append(r)
    result = {}
    for code, recs in by_code.items():
        with_origin = [r for r in recs if r[MDC_COL["origin"]]]
        if not with_origin:
            continue
        rep = max(with_origin, key=lambda r: (_int_or(r[MDC_COL["n_orbits"]]), _int_or(r[MDC_COL["ad_no"]])))
        result[code] = rep[MDC_COL["origin"]]
    return result


def build() -> None:
    _download(MDC_URL, CACHE_PATH)
    mdc_rows = _parse_mdc_rows(CACHE_PATH)
    parent_bodies = _parent_bodies_by_code(mdc_rows)

    out_rows = []
    for name, iau_number, code, los_deg, ra_deg, dec_deg, zhr in IMO_WORKING_LIST:
        peak_dt = astro.date_for_solar_longitude(los_deg, REFERENCE_YEAR)
        b_rise, b_decline = GMN_DECAY_RATES.get(code, (None, None))
        out_rows.append(dict(
            name=name, iau_code=code, iau_number=iau_number,
            radiant_ra_deg=round(ra_deg, 3), radiant_dec_deg=round(dec_deg, 3),
            peak_month=peak_dt.month, peak_day=peak_dt.day,
            peak_hour_utc=round(peak_dt.hour + peak_dt.minute / 60 + peak_dt.second / 3600, 2),
            parent_body=parent_bodies.get(code, ""),
            peak_zhr=zhr if zhr is not None else "",
            population_index=IMO_POPULATION_INDEX[code],
            b_rise=b_rise if b_rise is not None else "",
            b_decline=b_decline if b_decline is not None else "",
        ))

    out_rows.sort(key=lambda r: (r["peak_month"], r["peak_day"]))

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "name", "iau_code", "iau_number", "radiant_ra_deg", "radiant_dec_deg",
        "peak_month", "peak_day", "peak_hour_utc", "parent_body",
        "peak_zhr", "population_index", "b_rise", "b_decline",
    ]
    with OUT_PATH.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(out_rows)

    with_rate = sum(1 for r in out_rows if r["b_rise"] != "")
    with_parent = sum(1 for r in out_rows if r["parent_body"])
    print(
        f"wrote {len(out_rows)} showers to {OUT_PATH} "
        f"({with_rate} with a full rate estimate, {with_parent} with a known parent body)"
    )


if __name__ == "__main__":
    build()
