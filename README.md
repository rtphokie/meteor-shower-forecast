# meteor-shower-forecast

Location meteor-shower rate estimates, not the inflated, best case,
[Zenithal Hourly Rate](https://www.amsmeteors.org/glossary/#zenithal-hourly-rate)
(ZHR), but what a naked-eye observer at a specific place and time should
actually expect to see.

A shower's published ZHR assumes the
[radiant](https://www.amsmeteors.org/glossary/#radiant) is at the zenith
and the sky is dark enough to see magnitude 6.5 stars. Neither is usually
true, so this combines three corrections:

1. **Day-decayed ZHR** -- an asymmetric log-linear falloff around the
   shower's peak date (showers typically ramp up faster than they decline).
2. **Radiant altitude** -- computed from the observer's coordinates and
   time; a radiant near the horizon is heavily foreshortened
   (`rate ∝ sin(radiant altitude)`). A radiant up to ~10 degrees below the
   local horizon can still produce visible meteors (they ablate ~100km up,
   well above the radiant's own point on the horizon); that ~10 degrees is
   a derived geometric estimate, not a value from meteor-science
   literature -- see `forecast.HORIZON_DIP_DEG`.
3. **Sky brightness** -- the site's
   [Bortle class](https://en.wikipedia.org/wiki/Bortle_scale), combined
   with the actual moon illumination on the given date, degrades the
   naked-eye [limiting magnitude](https://en.wikipedia.org/wiki/Limiting_magnitude),
   which disproportionately hides a shower's fainter meteors according to
   its [population index](https://en.wikipedia.org/wiki/Population_index).

Both the Bortle class and the moon illumination are estimated
automatically (from coordinates and date/time, respectively) unless you
override them.

The rate-estimation formulas (the day-decay/radiant-altitude/sky-
brightness combination) are adapted from the
[DarkHours](https://github.com/mbeher2200/DarkHours) project (MIT
licensed). DarkHours has no PyPI package, so this reimplements the
relevant math directly (see the module docstrings in
`src/meteor_shower_forecast/` for exactly what was ported from where)
rather than vendoring its much larger, service-oriented codebase.

The ~38-shower catalog (names, IAU codes, radiants, activity, ZHR,
population index) is transcribed from the
[International Meteor Organization](https://www.imo.net/)'s Meteor
Shower Calendar, Table 5, "Working List of Visual Meteor Showers" --
the set of showers with a real, currently-maintained ZHR behind them.
Real activity decay rates -- required, alongside ZHR/population index,
for an actual rate estimate -- come from the
[Global Meteor Network](https://globalmeteornetwork.org/flux/)'s
operational flux-monitoring data (21 of the 38 showers have a credible
one). Parent bodies are cross-referenced from the
[IAU Meteor Data Center](https://www.iaumeteordatacenter.org/)'s List of
Established Showers; per its citation request, work using it should cite
Jenniskens et al. 2020 (*Planetary and Space Science*, 182, 104821) and
Jopek & Kanuchova 2017 (*Planetary and Space Science*, 143, 3-6) -- see
`catalog.py`'s docstring for both citations and exact source editions.

Runtime dependencies are kept minimal: [`bortlefinder`](https://pypi.org/project/bortlefinder/)
(coordinates -> Bortle class/SQM, extracted from this project's own
earlier coordinate-based Bortle estimation and published standalone) and
`timezonefinder` (coordinates -> IANA timezone, for local-time display).

## Install

```
pip install meteor-shower-forecast
```

## Usage

With just coordinates, the CLI reports the next upcoming shower, evaluated
at *this observer's* best upcoming viewing opportunity -- not the
shower's catalog peak. The catalog peak is a global, radiant-independent
instant (when Earth crosses the dust trail most densely); it says nothing
about whether the radiant is even above your horizon then. A bracketed
search across the shower's active window (restricted to astronomical
night, so it never picks a daytime moment) finds the actual best time.
`peak` is still reported alongside `best` for reference. Sky brightness
and moon illumination are estimated automatically, and dates/times are
shown in local time with the 3-letter timezone abbreviation (EDT, PST,
JST, ...) looked up from the coordinates:

```
meteor-shower-forecast --lat 35.7796 --lon -78.6382
```

```
Orionids (ORI)
  parent body:         1P/Halley
  best:                Wed, Oct 21, 2026, 5:47 AM EDT
  peak:                Wed, Oct 21, 2026, 11:13 AM EDT
  days from peak:      -0.2
  effective ZHR:       18.8
  radiant altitude:    70.0° (az 188.6°)
  sky brightness:      Bortle 8 (City sky, SQM 17.75) (estimated from coordinates)
  moon:                78% illuminated (estimated for this date/time)
  limiting magnitude:  3.67
  faint-meteor factor: 0.075
  ~ estimated rate:    1.3 meteors/hour
  note: moonlight brightens the sky by 0.10 mag/arcsec²
  note: +1.0 more meteors/hour than observing at the shower's catalog peak from this location
```

This matters most when the catalog peak itself is a bad time to look --
e.g. the Quadrantids' 2027 peak falls with the radiant just below the
horizon from this location (an estimated 3.4 meteors/hour), but a few
hours later, once the radiant has climbed well clear of the horizon, the
same shower produces over five times the rate:

```
meteor-shower-forecast --lat 35.0 --lon -78.6 --shower Quadrantids
```

```
Quadrantids (QUA)
  parent body:         2003 EH1
  best:                Mon, Jan 4, 2027, 5:51 AM EST
  peak:                Sun, Jan 3, 2027, 7:21 PM EST
  days from peak:      +0.4
  effective ZHR:       34.8
  radiant altitude:    56.4° (az 52.3°)
  ...
  ~ estimated rate:    18.8 meteors/hour
  note: +15.4 more meteors/hour than observing at the shower's catalog peak from this location
```

Pass `--year YYYY` to look at a different year's occurrence instead of the
next upcoming one (past years work too) -- it still searches for that
year's best time, unlike `--when`, which pins one exact instant (the two
are mutually exclusive):

```
meteor-shower-forecast --lat 35.0 --lon -78.6 --shower Quadrantids --year 2020
```

`--shower NAME` without `--when` works the same way: it searches that
shower's own upcoming active window for the best time. Pass `--when`
explicitly to evaluate one specific date and time instead (e.g. checking a
shower a day past peak) -- the field is then labeled `time`, not `best`,
since it's exactly what you asked for rather than a search result.

Override the automatic moon estimate, or pass 0 to ignore moonlight
entirely:

```
meteor-shower-forecast --lat 35.0 --lon -78.6 \
  --shower Perseids --when 2026-08-12T05:00:00 --moon-illumination 80
```

Compare every shower's realistic best-case rate, each at its own best
upcoming observing time, sorted chronologically by (catalog) peak date
(pass `--when` to instead compare what's active on one specific night):

```
meteor-shower-forecast --lat 35.0 --lon -78.6 --shower all
```

Pass `--bortle <1-9>` or `--sqm <value>` to override the coordinate-based
sky-brightness estimate with a known Bortle class or a direct
[sky quality meter](https://en.wikipedia.org/wiki/Sky_quality_meter)
reading (mag/arcsec²).

List the known showers (all ~38 on the IMO Working List, each with its
IAU code and parent body when known; look one up by name, prefix, or IAU
code -- e.g. `--shower SDA`):

```
meteor-shower-forecast --list-showers
```

```
QUA   Quadrantids                      peak 01-04  ZHR 80, parent 2003 EH1
GUM   gamma-Ursae Minorids             peak 01-18  ZHR 3, but no rate estimate (no decay-rate data)
...
PPU   pi-Puppids                       peak 04-23  no rate data (variable/outburst shower), parent 26P/Grigg-Skjellerup
```

### Showers without a rate estimate

A rate estimate needs ZHR, population index, *and* an activity decay
rate. IMO's calendar publishes the first two for every shower on its
Working List, but not decay rates -- those come from the Global Meteor
Network's flux-monitoring data (see catalog.py) for 21 showers with a
credible per-shower fit (`ShowerDef.has_rate_data`); those 21 are what
`--shower all` and the plain default compare. A few showers
(variable/outburst-only, e.g. alpha-Monocerotids) have no baseline ZHR at
all. Looking up any other shower by name/code still reports its radiant,
ZHR/population index (if known), parent body, and next peak date;
`estimated_rate_per_hour` (and the fields that feed it) are `None`.

### Estimating Bortle class from coordinates

By default, sky brightness is derived from `--lat`/`--lon` via the
[`bortlefinder`](https://pypi.org/project/bortlefinder/) package, which
reads a small, locally-fetched grid based on the Falchi et al. (2016)
*New World Atlas of Artificial Night Sky Brightness*
([doi.org/10.5880/GFZ.1.4.2016.001](https://doi.org/10.5880/GFZ.1.4.2016.001),
CC BY-NC 4.0). Fetch it once (downloads ~684MB, reduced to a ~19MB local
grid; needs the `build` extra):

```
pip install meteor-shower-forecast[build]
bortlefinder-fetch-grid
```

This derivation is approximate: satellite atlases measure *zenith*
brightness, while the Bortle scale is a subjective whole-sky rating
dominated by horizon light domes. David Lorenz's own validation against
paired dark-sky observations found real disagreement -- often a full
class, worse in the Bortle 5-7 range
([source](https://djlorenz.github.io/astronomy/lp/bortle.html)). Treat it
as a starting estimate, not a substitute for a real SQM meter reading --
every result marks whether its Bortle class was estimated or given
directly. See `bortlefinder`'s own README for more on its accuracy and
licensing (the fetched grid data is CC BY-NC 4.0, non-commercial; the
package code itself is MIT).

Moon illumination is likewise a low-precision approximation (mean lunar
elongation, no ephemeris) good to within a few percent -- see
`astro.moon_illumination_pct`.

The best-time search only considers
[astronomical night](https://en.wikipedia.org/wiki/Twilight#Astronomical_twilight)
(Sun below 18 degrees altitude), matching the fully-dark-sky assumption
already baked into the Bortle/SQM model, and never proposes a time before
now. If a
shower's whole active window at a given location never reaches
astronomical darkness (e.g. high-latitude summer), it falls back to the
least-bright moment found and says so in a note.

## Library

```python
from meteor_shower_forecast import estimate_shower_rate, estimate_all_showers, next_shower

estimate = estimate_shower_rate("Perseids", lat=35.0, lon=-78.6, bortle=4)
print(estimate.iau_code, estimate.parent_body, estimate.estimated_rate_per_hour)

for e in estimate_all_showers(lat=35.0, lon=-78.6, bortle=4):
    print(e.shower, e.estimated_rate_per_hour)

upcoming = next_shower()
print(upcoming.name)

# Any of the ~38 IMO-listed showers can be looked up by name/code;
# estimated_rate_per_hour is None where decay-rate data doesn't exist
# (alpha-Capricornids has a real ZHR but no published decay rate).
capricornids = estimate_shower_rate("CAP", lat=35.0, lon=-78.6, bortle=4)
print(capricornids.radiant_alt_deg, capricornids.estimated_rate_per_hour)

# reference anchors the best-time search to a given year (or any datetime)
# instead of now -- past years work too. Ignored if when_utc is passed.
from datetime import datetime, timezone
past = estimate_shower_rate(
    "Geminids", lat=35.0, lon=-78.6, bortle=4, reference=datetime(2020, 1, 1, tzinfo=timezone.utc)
)
print(past.when_utc, past.estimated_rate_per_hour)
```

## Development

```
git clone https://github.com/rtphokie/meteor-shower-forecast
cd meteor-shower-forecast
uv sync
```

The shower catalog (`src/meteor_shower_forecast/data/shower_catalog.csv`)
is generated, not hand-edited -- see `scripts/build_shower_catalog.py`'s
docstring to regenerate it from a newer IMO calendar or updated GMN data.
