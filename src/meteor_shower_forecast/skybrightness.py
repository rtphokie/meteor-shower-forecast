"""Meteor-specific sky-brightness helper: scattered-moonlight sky glow.

Bortle<->SQM conversion and naked-eye-limiting-magnitude are no longer
defined here -- they live in the `bortlefinder` package (PyPI), extracted
from this project's own earlier coordinate-based Bortle estimation. This
module keeps only what's specific to meteor rate estimation and isn't part
of bortlefinder's scope: how much a given moon illumination brightens the
sky, adapted from the DarkHours project (mbeher2200/DarkHours,
darkhours/moonlight.py), MIT licensed:

    Copyright (c) 2026 DarkHours contributors
"""

from __future__ import annotations

import math

# --- Simplified lunar sky-glow model --------------------------------------
#
# Full K&S(1991)/Winkler(2022) moon-target-observer geometry needs live moon
# ephemeris (RA/Dec/altitude/distance), which is out of scope for this
# project. When the caller supplies a moon illumination percentage, this
# reuses DarkHours' own site-wide reference geometry simplification
# (90 deg moon-target separation, 30 deg moon altitude -- see DarkHours
# moonlight.py `_KS_CREDIT_SEP_DEG`/`_KS_CREDIT_ALT_DEG`, "darkest
# accessible sky, representative mid-sky moon position") rather than a
# fabricated per-target geometry.

_KS_K_EXT = 0.172  # V-band extinction coefficient (mag/airmass), clear reference sky
_KS_NORM = 24130491.121213324  # kernel normalisation, ported from DarkHours moonlight.py
_HG_G = 0.8  # Henyey-Greenstein asymmetry parameter (aerosol forward-scatter)
_TAU_RAY = 0.1066
_TAU_ABS = 0.016
_MAG_PER_TAU = 2.5 * math.log10(math.e)
_AOD_REF = _KS_K_EXT / _MAG_PER_TAU - _TAU_RAY - _TAU_ABS

REFERENCE_MOON_SEP_DEG = 90.0
REFERENCE_MOON_ALT_DEG = 30.0


def _pathlength(alt_deg: float) -> float:
    z = math.radians(90.0 - max(0.0, alt_deg))
    return (1.0 - 0.96 * math.sin(z) ** 2) ** -0.5


def moon_delta_mag(
    illumination_pct: float,
    sky_sqm: float,
    target_alt_deg: float,
    sep_deg: float = REFERENCE_MOON_SEP_DEG,
    moon_alt_deg: float = REFERENCE_MOON_ALT_DEG,
) -> float:
    """Sky surface-brightness increase (mag/arcsec^2) from scattered
    moonlight, at reference moon geometry, for a given illumination
    fraction. Returns 0.0 when illumination is zero.
    """
    if illumination_pct <= 0 or moon_alt_deg <= 0:
        return 0.0

    illum = illumination_pct / 100.0
    alpha = math.degrees(math.acos(max(-1.0, min(1.0, 2.0 * illum - 1.0))))
    v_moon = -12.73 + 0.026 * alpha + 4e-9 * alpha**4
    i_moon = 10 ** (-0.4 * (v_moon + 16.57))

    tau_s = _TAU_RAY + _AOD_REF
    tau = tau_s + _TAU_ABS

    rho = math.radians(max(0.1, sep_deg))
    cos_rho = math.cos(rho)
    p_ray = 3.0 / (16.0 * math.pi) * (1.0 + cos_rho**2)
    p_mie = (1.0 - _HG_G**2) / (4.0 * math.pi * (1.0 + _HG_G**2 - 2.0 * _HG_G * cos_rho) ** 1.5)
    p = (_TAU_RAY * p_ray + _AOD_REF * p_mie) / tau_s

    x_t = _pathlength(target_alt_deg)
    x_m = _pathlength(moon_alt_deg)
    if abs(x_m - x_t) < 1e-6:
        kernel = tau * x_t * math.exp(-tau * x_t)
    else:
        kernel = x_t * (math.exp(-tau * x_t) - math.exp(-tau * x_m)) / (x_m - x_t)

    i_scatter = _KS_NORM * p * (tau_s / tau) * kernel * i_moon
    i_sky = 10 ** ((27.78 - sky_sqm) / 2.5)
    return 2.5 * math.log10(1.0 + i_scatter / i_sky)
