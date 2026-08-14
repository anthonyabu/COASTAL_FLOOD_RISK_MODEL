"""
Coastal hazard forcing, depth-damage functions, and expected annual damage.

PROVENANCE NOTICE
-----------------
  1. SURGE DISTRIBUTION - FITTED. A GEV fitted by maximum likelihood to 46 annual
     maxima of the tide-removed surge residual at NOAA CO-OPS station 8418150
     (Portland, Maine), 1980-2025, by code/fit_surge_gev_v2.R. Cross-checked against
     the extRemes package (agreement to four decimal places). The pipeline was
     validated independently: the sea level trend recovered from the removed datum
     offset is +2.83 mm/yr (p < 0.0001), comparable with the published rate for the
     gauge, and was not fitted to anything.

     NOT fitted: the tide amplitude (1.40 m) is assumed and should be checked against
     the published NOAA datums for the station. Sea level rise was removed before
     fitting, so this distribution is meteorological surge only and SLR must be added
     separately to still water level.

  2. ROAD DEPTH-DAMAGE CURVE - STILL ILLUSTRATIVE, and now the largest single source
     of uncertainty in the result. Published road depth-damage functions exist (HAZUS
     transportation modules and several European road-specific studies) and one must
     be substituted. The curve here has a plausible SHAPE - negligible damage until
     water is deep enough to reach the pavement structure, steep rise through
     0.3-1.0 m, saturating near 2-3 m - but the ordinates are not taken from a source
     and are not defensible on their own. Swapping between the two candidate curves
     in DD_VARIANTS changes expected annual damage by a factor of 7.6.

Provenance is printed at runtime so a reader cannot mistake one for the other.
"""

import numpy as np

# ---------------------------------------------------------------------------
# Tide + surge forcing
# ---------------------------------------------------------------------------

TIDE_AMPLITUDE = 1.40      # m, approx MHW above MSL at Portland (2.8 m range)
                           # VERIFY against the published NOAA datums for the station
TIDE_PERIOD = 12.42 * 3600.0
SURGE_DURATION = 12.0 * 3600.0
EVENT_DURATION = 18.0 * 3600.0   # 12 h surge + 6 h drainage

# Extreme value distribution for the annual maximum surge residual.
#
# Coles parameterisation: xi > 0 heavy (Frechet) tail, xi = 0 Gumbel,
# xi < 0 bounded (Weibull) tail.
# WARNING: scipy.stats.genextreme uses c = -xi. Mixing the two conventions
# silently flips the tail and corrupts every return level.
#
# These values are ILLUSTRATIVE until replaced by the output of
# code/fit_surge_gev.R, which fits them to NOAA CO-OPS tide gauge records.
GEV_LOC   = 0.6900         # NOAA station 8418150 (Portland, Maine), 1980-2025
GEV_SCALE = 0.1353         # n = 46 annual maxima of the tide-removed surge residual
GEV_SHAPE = -0.1319        # Coles convention; scipy genextreme uses c = -xi
GEV_IS_FITTED = True       # fitted by code/fit_surge_gev_v2.R, cross-checked in extRemes
#
# The negative shape implies a BOUNDED upper tail. Assuming Gumbel with the same
# location and scale would overstate the 100-year surge by 13.5% and the 500-year
# by 21.1%.
#
# Mean sea level rise was removed before fitting (the annual mean residual is
# subtracted year by year), so this is meteorological surge only. The gauge shows
# +2.83 mm/yr of sea level rise over this window; that must be added SEPARATELY to
# still water level, and is not included below.


def surge_for_return_period(T_years, loc=None, scale=None, shape=None):
    """
    Peak surge residual (m) for a given return period, GEV inverse CDF.

    The shape parameter matters more than it looks. Holding location and scale
    fixed, moving xi from 0.0 to +0.2 raises the 100-year surge by 40 %, and
    moving it to -0.2 lowers it by 22 %. Assuming Gumbel when the data support a
    heavy tail systematically understates rare events.
    """
    mu = GEV_LOC if loc is None else loc
    sig = GEV_SCALE if scale is None else scale
    xi = GEV_SHAPE if shape is None else shape

    y = -np.log(1.0 - 1.0 / np.asarray(T_years, dtype=float))
    if abs(xi) < 1e-8:
        return mu - sig * np.log(y)
    return mu + (sig / xi) * (y ** (-xi) - 1.0)


def water_level(t, surge_peak, surge_duration=None):
    """
    Still water level (m MSL) at time t (s).

    Astronomical high water is aligned with the surge peak at t = 6 h. This is the
    standard conservative assumption for a design event: it is the joint-occurrence
    case, not the expected case, and a full analysis would treat tide-surge timing
    as a joint probability rather than assuming coincidence.
    """
    Td = SURGE_DURATION if surge_duration is None else surge_duration
    t_peak = 0.5 * Td
    tide = TIDE_AMPLITUDE * np.cos(2 * np.pi * (t - t_peak) / TIDE_PERIOD)
    if t <= Td:
        surge = surge_peak * 0.5 * (1.0 - np.cos(2 * np.pi * t / Td))
    else:
        surge = 0.0
    return tide + surge


def peak_still_water_level(surge_peak):
    return TIDE_AMPLITUDE + surge_peak


# ---------------------------------------------------------------------------
# Depth-damage
# ---------------------------------------------------------------------------

# ILLUSTRATIVE road depth-damage curves (flexible pavement, 2-lane).
#
# Two variants are kept deliberately, because the choice between them turns out to
# dominate the expected-annual-damage result:
#
#   "no_threshold"  damage begins as soon as water reaches the carriageway.
#   "threshold"     damage is negligible until water is deep enough to affect the
#                   pavement structure and subgrade rather than merely covering the
#                   surface. Shallow nuisance flooding closes a road without
#                   structurally damaging it.
#
# The threshold variant is the default. The no-threshold variant is retained so the
# sensitivity of EAD to this single modelling choice can be quantified rather than
# buried.
DD_VARIANTS = {
    "no_threshold": (np.array([0.00, 0.10, 0.25, 0.50, 1.00, 2.00, 3.00]),
                     np.array([0.00, 0.05, 0.15, 0.32, 0.55, 0.85, 1.00])),
    "threshold":    (np.array([0.00, 0.20, 0.30, 0.50, 1.00, 2.00, 3.00]),
                     np.array([0.00, 0.00, 0.04, 0.15, 0.40, 0.75, 1.00])),
}

DEFAULT_VARIANT = "threshold"

# Backwards-compatible handles for plotting the default curve
DD_DEPTH, DD_RATIO = DD_VARIANTS[DEFAULT_VARIANT]


def damage_ratio(depth_m, variant=DEFAULT_VARIANT):
    """Fraction of replacement value lost, from peak inundation depth."""
    d, r = DD_VARIANTS[variant]
    return np.interp(np.asarray(depth_m, dtype=float), d, r, left=0.0, right=1.0)


def duration_multiplier(hours, enabled=False):
    """
    OPTIONAL and OFF BY DEFAULT.

    Prolonged submergence saturates the subgrade and increases damage beyond what
    peak depth alone implies. The physical mechanism is real, but the multiplier
    values below are invented, so this is disabled unless explicitly switched on for
    a sensitivity test. It must not contribute to any headline number.
    """
    if not enabled:
        return np.ones_like(np.asarray(hours, dtype=float))
    return np.interp(hours, [0.0, 6.0, 24.0, 72.0], [1.00, 1.00, 1.25, 1.40])


def assess_assets(assets, h_max, wet_duration_s, z=None, deck=None,
                  variant=DEFAULT_VARIANT):
    """
    Compute per-segment peak depth over the ASSET SURFACE, damage ratio and value.

    Depth over the asset is (peak water surface elevation - deck elevation), not the
    water depth over the bed. The two differ wherever the carriageway is not the
    ground surface - most importantly at bridges, where the channel bed beneath the
    opening can carry several metres of water while the deck above stays dry.

    Depth on a segment is taken as the MAXIMUM over its centreline cells, not the
    mean. For a linear asset like a road, functional loss is governed by the worst
    point on the segment - a single washed-out 20 m length closes the whole segment.
    Using the mean would systematically understate damage.
    """
    if z is None or deck is None:
        raise ValueError(
            "assess_assets requires both z (bed elevation) and deck (asset surface "
            "elevation). Passing only depth silently treats bridge decks as if they "
            "sat on the channel bed.")

    wse_max = h_max + z
    results = []
    for a in assets:
        for seg in a["segments"]:
            depths = np.array([
                max(0.0, wse_max[j, i] - deck[j, i]) for (j, i) in seg["cells"]])
            durs = np.array([wet_duration_s[j, i] for (j, i) in seg["cells"]]) / 3600.0
            d_peak = float(depths.max())
            d_mean = float(depths.mean())
            dur = float(durs.max())
            ratio = float(damage_ratio(d_peak, variant=variant))
            results.append({
                "asset": a["name"],
                "segment": seg["id"],
                "length_m": seg["length_m"],
                "value_usd": seg["value_usd"],
                "depth_peak_m": d_peak,
                "depth_mean_m": d_mean,
                "duration_h": dur,
                "damage_ratio": ratio,
                "damage_usd": ratio * seg["value_usd"],
            })
    return results


def expected_annual_damage(return_periods, damages):
    """
    EAD = integral of damage over annual exceedance probability.

    Trapezoidal integration in probability space. Two truncations are made explicit
    rather than hidden:

      - Above the largest return period sampled, damage is assumed constant at its
        value there. Because damage saturates once assets are fully lost, this is a
        mild assumption, but it is still an assumption.
      - Below the smallest return period sampled, damage is assumed to fall linearly
        to zero at p = 1. If damage is already zero at the smallest T, this term
        vanishes and the truncation is harmless.
    """
    T = np.asarray(return_periods, dtype=float)
    D = np.asarray(damages, dtype=float)
    order = np.argsort(T)
    T, D = T[order], D[order]

    p = 1.0 / T                       # annual exceedance probability, descending in T
    # integrate from p_min (rarest) to p_max (most frequent)
    p_asc = p[::-1]
    D_asc = D[::-1]

    ead = float(np.trapezoid(D_asc, p_asc))

    # tail beyond the rarest event sampled
    tail = float(D[-1] * p[-1])
    # segment between p_max and p = 1, damage tapering linearly to zero
    p_max = p_asc[-1]
    low = float(0.5 * D_asc[-1] * (1.0 - p_max))

    return {
        "ead_usd": ead + tail + low,
        "ead_core_usd": ead,
        "tail_usd": tail,
        "frequent_tail_usd": low,
    }


if __name__ == "__main__":
    src = ("FITTED - NOAA CO-OPS station 8418150 (Portland, Maine), 1980-2025, "
           "n = 46" if GEV_IS_FITTED else "ILLUSTRATIVE - not fitted to any gauge")
    fam = ("Gumbel" if abs(GEV_SHAPE) < 1e-8 else
           f"GEV, shape xi = {GEV_SHAPE:+.4f} "
           f"({'bounded' if GEV_SHAPE < 0 else 'heavy'} tail)")
    print(f"Surge hazard  [{src}]")
    print(f"  distribution: {fam}")
    print(f"  location {GEV_LOC:.4f} m, scale {GEV_SCALE:.4f} m")
    print(f"  tide amplitude {TIDE_AMPLITUDE:.2f} m (ASSUMED - verify against the "
          f"published NOAA datums)")
    print(f"  surge is METEOROLOGICAL ONLY; sea level rise is not included\n")
    print(f"{'T (yr)':>8} {'surge (m)':>11} {'peak SWL (m MSL)':>18} "
          f"{'if Gumbel':>11} {'error':>8}")
    for T in [1.1, 2, 5, 10, 25, 50, 100, 200, 500]:
        v = surge_for_return_period(T)
        g = surge_for_return_period(T, shape=0.0)
        print(f"{T:8.1f} {v:11.2f} {peak_still_water_level(v):18.2f} "
              f"{g:11.2f} {100*(g-v)/v:+7.1f}%")
    print("\n  'if Gumbel' holds location and scale fixed and sets xi = 0, showing\n"
          "  what assuming a Gumbel tail would cost at each return period.")

    print("\nRoad depth-damage curve  [ILLUSTRATIVE - shape is plausible,\n  ordinates are NOT sourced. Substitute a published curve before use.]")
    for d in [0.05, 0.1, 0.3, 0.5, 1.0, 1.5, 2.5]:
        print(f"  depth {d:4.2f} m -> damage ratio {damage_ratio(d):.3f}")
