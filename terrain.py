"""
Synthetic coastal domain with transportation assets.

SYNTHETIC DATA NOTICE
---------------------
The terrain, roads and unit costs in this module are constructed, not surveyed. They
are geometrically plausible for a low-lying north-east US coastal plain but they are
NOT a real site. Replace with a real DEM (e.g. USGS 3DEP / NOAA CoNED lidar) and a real
road network (e.g. state DOT centrelines) before any result is used for anything.

Domain layout, looking inland:

    x = 0 m        open sea boundary (west edge)
    x = 0-2000 m   coastal plain, bed rising from -1.2 m to +3.0 m MSL
    y = 0-1200 m   alongshore extent

Assets:
    A  Coastal Highway   shore-parallel embankment at x = 700 m,  crest 2.35 m
                         with a bridged creek opening where crest drops to 0.3 m
    B  Access Causeway   shore-perpendicular at y = 300 m, x = 200-700 m, crest 2.10 m
    C  Inland Arterial   shore-parallel at x = 1500 m, crest 2.55 m

A tidal creek runs inland at y ~ 600 m. It is the main conveyance path that lets surge
penetrate behind the Coastal Highway embankment. This is deliberate: it reproduces the
common real-world situation where a shore-parallel embankment is outflanked through a
drainage opening rather than overtopped along its length.
"""

import numpy as np

DX = 20.0            # cell size, m
NX = 100             # 2000 m cross-shore
NY = 60              # 1200 m alongshore

# Manning's n by surface type
N_MARSH = 0.10
N_PLAIN = 0.055
N_CREEK = 0.035
N_ROAD = 0.020

SEGMENT_CELLS = 5    # damage segment length = 5 cells = 100 m


def build_domain():
    """
    Return (z, manning, assets, x, y, deck).

    Two elevation fields are maintained deliberately:

      z     bed elevation used by the hydraulics. At a bridge this is the CHANNEL BED,
            because that is what governs conveyance through the opening.

      deck  asset surface elevation used for damage. At a bridge this is the DECK,
            which sits above the opening and is only damaged once water reaches it.

    Conflating the two makes every bridge appear catastrophically flooded in every
    event, because the model reads metres of water over the channel bed as metres of
    water over the carriageway. deck is NaN away from road cells.
    """
    x = np.arange(NX) * DX
    y = np.arange(NY) * DX

    # ---- base coastal plain ----
    # Concave-up profile: a flat marsh apron near the shore steepening inland. This is
    # more representative of a glaciated north-east coastline than a uniform ramp, and
    # it matters here because a uniform ramp floods the whole domain at once and
    # destroys any gradation between return periods.
    profile = -1.2 + 4.7 * (x / 2000.0) ** 1.4
    z = np.tile(profile, (NY, 1))

    # gentle deterministic undulation so the surface is not a perfect plane
    XX, YY = np.meshgrid(x, y)
    z += 0.15 * np.sin(2 * np.pi * XX / 900.0) * np.cos(2 * np.pi * YY / 700.0)

    manning = np.full_like(z, N_PLAIN)
    manning[:, :20] = N_MARSH                 # salt marsh in the first 400 m
    deck = np.full_like(z, np.nan)            # asset surface elevation, NaN off-road

    # ---- tidal creek, j = 29..31, running inland ----
    creek_j = slice(29, 32)
    z[creek_j, :] -= 1.0
    manning[creek_j, :] = N_CREEK

    # ---- assets ----
    assets = []

    # A: Coastal Highway, shore-parallel at x = 700 m (i = 34..36)
    a_i = slice(34, 37)
    z[:, a_i] = np.maximum(z[:, a_i], 2.35)
    deck[:, a_i] = 2.35                        # carriageway level, continuous
    # bridged opening over the creek: the channel bed stays low so flow passes
    # beneath, but the deck above it remains at carriageway level
    z[creek_j, a_i] = np.minimum(z[creek_j, a_i], 0.3)
    manning[:, a_i] = N_ROAD
    manning[creek_j, a_i] = N_CREEK
    assets.append(_make_asset(
        name="A - Coastal Highway",
        cells=[(j, 35) for j in range(NY)],      # centreline
        width_cells=3,
        cost_per_km=2_400_000.0,
        axis="alongshore",
    ))

    # B: Access Causeway, shore-perpendicular at y = 300 m (j = 14..16), x = 200-700 m
    b_j = slice(14, 17)
    b_i = slice(10, 36)
    z[b_j, b_i] = np.maximum(z[b_j, b_i], 2.10)
    deck[b_j, b_i] = z[b_j, b_i]
    manning[b_j, b_i] = N_ROAD
    assets.append(_make_asset(
        name="B - Access Causeway",
        cells=[(15, i) for i in range(10, 36)],
        width_cells=3,
        cost_per_km=3_200_000.0,
        axis="cross-shore",
    ))

    # C: Inland Arterial, shore-parallel at x = 1500 m (i = 74..76)
    c_i = slice(74, 77)
    z[:, c_i] = np.maximum(z[:, c_i], 2.55)
    deck[:, c_i] = z[:, c_i]
    manning[:, c_i] = N_ROAD
    assets.append(_make_asset(
        name="C - Inland Arterial",
        cells=[(j, 75) for j in range(NY)],
        width_cells=3,
        cost_per_km=2_000_000.0,
        axis="alongshore",
    ))

    return z, manning, assets, x, y, deck


def _make_asset(name, cells, width_cells, cost_per_km, axis):
    """Split a road centreline into fixed-length segments for asset-level damage."""
    segments = []
    for s in range(0, len(cells), SEGMENT_CELLS):
        chunk = cells[s:s + SEGMENT_CELLS]
        if len(chunk) < 2:
            continue
        length_m = len(chunk) * DX
        segments.append({
            "id": f"{name.split(' - ')[0]}{s // SEGMENT_CELLS + 1:02d}",
            "cells": chunk,
            "length_m": length_m,
            "value_usd": cost_per_km * length_m / 1000.0,
        })
    return {
        "name": name,
        "axis": axis,
        "width_cells": width_cells,
        "cost_per_km": cost_per_km,
        "segments": segments,
        "length_m": sum(s["length_m"] for s in segments),
        "value_usd": sum(s["value_usd"] for s in segments),
    }


def asset_mask(assets, shape):
    """Boolean mask of all road cells, for plotting."""
    mask = np.zeros(shape, dtype=bool)
    for a in assets:
        half = a["width_cells"] // 2
        for seg in a["segments"]:
            for (j, i) in seg["cells"]:
                if a["axis"] == "alongshore":
                    mask[j, max(0, i - half):i + half + 1] = True
                else:
                    mask[max(0, j - half):j + half + 1, i] = True
    return mask


if __name__ == "__main__":
    z, n, assets, x, y, deck = build_domain()
    print(f"Domain {NY} x {NX} cells at {DX:.0f} m  "
          f"({y[-1]+DX:.0f} m alongshore x {x[-1]+DX:.0f} m cross-shore)")
    print(f"Bed elevation range {z.min():+.2f} to {z.max():+.2f} m MSL\n")
    total = 0.0
    for a in assets:
        print(f"{a['name']:24s} {a['length_m']:6.0f} m  "
              f"{len(a['segments']):2d} segments  "
              f"exposure ${a['value_usd']/1e6:6.2f} M")
        total += a["value_usd"]
    print(f"{'TOTAL':24s} {'':6s}    {'':2s}           "
          f"          ${total/1e6:6.2f} M")
