"""
When is a static (bathtub) inundation map an acceptable substitute for a dynamic model?

The peak still water level is held fixed at the 100-year value. Only the DURATION of
the surge is varied. A bathtub map has no time dimension, so it returns the same answer
for every duration - the equilibrium answer. The dynamic model does not, because
filling a floodplain takes time and friction resists it.

The question this answers is a practical screening question: a state or municipal
agency doing a first-pass asset prioritisation will usually be handed static hazard
maps. This quantifies the conditions under which those maps are adequate and the
conditions under which they materially overstate exposure.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from terrain import build_domain, DX
from run_event import run_event
from analysis import bathtub
import damage as dmg

DURATIONS_H = [1.0, 2.0, 3.0, 6.0, 12.0, 24.0]
# Paths resolve relative to this file, so the code runs from any directory on any
# machine. Override the output location with the FLOOD_OUT environment variable.
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.environ.get("FLOOD_OUT", os.path.join(HERE, "outputs"))
os.makedirs(OUT, exist_ok=True)


def main():
    z, manning, assets, x, y, deck = build_domain()
    S = dmg.surge_for_return_period(100)
    swl = dmg.peak_still_water_level(S)

    bt = bathtub(z, swl)
    bt_area = (bt > 0.05).sum() * DX ** 2 / 1e4
    bt_assess = dmg.assess_assets(assets, bt - z + z * 0, np.zeros_like(z),
                                  z=np.zeros_like(z), deck=deck - z + deck * 0) \
        if False else None
    # bathtub damage: treat the bathtub water surface as the peak WSE
    bt_wse = np.where(bt > 0, swl, z)
    bt_assess = dmg.assess_assets(assets, bt_wse - z, np.zeros_like(z),
                                  z=z, deck=deck)
    bt_damage = sum(s["damage_usd"] for s in bt_assess)

    print(f"100-year peak still water level: {swl:.2f} m MSL")
    print(f"Static bathtub:  inundated {bt_area:6.1f} ha   "
          f"road damage ${bt_damage/1e6:.2f} M   (no duration dependence)\n")

    print(f"{'surge dur (h)':>14} {'area (ha)':>11} {'% of bathtub':>13} "
          f"{'damage ($M)':>12} {'% of bathtub':>13}")
    rows = []
    for Th in DURATIONS_H:
        r = run_event(S, z, manning, surge_duration=Th * 3600.0)
        area = (r["h_max"] > 0.05).sum() * DX ** 2 / 1e4
        a_ = dmg.assess_assets(assets, r["h_max"], r["wet_duration"], z=z, deck=deck)
        dmg_usd = sum(s["damage_usd"] for s in a_)
        rows.append((Th, area, dmg_usd, r["mass_pct_error"]))
        print(f"{Th:14.1f} {area:11.1f} {100*area/bt_area:12.1f}% "
              f"{dmg_usd/1e6:12.2f} {100*dmg_usd/max(bt_damage,1):12.1f}%")

    # ---- marsh roughness sensitivity ----
    # Salt marsh vegetation attenuates surge propagation. A static bathtub map cannot
    # represent this at all: it has no momentum equation, so roughness does not appear
    # anywhere in it. Any credible assessment of marsh as natural coastal protection
    # therefore requires a dynamic model.
    # Roughness is swept at several surge durations, not just the long one. Friction
    # controls the RATE of filling, not the equilibrium level, so it can only change
    # the answer while the event is still short relative to the basin filling time.
    # Sweeping roughness at 12 h alone would wrongly suggest roughness never matters.
    N_VALUES = [0.03, 0.06, 0.10, 0.18, 0.30]
    ROUGH_DURATIONS_H = [0.5, 1.0, 3.0]
    rough = {}
    print(f"\nRoughness x duration  (damage, USD million; "
          f"static bathtub = ${bt_damage/1e6:.2f} M for all cells)")
    print(f"{'marsh n':>10}" + "".join(f"{f'{t:g} h':>10}" for t in ROUGH_DURATIONS_H))
    for n_marsh in N_VALUES:
        m2 = manning.copy()
        m2[:, :20] = np.where(m2[:, :20] > 0.019, n_marsh, m2[:, :20])
        row = []
        for Th_ in ROUGH_DURATIONS_H:
            r = run_event(S, z, m2, surge_duration=Th_ * 3600.0)
            a_ = dmg.assess_assets(assets, r["h_max"], r["wet_duration"],
                                   z=z, deck=deck)
            row.append(sum(s["damage_usd"] for s in a_))
        rough[n_marsh] = row
        print(f"{n_marsh:10.2f}" + "".join(f"{v/1e6:10.2f}" for v in row))
    rough_rows = [(n, 0.0, rough[n][-1]) for n in N_VALUES]

    Th = np.array([r[0] for r in rows])
    A = np.array([r[1] for r in rows])
    D = np.array([r[2] for r in rows])
    Rn = np.array([r[0] for r in rough_rows])
    RD = np.array([r[2] for r in rough_rows])

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
    axes[0].axhline(bt_area, color="firebrick", ls="--", lw=2,
                    label="static bathtub (equilibrium)")
    axes[0].plot(Th, A, marker="o", lw=2, color="#08306b", label="dynamic model")
    axes[0].set_xlabel("Surge duration (hours)")
    axes[0].set_ylabel("Peak inundated area (ha)")
    axes[0].set_title("Inundation extent vs surge duration\n"
                      "(peak water level held constant at the 100-year value)",
                      fontsize=10)
    axes[0].legend(fontsize=9)
    axes[0].grid(alpha=0.3)

    axes[1].axhline(bt_damage / 1e6, color="firebrick", ls="--", lw=2,
                    label="static bathtub")
    axes[1].plot(Th, D / 1e6, marker="o", lw=2, color="#08306b", label="dynamic model")
    axes[1].set_xlabel("Surge duration (hours)")
    axes[1].set_ylabel("Road damage (USD million)")
    axes[1].set_title("Damage vs surge duration", fontsize=10)
    axes[1].legend(fontsize=9)
    axes[1].grid(alpha=0.3)

    axes[2].axhline(bt_damage / 1e6, color="firebrick", ls="--", lw=2,
                    label="static bathtub (roughness-independent)")
    cols = ["#238b45", "#41ab5d", "#a1d99b"]
    for k, Th_ in enumerate(ROUGH_DURATIONS_H):
        vals = np.array([rough[n][k] for n in N_VALUES])
        axes[2].plot(N_VALUES, vals / 1e6, marker="o", lw=2, color=cols[k],
                     label=f"dynamic, {Th_:g} h surge")
    axes[2].set_xlabel("Marsh Manning's n (first 400 m from shore)")
    axes[2].set_ylabel("Road damage (USD million)")
    axes[2].set_title("Damage vs marsh roughness\n"
                      "a static map cannot represent this axis at all", fontsize=10)
    axes[2].legend(fontsize=8)
    axes[2].grid(alpha=0.3)

    fig.suptitle("When does a static hazard map overstate exposure?", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(f"{OUT}/fig7_duration_sensitivity.png", dpi=160)
    plt.close(fig)

    with open(f"{OUT}/duration_sensitivity.csv", "w") as f:
        f.write("surge_duration_h,inundated_area_ha,pct_of_bathtub_area,"
                "damage_usd,pct_of_bathtub_damage,mass_balance_pct_error\n")
        for Th_, A_, D_, m_ in rows:
            f.write(f"{Th_},{A_:.1f},{100*A_/bt_area:.1f},{D_:.0f},"
                    f"{100*D_/max(bt_damage,1):.1f},{m_:.6f}\n")
        f.write(f"bathtub,{bt_area:.1f},100.0,{bt_damage:.0f},100.0,n/a\n")

    print(f"\nWritten to {OUT}/fig7_duration_sensitivity.png")


if __name__ == "__main__":
    main()
