"""
Coastal flood risk assessment for transportation assets.

Pipeline:
  1. Build the domain and asset inventory.
  2. Run a dynamic inundation simulation for a set of return periods.
  3. Compare the dynamic result against a connectivity-constrained bathtub map,
     which is the method this kind of study is usually screened with.
  4. Apply a depth-damage function per road segment.
  5. Integrate damage over exceedance probability to get expected annual damage.
  6. Write figures and a results table.
"""

import hashlib
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy import ndimage

from terrain import build_domain, asset_mask, DX
from run_event import run_event
import damage as dmg

RETURN_PERIODS = [1.5, 2, 3, 5, 10, 15, 25, 50, 100, 200, 350, 500]
# Paths resolve relative to this file, so the code runs from any directory on any
# machine. Override the output location with the FLOOD_OUT environment variable.
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.environ.get("FLOOD_OUT", os.path.join(HERE, "outputs"))
os.makedirs(OUT, exist_ok=True)

DEPTH_CMAP = LinearSegmentedColormap.from_list(
    "depth", ["#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b", "#041e42"])


# ---------------------------------------------------------------------------
def bathtub(z, swl):
    """
    Connectivity-constrained bathtub inundation: cells below the still water level
    that are hydraulically connected to the sea. This is the standard static
    screening method. It has no timing, no friction and no volume constraint, so it
    represents the equilibrium upper bound rather than what a finite-duration event
    can actually deliver.
    """
    below = z < swl
    labels, _ = ndimage.label(below)
    sea_labels = np.unique(labels[:, 0][below[:, 0]])
    connected = np.isin(labels, sea_labels[sea_labels > 0])
    return np.where(connected, swl - z, 0.0)


# ---------------------------------------------------------------------------
def main():
    z, manning, assets, x, y, deck = build_domain()
    roads = asset_mask(assets, z.shape)
    extent = [0, x[-1] + DX, 0, y[-1] + DX]

    print("Domain and assets")
    for a in assets:
        print(f"  {a['name']:24s} {a['length_m']:6.0f} m   "
              f"exposure ${a['value_usd']/1e6:5.2f} M")
    total_exposure = sum(a["value_usd"] for a in assets)
    print(f"  {'TOTAL EXPOSURE':24s} {'':6s}   ${total_exposure/1e6:5.2f} M\n")

    # ---- run the ensemble (cached: the hydraulics do not depend on the damage curve) ----
    cache = os.path.join(HERE, "_ensemble_cache.npz")

    # The cache stores simulated water depth fields. Those fields depend on the
    # terrain, the roughness, the surge distribution and the tide, so the cache is
    # only valid for the exact parameters that produced it. Keying on the filename
    # alone silently mixes stale hydraulics with new terrain - which produces a
    # result that is neither the old one nor the new one, and looks plausible.
    # The fingerprint below makes that failure loud instead of silent.
    fingerprint = hashlib.sha256(b"".join([
        z.tobytes(), manning.tobytes(), deck[~np.isnan(deck)].tobytes(),
        np.array([dmg.GEV_LOC, dmg.GEV_SCALE, dmg.GEV_SHAPE,
                  dmg.TIDE_AMPLITUDE, dmg.TIDE_PERIOD, dmg.SURGE_DURATION,
                  dmg.EVENT_DURATION]).tobytes(),
        np.array(RETURN_PERIODS, dtype=float).tobytes(),
    ])).hexdigest()[:16]

    runs = {}
    reuse = False
    if os.path.exists(cache):
        blob = np.load(cache, allow_pickle=True)
        stored = str(blob["fingerprint"]) if "fingerprint" in blob else "<none>"
        if stored == fingerprint:
            reuse = True
        else:
            print(f"Cache found but parameters have changed "
                  f"(cache {stored}, current {fingerprint}).")
            print("Discarding it and re-running the hydraulics.")
            os.remove(cache)

    if reuse:
        print(f"Loading cached simulations (fingerprint {fingerprint})")
        for T in RETURN_PERIODS:
            runs[T] = {k: blob[f"{T}_{k}"] for k in
                       ("h_max", "wet_duration", "bathtub")}
            runs[T].update(blob[f"{T}_meta"].item())
    else:
        print(f"Running {len(RETURN_PERIODS)} events")
        store = {}
        for T in RETURN_PERIODS:
            S = dmg.surge_for_return_period(T)
            print(f"  T = {T:5.1f} yr  surge {S:.2f} m  "
                  f"peak SWL {dmg.peak_still_water_level(S):.2f} m", end="", flush=True)
            r = run_event(S, z, manning)
            r["T"] = T
            r["swl_peak"] = dmg.peak_still_water_level(S)
            r["bathtub"] = bathtub(z, r["swl_peak"])
            runs[T] = r
            print(f"   mass err {r['mass_pct_error']:+.4f} %")
            for k in ("h_max", "wet_duration", "bathtub"):
                store[f"{T}_{k}"] = r[k]
            store[f"{T}_meta"] = {
                "T": T, "swl_peak": r["swl_peak"], "surge_peak": r["surge_peak"],
                "series": r["series"], "mass_pct_error": r["mass_pct_error"],
                "clip_correction_m3": r["clip_correction_m3"],
                "n_steps": r["n_steps"], "runtime_s": r["runtime_s"],
            }
        store["fingerprint"] = fingerprint
        np.savez_compressed(cache, **store)

    # ---- damage curve sensitivity ----
    Ts = np.array(RETURN_PERIODS, dtype=float)
    print("\nDamage curve sensitivity (identical hydraulics, different damage curve)")
    sens = {}
    for variant in dmg.DD_VARIANTS:
        Dv = []
        for T in RETURN_PERIODS:
            a_ = dmg.assess_assets(assets, runs[T]["h_max"], runs[T]["wet_duration"],
                                   z=z, deck=deck, variant=variant)
            Dv.append(sum(s["damage_usd"] for s in a_))
        e = dmg.expected_annual_damage(Ts, np.array(Dv))
        sens[variant] = {"damages": Dv, "ead": e}
        print(f"  {variant:14s}  EAD ${e['ead_usd']:>10,.0f} /yr   "
              f"({100*e['ead_usd']/total_exposure:5.2f} % of exposure per year)   "
              f"100-yr damage ${Dv[RETURN_PERIODS.index(100)]/1e6:.2f} M")
    r_lo = sens["threshold"]["ead"]["ead_usd"]
    r_hi = sens["no_threshold"]["ead"]["ead_usd"]
    print(f"  -> EAD changes by a factor of {r_hi/r_lo:.2f} from this choice alone, "
          f"while 100-yr damage changes by only "
          f"{sens['no_threshold']['damages'][RETURN_PERIODS.index(100)]/sens['threshold']['damages'][RETURN_PERIODS.index(100)]:.2f}x")

    # ---- attach the default-variant assessment ----
    for T in RETURN_PERIODS:
        runs[T]["assessment"] = dmg.assess_assets(
            assets, runs[T]["h_max"], runs[T]["wet_duration"], z=z, deck=deck)
        runs[T]["total_damage"] = sum(s["damage_usd"] for s in runs[T]["assessment"])

    # ---- mass balance audit ----
    worst_mass = max(abs(r["mass_pct_error"]) for r in runs.values())
    worst_clip = max(r["clip_correction_m3"] for r in runs.values())
    print(f"\nMass balance audit: worst error {worst_mass:.5f} %, "
          f"worst clipping correction {worst_clip:.3e} m3")

    # ---- expected annual damage ----
    D = np.array([runs[T]["total_damage"] for T in RETURN_PERIODS])
    ead = dmg.expected_annual_damage(Ts, D)
    print(f"\nExpected annual damage  ${ead['ead_usd']:,.0f} /yr")
    print(f"  core integral         ${ead['ead_core_usd']:,.0f}")
    print(f"  rare tail (T > {Ts.max():.0f} yr)  ${ead['tail_usd']:,.0f}")
    print(f"  frequent tail (T < {Ts.min():.1f} yr) ${ead['frequent_tail_usd']:,.0f}")

    # per-asset EAD
    per_asset_ead = {}
    for a in assets:
        Da = np.array([
            sum(s["damage_usd"] for s in runs[T]["assessment"] if s["asset"] == a["name"])
            for T in RETURN_PERIODS])
        per_asset_ead[a["name"]] = dmg.expected_annual_damage(Ts, Da)["ead_usd"]

    print("\nEAD by asset")
    for name, v in sorted(per_asset_ead.items(), key=lambda kv: -kv[1]):
        exp = next(a["value_usd"] for a in assets if a["name"] == name)
        print(f"  {name:24s} ${v:>10,.0f} /yr   "
              f"({100*v/exp:.2f} % of asset value per year)")

    # ---- figures ----
    fig_domain(z, roads, assets, extent)
    fig_hazard(runs, roads, extent, [10, 100, 500])
    fig_bathtub_compare(runs[100], z, roads, extent)
    fig_hydrograph(runs, [10, 100, 500])
    fig_damage_curves(runs, assets, Ts, D, ead, per_asset_ead, sens)
    fig_segment_map(runs[100], assets, z, extent)

    # ---- results table ----
    write_results(runs, assets, ead, per_asset_ead, total_exposure, sens)
    print(f"\nOutputs written to {OUT}")


# ---------------------------------------------------------------------------
def fig_domain(z, roads, assets, extent):
    fig, ax = plt.subplots(figsize=(11, 5.2))
    norm = TwoSlopeNorm(vmin=z.min(), vcenter=0.0, vmax=z.max())
    im = ax.imshow(z, origin="lower", extent=extent, cmap="terrain", norm=norm,
                   aspect="auto")
    ax.contour(z, levels=[0], origin="lower", extent=extent,
               colors="k", linewidths=0.8, linestyles="--")
    ov = np.ma.masked_where(~roads, np.ones_like(z))
    ax.imshow(ov, origin="lower", extent=extent, cmap="autumn_r",
              alpha=0.95, aspect="auto", vmin=0, vmax=1)
    ax.annotate("A  Coastal Highway\ncrest 2.35 m", (700, 1080), color="k",
                fontsize=9, ha="center",
                bbox=dict(fc="white", alpha=0.85, ec="none", pad=2))
    ax.annotate("B  Access Causeway\ncrest 2.10 m", (430, 190), color="k", fontsize=9,
                ha="center", bbox=dict(fc="white", alpha=0.85, ec="none", pad=2))
    ax.annotate("C  Inland Arterial\ncrest 2.55 m", (1500, 1080), color="k", fontsize=9,
                ha="center", bbox=dict(fc="white", alpha=0.85, ec="none", pad=2))
    ax.annotate("tidal creek\n(bridged opening)", (760, 600), color="k", fontsize=8,
                ha="left", va="center",
                bbox=dict(fc="white", alpha=0.85, ec="none", pad=2))
    ax.set_xlabel("Cross-shore distance from open boundary (m)")
    ax.set_ylabel("Alongshore distance (m)")
    ax.set_title("Model domain: bed elevation, transportation assets and drainage\n"
                 "SYNTHETIC TERRAIN - not a surveyed site", fontsize=11)
    fig.colorbar(im, ax=ax, label="Bed elevation (m MSL)", pad=0.02)
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig1_domain.png", dpi=160)
    plt.close(fig)


def fig_hazard(runs, roads, extent, Ts):
    fig, axes = plt.subplots(len(Ts), 1, figsize=(11, 11), sharex=True)
    vmax = max(runs[T]["h_max"].max() for T in Ts)
    for ax, T in zip(axes, Ts):
        r = runs[T]
        d = np.ma.masked_where(r["h_max"] < 0.05, r["h_max"])
        ax.imshow(np.zeros_like(r["h_max"]), origin="lower", extent=extent,
                  cmap="Greys", vmin=0, vmax=1, aspect="auto")
        im = ax.imshow(d, origin="lower", extent=extent, cmap=DEPTH_CMAP,
                       vmin=0, vmax=vmax, aspect="auto")
        ov = np.ma.masked_where(~roads, np.ones_like(r["h_max"]))
        ax.imshow(ov, origin="lower", extent=extent, cmap="autumn_r",
                  alpha=0.9, aspect="auto", vmin=0, vmax=1)
        ax.set_ylabel("Alongshore (m)")
        ax.set_title(f"{T:.0f}-year event   peak SWL {r['swl_peak']:.2f} m MSL   "
                     f"road damage ${r['total_damage']/1e6:.2f} M", fontsize=10)
        fig.colorbar(im, ax=ax, label="Peak depth (m)", pad=0.02)
    axes[-1].set_xlabel("Cross-shore distance (m)")
    fig.suptitle("Modelled peak inundation depth by return period", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(f"{OUT}/fig2_hazard_maps.png", dpi=160)
    plt.close(fig)


def fig_bathtub_compare(r, z, roads, extent):
    diff = r["bathtub"] - r["h_max"]
    fig, axes = plt.subplots(3, 1, figsize=(11, 11), sharex=True)

    for ax, field, title, cmap, kw in [
        (axes[0], r["bathtub"], "Connectivity-constrained bathtub (static)",
         DEPTH_CMAP, dict(vmin=0, vmax=max(r["bathtub"].max(), r["h_max"].max()))),
        (axes[1], r["h_max"], "Dynamic local-inertial model (this work)",
         DEPTH_CMAP, dict(vmin=0, vmax=max(r["bathtub"].max(), r["h_max"].max()))),
    ]:
        d = np.ma.masked_where(field < 0.05, field)
        ax.imshow(np.zeros_like(z), origin="lower", extent=extent, cmap="Greys",
                  vmin=0, vmax=1, aspect="auto")
        im = ax.imshow(d, origin="lower", extent=extent, cmap=cmap, aspect="auto", **kw)
        ov = np.ma.masked_where(~roads, np.ones_like(z))
        ax.imshow(ov, origin="lower", extent=extent, cmap="autumn_r", alpha=0.9,
                  aspect="auto", vmin=0, vmax=1)
        ax.set_title(title, fontsize=10)
        ax.set_ylabel("Alongshore (m)")
        fig.colorbar(im, ax=ax, label="Peak depth (m)", pad=0.02)

    dm = np.ma.masked_where(np.abs(diff) < 0.02, diff)
    im = axes[2].imshow(dm, origin="lower", extent=extent, cmap="RdBu_r",
                        norm=TwoSlopeNorm(vcenter=0, vmin=-abs(diff).max(),
                                          vmax=abs(diff).max()), aspect="auto")
    axes[2].set_title("Bathtub minus dynamic  (red = static method overpredicts depth)",
                      fontsize=10)
    axes[2].set_ylabel("Alongshore (m)")
    axes[2].set_xlabel("Cross-shore distance (m)")
    fig.colorbar(im, ax=axes[2], label="Depth difference (m)", pad=0.02)

    over = (diff > 0.05).sum() / diff.size * 100
    fig.suptitle(f"Static vs dynamic inundation, 100-year event\n"
                 f"Static method overpredicts depth by >0.05 m over {over:.0f} % "
                 f"of the domain", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(f"{OUT}/fig3_bathtub_vs_dynamic.png", dpi=160)
    plt.close(fig)


def fig_hydrograph(runs, Ts):
    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    colors = ["#4292c6", "#2171b5", "#08306b"]
    for T, c in zip(Ts, colors):
        s = runs[T]["series"]
        axes[0].plot(s["t_h"], s["swl_m"], color=c, lw=2, label=f"{T:.0f}-yr")
        axes[1].plot(s["t_h"], s["wet_area_ha"], color=c, lw=2, label=f"{T:.0f}-yr")
    for crest, lbl, ls in [(2.35, "Coastal Highway crest", "--"),
                           (2.10, "Causeway crest", ":"),
                           (2.55, "Inland Arterial crest", "-.")]:
        axes[0].axhline(crest, color="firebrick", ls=ls, lw=1.1, label=lbl)
    axes[0].set_ylabel("Boundary still water level (m MSL)")
    axes[0].legend(fontsize=8, ncol=2)
    axes[0].grid(alpha=0.3)
    axes[0].set_title("Forcing: tide + surge at the open boundary", fontsize=10)
    axes[1].set_ylabel("Inundated area, depth > 0.05 m (ha)")
    axes[1].set_xlabel("Time from start of event (hours)")
    axes[1].grid(alpha=0.3)
    axes[1].legend(fontsize=8)
    axes[1].set_title("Response: inundated area", fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig4_hydrograph.png", dpi=160)
    plt.close(fig)


def fig_damage_curves(runs, assets, Ts, D, ead, per_asset_ead, sens):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))

    # depth-damage curves, both variants
    dd = np.linspace(0, 3.2, 300)
    style = {"threshold": ("firebrick", "-"), "no_threshold": ("#2171b5", "--")}
    for variant in dmg.DD_VARIANTS:
        c, ls = style[variant]
        e = sens[variant]["ead"]["ead_usd"]
        axes[0].plot(dd, dmg.damage_ratio(dd, variant=variant), color=c, ls=ls, lw=2,
                     label=f"{variant}\nEAD ${e/1e3:,.0f} k/yr")
        d_, r_ = dmg.DD_VARIANTS[variant]
        axes[0].scatter(d_, r_, color=c, s=20, zorder=3)
    axes[0].set_xlabel("Peak inundation depth (m)")
    axes[0].set_ylabel("Damage ratio (fraction of replacement value)")
    axes[0].set_title("Road depth-damage function\nILLUSTRATIVE - not sourced curves",
                      fontsize=10)
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.3)

    # damage vs return period
    for a in assets:
        Da = np.array([
            sum(s["damage_usd"] for s in runs[T]["assessment"] if s["asset"] == a["name"])
            for T in Ts])
        axes[1].semilogx(Ts, Da / 1e6, marker="o", ms=4, lw=1.8,
                         label=a["name"].split(" - ")[1])
    axes[1].semilogx(Ts, D / 1e6, marker="s", ms=5, lw=2.4, color="k", label="Total")
    axes[1].set_xlabel("Return period (years)")
    axes[1].set_ylabel("Damage (USD million)")
    axes[1].set_title("Damage vs return period", fontsize=10)
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.3, which="both")

    # exceedance / EAD
    p = 1.0 / Ts
    axes[2].plot(p, D / 1e6, marker="o", ms=4, lw=2, color="#08306b")
    axes[2].fill_between(p, 0, D / 1e6, alpha=0.25, color="#6baed6")
    axes[2].set_xlabel("Annual exceedance probability")
    axes[2].set_ylabel("Damage (USD million)")
    axes[2].set_title(f"Risk curve\nshaded area = EAD = "
                      f"${ead['ead_usd']/1e3:,.0f} k / yr", fontsize=10)
    axes[2].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig5_damage_and_risk.png", dpi=160)
    plt.close(fig)


def fig_segment_map(r, assets, z, extent):
    fig, ax = plt.subplots(figsize=(11, 5.2))
    ax.imshow(z, origin="lower", extent=extent, cmap="Greys_r", alpha=0.35,
              aspect="auto")
    lookup = {(s["asset"], s["segment"]): s for s in r["assessment"]}
    cmap = plt.get_cmap("YlOrRd")
    for a in assets:
        for seg in a["segments"]:
            rec = lookup[(a["name"], seg["id"])]
            xs = [i * DX for (j, i) in seg["cells"]]
            ys = [j * DX for (j, i) in seg["cells"]]
            ax.plot(xs, ys, lw=6, color=cmap(rec["damage_ratio"]),
                    solid_capstyle="butt")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 1))
    fig.colorbar(sm, ax=ax, label="Damage ratio", pad=0.02)
    ax.set_xlabel("Cross-shore distance (m)")
    ax.set_ylabel("Alongshore distance (m)")
    ax.set_title("Segment-level damage, 100-year event\n"
                 "Each segment is 100 m; damage driven by the worst cell on the segment",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig6_segment_damage.png", dpi=160)
    plt.close(fig)


def write_results(runs, assets, ead, per_asset_ead, total_exposure, sens):
    lines = ["return_period_yr,peak_swl_m,inundated_area_ha,max_depth_m,"
             "damage_usd,damage_pct_of_exposure,mass_balance_pct_error"]
    for T, r in runs.items():
        area = (r["h_max"] > 0.05).sum() * DX ** 2 / 10_000.0
        lines.append(f"{T},{r['swl_peak']:.3f},{area:.1f},{r['h_max'].max():.3f},"
                     f"{r['total_damage']:.0f},"
                     f"{100*r['total_damage']/total_exposure:.2f},"
                     f"{r['mass_pct_error']:.6f}")
    with open(f"{OUT}/results_by_return_period.csv", "w") as f:
        f.write("\n".join(lines) + "\n")

    seg_lines = ["return_period_yr,asset,segment,length_m,value_usd,"
                 "peak_depth_m,duration_h,damage_ratio,damage_usd"]
    for T, r in runs.items():
        for s in r["assessment"]:
            seg_lines.append(
                f"{T},{s['asset']},{s['segment']},{s['length_m']:.0f},"
                f"{s['value_usd']:.0f},{s['depth_peak_m']:.3f},{s['duration_h']:.2f},"
                f"{s['damage_ratio']:.4f},{s['damage_usd']:.0f}")
    with open(f"{OUT}/results_by_segment.csv", "w") as f:
        f.write("\n".join(seg_lines) + "\n")

    summary = {
        "total_exposure_usd": total_exposure,
        "expected_annual_damage_usd": ead["ead_usd"],
        "ead_components": ead,
        "ead_by_asset_usd": per_asset_ead,
        "damage_curve_sensitivity": {
            k: {"ead_usd": v["ead"]["ead_usd"],
                "damage_by_return_period_usd": v["damages"]}
            for k, v in sens.items()},
        "worst_mass_balance_pct_error":
            max(abs(r["mass_pct_error"]) for r in runs.values()),
        "caveats": [
            "Terrain, road network and unit costs are synthetic.",
            "Surge distribution (Gumbel) is illustrative, not gauge-fitted.",
            "Depth-damage curve is illustrative; shape is plausible, ordinates are not sourced.",
            "Tide and surge peaks assumed coincident (conservative joint-occurrence case).",
            "Wave setup, wave overtopping, rainfall and groundwater are not modelled.",
            "Local-inertial scheme omits advective acceleration; valid for floodplain "
            "flow, not for supercritical or strongly advective flow.",
        ],
    }
    with open(f"{OUT}/summary.json", "w") as f:
        json.dump(summary, f, indent=2)


if __name__ == "__main__":
    main()
