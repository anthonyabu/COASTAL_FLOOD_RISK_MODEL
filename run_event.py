"""Run a single coastal surge event over the domain and return the hazard fields."""

import time
import numpy as np

from solver import LocalInertialSolver
from terrain import build_domain, DX
import damage as dmg


def run_event(surge_peak, z=None, manning=None, verbose=False,
              t_end=None, report_every=1800.0, surge_duration=None):
    """
    Drive the solver with a tide + surge boundary condition on the west edge.

    Returns a dict with the peak-depth field, wet-duration field, mass-balance audit
    and a time series of the boundary water level and inundated area.
    """
    if z is None or manning is None:
        z, manning, _, _, _, _ = build_domain()

    Td = dmg.SURGE_DURATION if surge_duration is None else surge_duration
    t_end = (Td + 6 * 3600.0) if t_end is None else t_end
    s = LocalInertialSolver(z, dx=DX, manning=manning, cfl=0.5)

    series = {"t_h": [], "swl_m": [], "wet_area_ha": [], "volume_m3": []}
    next_report = 0.0
    n_steps = 0
    t0 = time.time()

    while s.t < t_end:
        s.set_sea_level(dmg.water_level(s.t, surge_peak, Td))
        s.step()
        n_steps += 1
        if s.t >= next_report:
            wet = (s.h > 0.05).sum() * DX ** 2 / 10_000.0
            series["t_h"].append(s.t / 3600.0)
            series["swl_m"].append(dmg.water_level(s.t, surge_peak, Td))
            series["wet_area_ha"].append(wet)
            series["volume_m3"].append(float(s.h.sum()) * DX ** 2)
            next_report += report_every

    stored, vin, corr, pct = s.mass_balance()
    runtime = time.time() - t0

    if verbose:
        print(f"  surge {surge_peak:.2f} m | {n_steps:6d} steps | {runtime:5.1f} s | "
              f"mass err {pct:+.4f} % | clip {corr:.3e} m3")

    return {
        "surge_peak": surge_peak,
        "h_max": s.h_max.copy(),
        "wet_duration": s.wet_duration.copy(),
        "series": {k: np.array(v) for k, v in series.items()},
        "mass_pct_error": pct,
        "clip_correction_m3": corr,
        "n_steps": n_steps,
        "runtime_s": runtime,
    }


if __name__ == "__main__":
    z, manning, assets, x, y, deck = build_domain()
    S = dmg.surge_for_return_period(100)
    print(f"Test run: 100-year event, surge {S:.2f} m, "
          f"peak SWL {dmg.peak_still_water_level(S):.2f} m MSL")
    r = run_event(S, z, manning, verbose=True)
    print(f"  max depth anywhere {r['h_max'].max():.2f} m")
    print(f"  peak inundated area {r['series']['wet_area_ha'].max():.1f} ha")
