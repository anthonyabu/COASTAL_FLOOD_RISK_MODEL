"""
Verification of the local-inertial solver against cases with known answers.

T1  Closed basin, flat bed, initial mound of water.
    Expected: relaxes to a horizontal free surface; volume exactly conserved.

T2  Closed basin, tilted bed, initial mound.
    Expected: relaxes to a horizontal free surface (water surface elevation uniform
    over all wet cells); volume conserved.

T3  Sloping shore with the west boundary held at a fixed sea level.
    Expected: fills until the water surface equals the sea level and inflow stops;
    the shoreline sits where bed elevation = sea level.

T4  Dry-bed dam break (Ritter-type), 1D.
    Local-inertial schemes do NOT reproduce the exact Ritter solution because the
    advective term is dropped, but the front position should be the right order and
    the solution must stay bounded and monotone. This is a stability check, not an
    accuracy check.
"""

import numpy as np
from solver import LocalInertialSolver


def run(solver, t_end, sea=None):
    while solver.t < t_end:
        if sea is not None:
            solver.set_sea_level(sea)
        solver.step()


def t1_flat_basin():
    z = np.zeros((40, 40))
    s = LocalInertialSolver(z, dx=10.0, manning=0.03)
    s.h[15:25, 15:25] = 2.0
    v0 = s.h.sum() * s.dx ** 2
    run(s, 3000.0)
    v1 = s.h.sum() * s.dx ** 2
    wse = s.h + s.z
    print(f"T1 flat basin      volume err {100*(v1-v0)/v0:+8.4f} %   "
          f"WSE range {wse.max()-wse.min():.5f} m   mean depth {s.h.mean():.4f} m")


def t2_tilted_basin():
    x = np.arange(40) * 10.0
    z = np.tile(0.002 * x, (40, 1))          # 0.2 % slope
    s = LocalInertialSolver(z, dx=10.0, manning=0.03)
    s.h[:, :10] = np.maximum(1.5 - z[:, :10], 0.0)
    v0 = s.h.sum() * s.dx ** 2
    run(s, 6000.0)
    v1 = s.h.sum() * s.dx ** 2
    wet = s.h > 0.02
    wse = (s.h + s.z)[wet]
    print(f"T2 tilted basin    volume err {100*(v1-v0)/v0:+8.4f} %   "
          f"WSE range over wet cells {wse.max()-wse.min():.5f} m")


def t3_fixed_sea_level():
    x = np.arange(60) * 20.0
    z = np.tile(-1.0 + 0.0025 * x, (20, 1))   # -1.0 m at x=0 rising to +1.95 m
    sea = 0.60
    s = LocalInertialSolver(z, dx=20.0, manning=0.05)
    run(s, 40000.0, sea=sea)
    wet = s.h > 0.02
    wse = (s.h + s.z)[wet]
    # analytic shoreline: bed elevation == sea level
    x_shore_exact = (sea + 1.0) / 0.0025
    wet_cols = np.where(wet.any(axis=0))[0]
    x_shore_model = x[wet_cols[-1]]
    stored, vin, corr, pct = s.mass_balance()
    print(f"T3 fixed sea level WSE range {wse.max()-wse.min():.5f} m  "
          f"(target {sea:.2f}, model mean {wse.mean():.4f})")
    print(f"                   shoreline exact {x_shore_exact:.0f} m, "
          f"model {x_shore_model:.0f} m (cell size 20 m)")
    print(f"                   mass balance err {pct:+.4f} %, "
          f"clip correction {corr:.4f} m3 of {vin:.1f} m3 inflow")


def t4_dam_break():
    z = np.zeros((5, 200))
    s = LocalInertialSolver(z, dx=5.0, manning=0.01, cfl=0.4)
    s.h[:, :100] = 1.0
    v0 = s.h.sum() * s.dx ** 2
    run(s, 30.0)
    v1 = s.h.sum() * s.dx ** 2
    front = np.where(s.h[2, :] > 0.01)[0]
    x_front = front[-1] * 5.0
    # Ritter dry-bed front celerity = 2*sqrt(g*h0); local-inertial underpredicts this
    x_ritter = 500.0 + 2.0 * np.sqrt(9.81 * 1.0) * 30.0
    print(f"T4 dam break       volume err {100*(v1-v0)/v0:+8.4f} %   "
          f"front at {x_front:.0f} m (Ritter {x_ritter:.0f} m)")
    print(f"                   depth range {s.h.min():.4f} to {s.h.max():.4f} m "
          f"(must stay within 0 to 1)")


if __name__ == "__main__":
    print("Solver verification\n" + "-" * 70)
    t1_flat_basin()
    t2_tilted_basin()
    t3_fixed_sea_level()
    t4_dam_break()
