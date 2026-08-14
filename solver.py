"""
2D inundation solver - local inertial approximation of the shallow water equations.

Formulation follows the local-inertial (inertial-diffusive) scheme, which drops the
advective acceleration term from the full SWE momentum equation. This is the standard
approach for floodplain / coastal overland flow where advection is small relative to
gravity, friction and local acceleration. It is the scheme underlying LISFLOOD-FP and
several operational flood forecasting models.

Momentum (per unit width), semi-implicit in friction:

    q^(n+1) = [ q^n  -  g * h_flow * dt * d(h+z)/dx ]
              --------------------------------------------------
              [ 1  +  g * dt * n^2 * |q^n| / h_flow^(7/3) ]

Continuity:

    dh/dt = ( sum q_in - sum q_out ) / dx

Discretisation: staggered (Arakawa C) grid. Depths and bed elevations at cell centres,
unit-width discharges at cell faces. This avoids the chequerboard pressure oscillation
that co-located variables produce.

Stability: adaptive timestep from the gravity-wave CFL condition,
    dt = alpha * dx / sqrt(g * h_max),  alpha < 1
"""

import numpy as np

G = 9.81


class LocalInertialSolver:
    """
    Parameters
    ----------
    z : (ny, nx) float array
        Bed elevation, metres relative to datum (mean sea level).
    dx : float
        Cell size, metres. Square cells assumed.
    manning : float or (ny, nx) array
        Manning's n roughness coefficient.
    cfl : float
        CFL safety factor for the adaptive timestep.
    depth_tol : float
        Flux depth below which a face is treated as dry (metres). Guards the
        h_flow^(7/3) division in the friction term.
    """

    def __init__(self, z, dx, manning=0.06, cfl=0.5, depth_tol=1e-3,
                 dt_max=10.0, dt_min=0.01):
        self.z = np.asarray(z, dtype=float)
        self.ny, self.nx = self.z.shape
        self.dx = float(dx)
        self.n = np.broadcast_to(np.asarray(manning, dtype=float), self.z.shape).copy()
        self.cfl = float(cfl)
        self.depth_tol = float(depth_tol)
        self.dt_max = float(dt_max)
        self.dt_min = float(dt_min)

        # State
        self.h = np.zeros_like(self.z)          # water depth at cell centres
        self.qx = np.zeros((self.ny, self.nx + 1))  # unit discharge, x-faces (m2/s)
        self.qy = np.zeros((self.ny + 1, self.nx))  # unit discharge, y-faces (m2/s)

        # Diagnostics
        self.t = 0.0
        self.h_max = np.zeros_like(self.z)      # running max depth (hazard map)
        self.wet_duration = np.zeros_like(self.z)  # seconds above wet_threshold
        self.volume_in = 0.0                    # cumulative inflow through boundaries
        self.mass_correction = 0.0              # volume added by negative-depth clipping

        # Boundary: west edge sea level (None = closed)
        self._sea_level = None

    # ------------------------------------------------------------------
    def set_sea_level(self, eta):
        """Set water surface elevation (m) of the open sea ghost cells on the west edge."""
        self._sea_level = float(eta)

    # ------------------------------------------------------------------
    def _face_flux(self, q, h_a, z_a, h_b, z_b, n_face):
        """
        Local-inertial momentum update on a face between cell a (upstream side in -index
        direction) and cell b.

        h_flow is the depth available for conveyance across the face: the difference
        between the higher water surface and the higher bed. This is what makes the
        scheme handle wet/dry fronts and flow over embankments correctly.
        """
        wse_a = h_a + z_a
        wse_b = h_b + z_b

        h_flow = np.maximum(wse_a, wse_b) - np.maximum(z_a, z_b)
        wet = h_flow > self.depth_tol

        # Safe denominator base - only used where wet
        h_safe = np.where(wet, h_flow, 1.0)

        slope = (wse_b - wse_a) / self.dx

        num = q - G * h_safe * self.dt * slope
        den = 1.0 + G * self.dt * n_face ** 2 * np.abs(q) / h_safe ** (7.0 / 3.0)

        q_new = np.where(wet, num / den, 0.0)
        return q_new

    # ------------------------------------------------------------------
    def _compute_dt(self):
        h_max = self.h.max()
        if h_max <= self.depth_tol:
            return self.dt_max
        dt = self.cfl * self.dx / np.sqrt(G * h_max)
        return float(np.clip(dt, self.dt_min, self.dt_max))

    # ------------------------------------------------------------------
    def step(self, dt=None, wet_threshold=0.05):
        """Advance one timestep. Returns the timestep actually used."""
        self.dt = self._compute_dt() if dt is None else float(dt)

        # ---- interior x-faces (between column i-1 and column i) ----
        n_face_x = 0.5 * (self.n[:, :-1] + self.n[:, 1:])
        self.qx[:, 1:-1] = self._face_flux(
            self.qx[:, 1:-1],
            self.h[:, :-1], self.z[:, :-1],
            self.h[:, 1:],  self.z[:, 1:],
            n_face_x,
        )

        # ---- interior y-faces ----
        n_face_y = 0.5 * (self.n[:-1, :] + self.n[1:, :])
        self.qy[1:-1, :] = self._face_flux(
            self.qy[1:-1, :],
            self.h[:-1, :], self.z[:-1, :],
            self.h[1:, :],  self.z[1:, :],
            n_face_y,
        )

        # ---- west open boundary: ghost cell held at the sea level ----
        if self._sea_level is not None:
            z_ghost = self.z[:, 0]
            h_ghost = np.maximum(self._sea_level - z_ghost, 0.0)
            self.qx[:, 0] = self._face_flux(
                self.qx[:, 0],
                h_ghost, z_ghost,
                self.h[:, 0], self.z[:, 0],
                self.n[:, 0],
            )
        else:
            self.qx[:, 0] = 0.0

        # East / north / south edges closed (qx[:, -1], qy[0, :], qy[-1, :] stay zero)

        # ---- volume-based flux limiter ----
        # A cell cannot export more water than it holds. Compute each cell's total
        # outflow over this timestep, and if it exceeds the stored volume, scale every
        # outgoing face flux from that cell by the same factor. Scaling per source cell
        # (rather than clipping depths afterwards) keeps the scheme exactly conservative
        # and guarantees h >= 0 without adding spurious mass.
        outflow = (
            np.maximum(self.qx[:, 1:], 0.0)        # out through the east face
            + np.maximum(-self.qx[:, :-1], 0.0)    # out through the west face
            + np.maximum(self.qy[1:, :], 0.0)      # out through the south face
            + np.maximum(-self.qy[:-1, :], 0.0)    # out through the north face
        ) * self.dx * self.dt

        available = self.h * self.dx ** 2
        factor = np.where(outflow > 1e-12,
                          np.minimum(1.0, available / np.maximum(outflow, 1e-12)),
                          1.0)

        # Each face is scaled by the factor of whichever cell the water is leaving.
        self.qx[:, 1:-1] *= np.where(self.qx[:, 1:-1] > 0, factor[:, :-1], factor[:, 1:])
        self.qy[1:-1, :] *= np.where(self.qy[1:-1, :] > 0, factor[:-1, :], factor[1:, :])
        # West boundary face: inflow from the sea is unlimited, outflow is limited.
        self.qx[:, 0] = np.where(self.qx[:, 0] < 0,
                                 self.qx[:, 0] * factor[:, 0],
                                 self.qx[:, 0])

        # ---- continuity ----
        dh = (self.dt / self.dx) * (
            (self.qx[:, :-1] - self.qx[:, 1:])
            + (self.qy[:-1, :] - self.qy[1:, :])
        )
        self.h += dh

        # Track boundary inflow volume for the mass audit
        self.volume_in += self.qx[:, 0].sum() * self.dx * self.dt

        # With the limiter active this should never fire. Kept as an assertion-style
        # diagnostic: any non-zero value here means the limiter has a hole in it.
        neg = self.h < 0.0
        if neg.any():
            self.mass_correction += float(-self.h[neg].sum()) * self.dx ** 2
            self.h[neg] = 0.0

        # ---- diagnostics ----
        np.maximum(self.h_max, self.h, out=self.h_max)
        self.wet_duration += self.dt * (self.h > wet_threshold)

        self.t += self.dt
        return self.dt

    # ------------------------------------------------------------------
    def mass_balance(self):
        """Return (stored volume, boundary inflow, clipping correction, % error)."""
        stored = float(self.h.sum()) * self.dx ** 2
        expected = self.volume_in + self.mass_correction
        if abs(expected) < 1e-9:
            pct = 0.0
        else:
            pct = 100.0 * (stored - expected) / expected
        return stored, self.volume_in, self.mass_correction, pct
