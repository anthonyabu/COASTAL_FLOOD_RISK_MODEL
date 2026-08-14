# Coastal flood risk assessment for transportation infrastructure

A 2D hydrodynamic inundation model built from the governing equations, coupled to
asset-level depth-damage functions and integrated to expected annual damage.

Author: Anthony Ariu Abu
Status: self-directed exercise, synthetic domain, unpublished

---

## 1. What this is and what it is not

This is a working demonstration of the full chain from governing equations to a risk
number: storm surge forcing → 2D shallow water hydrodynamics → peak inundation depth
over specific assets → depth-damage function → expected annual damage → asset ranking.

It is **not** a study of a real place. The terrain, road network, unit costs, surge
distribution and depth-damage curve are all synthetic or illustrative. They are
labelled as such in the code and in every output file. Nothing here should be quoted
as a result about anywhere.

What it does demonstrate is the method, and three findings that fall out of it which
are general rather than site-specific (Section 6).

---

## 2. Hydrodynamics

### Governing equations

The **local-inertial** (inertial-diffusive) approximation of the 2D shallow water
equations. The advective acceleration term is dropped from the momentum equation,
which is standard for floodplain and coastal overland flow where advection is small
relative to gravity, friction and local acceleration. This is the scheme underlying
LISFLOOD-FP and several operational flood models.

Momentum, per unit width, semi-implicit in friction:

```
                q^n  -  g · h_flow · Δt · ∂(h+z)/∂x
q^(n+1)  =  ─────────────────────────────────────────────
             1  +  g · Δt · n² · |q^n| / h_flow^(7/3)
```

Continuity:

```
∂h/∂t = ( Σ q_in − Σ q_out ) / Δx
```

The conveyance depth at a face is

```
h_flow = max(h_i + z_i, h_j + z_j) − max(z_i, z_j)
```

which is what allows the scheme to handle wet/dry fronts and flow over embankments
without special-casing.

### Discretisation and stability

- Staggered (Arakawa C) grid: depths and bed elevations at cell centres, unit-width
  discharges at cell faces. Avoids the chequerboard oscillation that co-located
  variables produce.
- Adaptive timestep from the gravity-wave CFL condition,
  `Δt = α · Δx / √(g · h_max)` with α = 0.5.
- **Volume-based flux limiter.** A cell cannot export more water than it holds. Each
  cell's total outflow is computed before the continuity update; if it exceeds the
  stored volume, every outgoing face flux from that cell is scaled by the same factor.
  This was added after an initial implementation used post-hoc depth clipping instead
  and lost 30% mass conservation on a dry-bed spreading test.

### Domain

100 × 60 cells at 20 m (2.0 km cross-shore × 1.2 km alongshore). Concave-up coastal
plain from −1.2 m MSL at the shore to +5.3 m at 2 km, with a tidal creek providing
conveyance inland. Manning's n: 0.02 road, 0.035 creek, 0.055 plain, 0.10 marsh.

### Forcing

Semi-diurnal tide (1.40 m amplitude) plus a raised-cosine surge, with astronomical
high water aligned to the surge peak. This is the conservative joint-occurrence case,
not the expected case — a full analysis would treat tide/surge timing as a joint
probability rather than assuming coincidence.

Surge magnitude from a GEV inverse CDF **fitted to NOAA CO-OPS station 8418150
(Portland, Maine), 1980-2025, n = 46 annual maxima** — see Section 3a. Twelve return
periods from 1.5 to 500 years.

---

## 3. Verification

The solver was tested against cases with known answers before being used.

| Test | Expected | Result |
|---|---|---|
| T1 Closed flat basin, initial mound | Level free surface, exact mass conservation | Volume error 0.0000%, water surface flat to 2.6 mm, mean depth exactly 0.125 m as predicted |
| T2 Closed tilted basin | Level free surface, exact conservation | Volume error 0.0000% |
| T3 Sloping shore, fixed sea level | Fills to sea level, shoreline where bed = SWL | Converged to 0.6001 m against a 0.600 m target; shoreline within one cell of analytic; zero clipping correction |
| T4 Dry-bed dam break | Bounded, monotone, conservative | Volume error 0.0000%, depths stay within [0, h₀] |

T4 front position underpredicts the Ritter analytical solution (570 m vs 688 m at
t = 30 s). This is expected and not a defect: dropping advective acceleration is
known to retard the dry-bed front. It bounds where the scheme is appropriate — it is
valid for floodplain flow, not for supercritical or strongly advective flow.

Mass balance was audited on every production run. Worst error across all 12 return
periods: **0.00000%**, with zero flux-limiter correction.

---

## 3a. Extreme value analysis of the surge forcing

The surge distribution was originally invented. It is now fitted to observed water
levels (`code/fit_surge_gev_v2.R`).

**Method.** Hourly observed water level minus hourly astronomical prediction gives the
residual. NOAA predictions reference the 1983-2001 tidal epoch, so the residual carries
accumulated sea level rise as well as meteorology; the annual mean residual is therefore
removed year by year, leaving meteorological surge. Annual maxima are extracted with an
80% data-coverage filter, and a GEV is fitted by maximum likelihood in the Coles
parameterisation.

**Validation.** Two independent checks:

- The slope of the *removed* offset is a local sea level rise estimate: **+2.83 mm/yr
  (p < 0.0001)**. This was not fitted to anything — it falls out of the datum handling —
  and can be compared against NOAA's published trend for the station.
- The hand-rolled MLE was cross-checked against the `extRemes` package, which returned
  **+0.6900, +0.1353, -0.1319 — identical to four decimal places**.

**Result.**

| | Illustrative (v1) | Fitted (v2) |
|---|---|---|
| Location | 1.05 m | 0.6900 m |
| Scale | 0.38 m | 0.1353 m |
| Shape | 0 (Gumbel assumed) | −0.1319 (bounded tail) |
| 100-yr surge | 2.80 m | **1.16 m** |

The invented distribution overstated the 100-year surge by a factor of 2.4, almost
entirely through a scale parameter nearly three times too large.

**Stationarity.** No trend in annual maximum surge (+1.09 mm/yr, p = 0.53), so a
stationary GEV is defensible over this window. Sea level rise was removed before fitting
and must be added separately to still water level — it is not in this distribution.

**The shape parameter matters.** ξ = −0.13 implies a bounded upper tail. Assuming
Gumbel with the same location and scale would overstate the 100-year surge by 13.5% and
the 500-year by 21.1%.

**Domain rescaling.** Road crest elevations in the synthetic terrain were originally
chosen against the invented 4.2 m still water level. Under the fitted distribution the
100-year level is 2.56 m MSL, so the crests were rescaled (highway 3.5 → 2.35 m,
causeway 2.4 → 2.10 m, inland arterial 3.8 → 2.55 m) along with the ground profile.
Without this the assets would never flood at any return period.

---

## 4. Damage assessment

### Asset representation

Three road assets, split into 100 m segments (29 segments, $6.88 M total exposure):

| Asset | Length | Crest | Exposure |
|---|---|---|---|
| A Coastal Highway (shore-parallel embankment, bridged creek opening) | 1200 m | 2.35 m | $2.88 M |
| B Access Causeway (shore-perpendicular, low) | 500 m | 2.10 m | $1.60 M |
| C Inland Arterial | 1200 m | 2.55 m | $2.40 M |

Segment depth is the **maximum** over the segment's cells, not the mean. For a linear
asset, functional loss is governed by the worst point — a single washed-out 20 m
length closes the whole segment. Using the mean would systematically understate damage.

### Two elevation fields

The model maintains bed elevation `z` and asset surface elevation `deck` separately.
At a bridge these differ by several metres: `z` is the channel bed, because that is
what governs conveyance through the opening, while `deck` is the carriageway above it.
Damage is assessed against `(peak water surface elevation − deck)`, never against
water depth over the bed. See Section 6.2 for why this matters.

### Depth-damage function

Two illustrative curves are retained deliberately, because the choice between them
turns out to be a first-order driver of the result:

- `no_threshold` — damage begins as soon as water reaches the carriageway
- `threshold` — damage negligible until water is deep enough to affect the pavement
  structure and subgrade rather than merely covering the surface (default)

**Neither is sourced.** Published road depth-damage functions exist and one must be
substituted before any result is used. A duration multiplier is implemented but
disabled by default, because the physical mechanism (subgrade saturation) is real but
the multiplier values would be invented.

### Expected annual damage

EAD is the integral of damage over annual exceedance probability, by trapezoidal rule.
Two truncations are made explicit rather than hidden: damage is held constant above the
rarest event sampled, and tapered linearly to zero below the most frequent. Both terms
are reported separately in the output so their contribution is visible.

---

## 5. Results

**Expected annual damage: $13,210 /yr against $6.88 M exposure (0.19 %/yr),** using
gauge-fitted surge statistics.

| Asset | EAD | % of asset value per year |
|---|---|---|
| B Access Causeway | $12,388 | 0.77 % |
| A Coastal Highway | $822 | 0.03 % |
| C Inland Arterial | $0 | 0.00 % |

The prioritisation result is that the **causeway carries 94 % of the risk while
representing 23 % of the exposure**. It is the cheapest asset and the one to act on
first. A ranking by asset value, or by peak-event damage alone, would have pointed at
the Coastal Highway — the highway takes more damage in a rare event, but the causeway
is damaged in events that recur every few years, and expected annual damage is
dominated by frequency rather than severity.

Substituting the fitted distribution for the invented one reduced EAD by a factor of
15. The asset *ranking* was unchanged. This is the useful pattern: relative
prioritisation is far more robust to hazard misspecification than absolute loss is,
which matters because prioritisation is usually the decision actually being made.

---

## 6. Three findings worth stating

### 6.1 EAD is far more sensitive to the damage curve than to the hazard

Holding the hydraulics completely fixed and changing only the depth-damage curve
changed EAD by a factor of **7.6**, while 100-year event damage changed by 3.9.

This effect grew sharply once real surge statistics were used. Under the invented
distribution the factor was 1.8; under the fitted one it is 7.6. The reason is that the
true surge range is narrow, so inundation depths over the assets sit mostly in the
0–0.3 m band — precisely where the two candidate damage curves disagree most.

The conclusion is that the vulnerability model, not the hazard model, is the binding
constraint on this estimate. Effort spent refining the hydraulics past a certain point
is misallocated while the damage curve remains uncalibrated.

### 6.2 Conflating hydraulic and structural elevation inflates risk by ~3×

The first working version assessed damage from water depth over the bed. At the bridged
creek crossing this read several metres of water in the channel as several metres of
water on the carriageway, so the highway appeared severely damaged in events whose
surge never approached its crest. Separating `z` from `deck` moved EAD from $539,930/yr
to $196,827/yr — a factor of 2.7 from one modelling distinction, with no change to the
hydraulics.

This is a general trap for network-scale assessments where assets are rasterised onto a
DEM: bridges, culverts and elevated structures are exactly the assets where the two
elevations diverge most, and exactly the assets that matter most for network
connectivity.

### 6.4 In a macrotidal, modest-surge setting, tide timing dominates storm size

Across the entire 1.5-to-500-year range the peak still water level spans only **0.59 m**
(2.08 to 2.66 m MSL), because the tidal amplitude (1.40 m) exceeds the surge at every
return period. Whether an asset floods is therefore governed largely by whether a storm
coincides with high water, not by how large the storm is.

This model assumes coincidence, which is the conservative design case. A defensible
assessment for this coastline would have to treat tide-surge timing as a joint
probability rather than assuming the worst alignment, and that choice would plausibly
matter more to the answer than any refinement of the hydraulic scheme.

### 6.3 Static bathtub maps are adequate here — but only because the basin equilibrates fast

Holding peak water level at the 100-year value and varying only surge duration:

| Surge duration | Damage vs static bathtub |
|---|---|
| 1 h | 74 % |
| 2 h and longer | 100–108 % |

And roughness only matters while the basin is still filling:

| Marsh Manning's n | 0.5 h surge | 1 h surge | 3 h surge |
|---|---|---|---|
| 0.03 | $1.59 M | $1.72 M | $2.49 M |
| 0.30 | $0.66 M | $0.96 M | $2.01 M |
| **Reduction** | **58 %** | **44 %** | **19 %** |

Friction controls the *rate* of filling, not the equilibrium level. For a 2 km
low-relief plain and a multi-hour surge, the static map is a good approximation and the
dynamic model buys little. For short, fast-moving events it overstates damage by ~26 %,
and it cannot represent marsh attenuation at all, since roughness appears nowhere in a
bathtub calculation.

The practical implication: a screening study using static maps should state the surge
durations over which that assumption holds. Any assessment of marsh or living
shorelines as natural coastal protection requires a dynamic model by construction.

---

## 7. Known limitations

- Terrain, road network, unit costs are synthetic; road crests were rescaled to sit in
  a plausible relationship to the fitted still water levels
- Surge distribution is fitted to a real gauge, but tide amplitude (1.40 m) is assumed
  and should be checked against the published NOAA datums for the station
- Depth-damage curve shape is plausible, ordinates are not sourced — and per Section 6.1
  this is now the dominant source of uncertainty in the result
- Sea level rise is removed from the surge distribution and is NOT added back into still
  water level anywhere in this model
- Tide and surge peaks assumed coincident (conservative, not expected)
- Wave setup, wave overtopping, rainfall and groundwater are not modelled
- Local-inertial scheme omits advective acceleration — invalid for supercritical flow
- No traffic disruption or network-detour cost; direct asset damage only, which for
  transportation infrastructure is usually the smaller share of total economic impact
- Single deterministic run per return period; no uncertainty propagation through the
  damage function

---

## 8. Files

| File | Contents |
|---|---|
| `solver.py` | Local-inertial 2D shallow water solver with flux limiter |
| `test_solver.py` | Four verification cases with known answers |
| `terrain.py` | Synthetic domain, dual elevation fields, asset inventory |
| `damage.py` | Surge hazard, depth-damage curves, EAD integration |
| `run_event.py` | Single-event driver |
| `analysis.py` | Ensemble across return periods, damage curve sensitivity, figures |
| `duration_study.py` | Static vs dynamic, duration and roughness sensitivity |
| `fit_surge_gev_v2.R` | Extreme value analysis of NOAA CO-OPS gauge records (R) |

Reproduce with `python3 analysis.py && python3 duration_study.py`.
Runtime approximately 5 minutes for the 12-event ensemble on a single core.
