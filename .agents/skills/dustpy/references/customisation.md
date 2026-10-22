# DustPy Customisation Reference

Patterns below are taken from the official DustPy documentation/examples
(v1.0.9) and the source code. All code is CGS.

**DustPyLib first.** Several of these patterns have ready-made, well-tested
implementations in the companion library **DustPyLib** (`pip install dustpylib`)
— grid refinement (§2), planetary gaps (§6), planetesimal formation (§8), and
dust backreaction (§13). Check `references/dustpylib.md` and prefer the
DustPyLib version when the package is installed; the patterns below are the
DustPy-only fallbacks (and the general mechanisms DustPyLib plugs into).

Table of contents:
1. [Replacing a physics updater](#1-replacing-a-physics-updater)
2. [Custom radial grid / grid refinement](#2-custom-radial-grid--grid-refinement)
3. [Adding custom fields and groups](#3-adding-custom-fields-and-groups)
4. [Turning off physics modules](#4-turning-off-physics-modules)
5. [External sources and sinks](#5-external-sources-and-sinks)
6. [Planetary gaps (Kanagawa 2017 torque)](#6-planetary-gaps-kanagawa-2017-torque)
7. [Ice lines / spatially varying fragmentation velocity](#7-ice-lines--spatially-varying-fragmentation-velocity)
8. [Planetesimal formation sink](#8-planetesimal-formation-sink)
9. [Systoles and diastoles](#9-systoles-and-diastoles)
10. [Boundary conditions](#10-boundary-conditions)
11. [Changing the dust integrator](#11-changing-the-dust-integrator)
12. [Snapshot and output customisation](#12-snapshot-and-output-customisation)
13. [Dust backreaction (DustPyLib)](#13-dust-backreaction-dustpylib)

---

## 1. Replacing a physics updater

Every field has an `.updater` attribute — a simframe `Heartbeat` wrapping a
callable `f(sim) -> np.ndarray` (or `None`). Replacing it swaps the physics
for that field; setting it to `None` freezes the field at its current value.

```python
def my_temperature(sim):
    """Power-law midplane temperature profile."""
    T0 = 150.    # [K] at r0
    r0 = 1. * c.au
    q  = -0.5
    return T0 * (sim.grid.r / r0) ** q

sim = Simulation()
sim.gas.T.updater = my_temperature   # before initialize()
sim.initialize()
```

Rules:
- Signature is `f(sim)`; return an array matching the field shape.
- Updaters run **every time step** — keep them fast (vectorised NumPy).
- Two timing options:
  - **Before `initialize()`**: the custom updater is used from the start
    (default updaters are only assigned if the field is still `None`).
  - **After `initialize()`**: assign the value or the updater, then call
    `sim.update()` so all dependent fields are recomputed.
- If the field already has an updater and you assign a *new* updater
  directly (e.g. `sim.gas.T.updater = f`), systoles/diastoles attached to the
  old heartbeat are dropped — re-attach them if needed.

Example (from the docs, custom Keplerian-dependent quantity):

```python
sim = Simulation()
sim.gas.T.updater = my_temperature
sim.initialize()
sim.update()          # optional if you also changed values directly
```

---

## 2. Custom radial grid / grid refinement

Assign the cell *interfaces* `sim.grid.ri` (a plain `np.ndarray`, shape
`Nr+1`) before `sim.initialize()`; cell centres and annulus areas are derived
automatically (you can also call `sim.makegrids()` explicitly).

**Preferred (DustPyLib):** `refine_radial_local` refines the grid locally
around one or more locations (keeps the grid log-spaced):

```python
import numpy as np
from dustpy import Simulation
from dustpylib.grid.refinement import refine_radial_local

sim = Simulation()
ri = np.geomspace(sim.ini.grid.rmin, sim.ini.grid.rmax, sim.ini.grid.Nr)
ri = refine_radial_local(ri, 10. * c.au, num=3)    # refine around 10 AU
sim.grid.ri = ri
sim.initialize()
# sim.grid.Nr is now > 100
```

See `references/dustpylib.md` §1 for details (several refinement sites,
`num` levels, etc.).

Fallback (DustPy only, no DustPyLib) — hand-rolled local refinement:

```python
def refinegrid(ri, r0, num=3):
    """Double the resolution in `num` cells around r0, recursively."""
    if num == 0:
        return ri
    ind = np.argmin(r0 > ri) - 1
    indl, indr = ind - num, ind + num + 1
    ril, rir = ri[:indl], ri[indr:]
    N = (2 * num + 1) * 2
    rim = np.empty(N)
    for i in range(0, N, 2):
        j = ind - num + int(i / 2)
        rim[i]     = ri[j]
        rim[i + 1] = 0.5 * (ri[j] + ri[j + 1])
    return refinegrid(np.concatenate([ril, rim, rir]), r0, num=num - 1)

ri = refinegrid(ri, 10. * c.au, num=3)
sim.grid.ri = ri
sim.initialize()
```

Constraints:
- The **mass grid must stay strictly logarithmic** — never assign
  `sim.grid.m`; use `sim.ini.grid.mmin`, `mmax`, `Nmbpd` instead.
- If you change the mass grid, coagulation matrices must be recomputed,
  i.e. re-run `initialize()` on a fresh simulation.

---

## 3. Adding custom fields and groups

Use simframe's group/field machinery (pattern from the planetary-gaps
example):

```python
sim = Simulation()
sim.initialize()

sim.addgroup("planets")
sim.planets.addgroup("jupiter")
sim.planets.addfield("a", 5.2 * c.au, description="Semi-major axis [cm]")
sim.planets.addfield("M", 1. * c.M_jup, description="Mass [g]")

# Make a planet grow linearly in time (updater pattern)
def Mjup(sim):
    Mini, Mfin = 1. * c.M_earth, 1. * c.M_jup
    t0, t1 = 0., 1.e5 * c.year
    return np.clip((Mfin - Mini) / (t1 - t0) * (sim.t - t0) + Mini, Mini, Mfin)

sim.planets.jupiter.M.updater.updater = Mjup   # keep existing heartbeat
sim.planets.jupiter.updater = ["M"]
sim.planets.updater = ["jupiter"]
```

Note: assigning to `field.updater.updater` (rather than `field.updater =`)
preserves the heartbeat's systoles/diastoles.

---

## 4. Turning off physics modules

All patterns from the official docs (6_dust_evolution / 7_gas_evolution),
applied **after** `sim.initialize()`.

### Fragmentation off (sticking only)
```python
sim.dust.p.stick = 1.
sim.dust.p.frag  = 0.
sim.dust.p.updater = None
```
(Alternative: set `sim.dust.v.frag` to a huge value — but probabilities are
still evaluated every step, wasting CPU.)

### Coagulation off (hydrodynamics only)
```python
sim.dust.p.stick = 0.
sim.dust.p.frag  = 0.
sim.dust.p.updater     = None
sim.dust.v.rel.updater = None
sim.dust.v.frag.updater = None
sim.dust.kernel.updater = None
```

### Dust drift off
```python
sim.dust.v.rad = 0.
sim.dust.v.rad.updater = None
```

### Dust diffusion off
```python
sim.dust.D = 0.
sim.dust.D.updater = None
```

### Dust hydrodynamics off (drift + diffusion)
```python
sim.dust.v.rad = 0.
sim.dust.v.rad.updater = None
sim.dust.D = 0.
sim.dust.D.updater = None
```

### Gas hydrodynamics off (keep external sources)
Set the viscosity to zero and freeze it (only valid without dust
backreaction, i.e. `A = 1`, `B = 0`):
```python
sim.gas.nu[:] = 0.
sim.gas.nu.updater = None
```
Alternative: `sim.gas.alpha[:] = 0.; sim.gas.alpha.updater = None` — but
`alpha` feeds the dust mixing parameters too.

### Entire dust evolution off
Remove the dust integration instruction (index 0 of the default integrator):
```python
del sim.integrator.instructions[0]
```
Unset the updaters of dust fields you do not need, otherwise they are
recomputed every step (and dust backreaction can still affect the gas).

### Entire gas evolution off
```python
del sim.integrator.instructions[1]
```
The gas velocities are still computed unless you also freeze `nu`/`alpha` —
freeze them if dust still reads `gas.v.rad`.

---

## 5. External sources and sinks

`sim.dust.S.ext` (shape `(Nr, Nm)`) and `sim.gas.S.ext` (shape `(Nr,)`),
units g cm⁻² s⁻¹. Negative values are sinks. Dust boundary cells
(`[0]`, `[-1]`) are zeroed by the implicit preparator; gas boundary cells
are overwritten by the boundary conditions.

```python
def dust_infall(sim):
    """Constant infall rate inside 50 AU."""
    src = np.zeros_like(sim.dust.Sigma)
    mask = sim.grid.r < 50. * c.au
    src[mask, :] = 1.e-20 / c.year       # g cm^-2 s^-1 per mass bin (spread as you like)
    return src

sim.dust.S.ext.updater = dust_infall
```

Gas disk wind / photoevaporation:
```python
def disk_wind(sim):
    """Lose gas at a fixed rate per orbital time."""
    return -1.e-5 * sim.gas.Sigma * sim.grid.OmegaK   # [g/cm^2/s]

sim.gas.S.ext.updater = disk_wind
```

Set the field once and freeze it instead of an updater if the source is
constant:
```python
sim.gas.S.ext[:] = 0.
sim.gas.S.ext.updater = None
```

---

## 6. Planetary gaps (Kanagawa 2017 torque)

**DustPyLib first:** `dustpylib.substructures.gaps` provides `kanagawa2017`
and `duffell2020` gap-profile functions (drop-in for the `Kanagawa2017_gap_profile`
defined below), and the official DustPyLib gap example uses a different
mechanism — imposing the *inverse* gap profile on `sim.gas.alpha` (since
`ν·Σ = const` in steady state), which keeps the gas evolving. See
`references/dustpylib.md` §2 for the full pattern. Use that when DustPyLib
is installed; the torque-based pattern below is the DustPy-only approach.

The gap enters through `sim.gas.torque.Lambda` (specific angular momentum
injection rate, cm² s⁻²); the effective velocity `v = 2Λ/(Ω_K r)` is
already wired into the gas Jacobian (`std.gas.vtorque`). Pattern from the
official example (multiple planets, grid refinement, growing planets):

```python
import numpy as np
from dustpy import Simulation
import dustpy.constants as c

def Kanagawa2017_gap_profile(r, a, q, h, alpha):
    """Gap perturbation f = Sigma_g(r) / Sigma_0 (Kanagawa et al. 2017)."""
    ret = np.ones_like(r)
    dist = np.abs(r - a) / a
    K   = q**2 / (h**5 * alpha)
    Kp  = q**2 / (h**3 * alpha)
    Kp4 = Kp**0.25
    SigMin  = 1. / (1 + 0.04 * K)
    SigGap  = 4 / Kp4 * dist - 0.32
    dr1 = (0.25 * SigMin + 0.08) * Kp4
    dr2 = 0.33 * Kp4
    mask1 = np.logical_and(dr1 < dist, dist < dr2)
    ret = np.where(mask1, SigGap, ret)
    ret = np.where(dist < dr1, SigMin, ret)
    return ret

sim = Simulation()
# refine the grid around the planets (see §2)
ri = np.geomspace(1., 1.e3, 100) * c.au
ri = refinegrid(ri, 5.2 * c.au)
ri = refinegrid(ri, 9.6 * c.au)
sim.grid.ri = ri
sim.initialize()

# planets as custom fields (see §3)
sim.addgroup("planets")
sim.planets.addgroup("jupiter")
sim.planets.jupiter.addfield("a", 5.2 * c.au, description="Semi-major axis [cm]")
sim.planets.jupiter.addfield("M", 1. * c.M_jup, description="Mass [g]")
sim.planets.updater = ["jupiter"]

from scipy.interpolate import interp1d

def Lambda_torque(sim):
    f = np.ones(sim.grid.Nr)
    fi = np.ones(sim.grid.Nr + 1)
    for name, planet in sim.planets:
        q = planet.M / sim.star.M
        _h = interp1d(sim.grid.r, sim.gas.Hp / sim.grid.r)
        h = _h(planet.a)
        _alpha = interp1d(sim.grid.r, sim.gas.alpha)
        alpha = _alpha(planet.a)
        f  *= Kanagawa2017_gap_profile(sim.grid.r,  planet.a, q, h, alpha)
        fi *= Kanagawa2017_gap_profile(sim.grid.ri, planet.a, q, h, alpha)
    grad = np.diff(np.log10(fi)) / np.diff(np.log10(sim.grid.ri))
    return -1.5 * sim.grid.OmegaK * sim.gas.nu * (1. / (2. * f) - grad - 0.5)

sim.gas.torque.Lambda[:] = Lambda_torque(sim)
```

Notes:
- Compute the gap-gradient on the **cell interfaces** (`ri`), as above.
- If planets grow, give their mass fields a time-dependent updater
  (`sim.planets.jupiter.M.updater.updater = Mjup`, see §3).
- A growing planet needs a growing torque: either recompute `Lambda` in a
  systole/diastole, or make `Lambda_torque` depend on the live planet mass
  (it does, via `planet.M`).

---

## 7. Ice lines / spatially varying fragmentation velocity

Make `sim.dust.v.frag` depend on the local gas temperature (Pinilla et al.
2017 style; pattern from the official example):

```python
sim = Simulation()
sim.initialize()

def vfrag(sim):
    vfrag = np.where(sim.gas.T < 150., 1000., 100.)   # water ice line at 150 K
    vfrag = np.where(sim.gas.T <  80.,  700., vfrag)  # NH3 at 80 K
    vfrag = np.where(sim.gas.T <  44.,  100., vfrag)  # CO2 at 44 K
    return vfrag                                      # [cm/s]

sim.dust.v.frag = vfrag(sim)   # static profile: set once after initialize
sim.update()                   # recompute dependent fields
```

If your temperature profile **evolves** (custom `gas.T.updater`), assign the
updater instead so it is recomputed each step:
```python
sim.dust.v.frag.updater = vfrag
```

---

## 8. Planetesimal formation sink

**DustPyLib first:** `dustpylib.planetesimals.formation` provides the
published prescriptions as ready-made `S.ext` source terms —
`drazkowska2016` (pebble threshold, per-orbit), `schoonenberg2018` (all-
particle threshold, per settling time) and `miller2021` (smooth tanh
transition) — plus the pattern for tracking the converted planetesimal mass
with a differentiator + explicit-Euler instruction. Use them when DustPyLib
is installed; see `references/dustpylib.md` §3 (parameters, citations).

Minimal fallback (DustPy only): convert dust into planetesimals where a
threshold (e.g. streaming-instability criterion `eps > eps_crit`) is
exceeded, using a negative external source:

```python
eps_crit = 1.0

def planetesimal_sink(sim):
    """Remove dust where eps > eps_crit, on one orbital timescale."""
    sink = np.zeros_like(sim.dust.Sigma)
    eps  = sim.dust.eps                       # (Nr,)
    over = eps > eps_crit
    if np.any(over):
        tau = 1. / sim.grid.OmegaK            # one orbital period [s]
        sink[over, :] = -sim.dust.Sigma[over, :] / tau[over, np.newaxis]
    return sink

sim.dust.S.ext.updater = planetesimal_sink
```

Track the converted mass over time with a diastole on a custom field (see §9
and §3), e.g. accumulate `∫ sink · A dt` per snapshot.

---

## 9. Systoles and diastoles

Every updater is a simframe `Heartbeat` with three parts:
- **systole** — called *before* the updater,
- **updater** — the update function,
- **diastole** — called *after* the updater.

Both take `sim` and return nothing (they act on the frame). Assign functions
directly to the heartbeat attributes:

```python
def update_snowline(sim):
    """Runs after gas temperature is updated (diastole of gas.T)."""
    sim.grid.rsnow.update()    # e.g. a custom field added via addgroup/addfield

sim.gas.T.updater.diastole = update_snowline
```

Typical uses:
- per-timestep state updates that must respect field update order
  (diastole of the field they depend on, systole if they must run first);
- bookkeeping per update (masses, planet growth, torque re-evaluation).

Caveats:
- Assigning a *new updater* to a field (`sim.x.updater = f`) replaces the
  heartbeat and drops its systoles/diastoles. To swap only the function,
  assign to `sim.x.updater.updater = f`.
- Assigning lists to group updaters always resets their systoles/diastoles
  to `None`.

---

## 10. Boundary conditions

Boundaries are `dustpy.utils.Boundary` objects operating on the 3 outermost
cells. Defaults (set in `initialize()`):

| | Gas | Dust |
|---|---|---|
| Inner | constant gradient | constant gradient |
| Outer | fixed value = `SigmaFloor` | fixed value = 0.1 × `SigmaFloor` |

Available conditions (`setcondition`): `None` (no BC), `"const_grad"`,
`"const_pow"`, `"const_val"`, `"val"`, `"grad"`, `"pow"` (the last three
take a value).

```python
# After initialize():
sim.gas.boundary.inner.setcondition("const_val")     # closed inner wall
sim.dust.boundary.outer.setcondition("val", 1.e-6)   # fixed value
```

For implicit dust integration the Jacobian already honours the boundary
condition type (see `std.dust.jacobian`), so changing the condition is
enough — no Jacobian changes required.

---

## 11. Changing the dust integrator

The default is an implicit 1st-order direct-solver (sparse LU) for
`dust.Sigma`. Switch to the explicit 5th-order adaptive Cash-Karp scheme
when the implicit matrix is weakly conditioned:

```python
sim.setdustintegrator(scheme="explicit", method="cash-karp")
```

The function also swaps the dust update order to the explicit one
(`... Fi ...` and `S.updater = ["coag", "hyd", "ext", "tot"]`) and the
preparator/finalizer. Call it **after** `initialize()` (it needs the grids).
Explicit options: `scheme="explicit", method="cash-karp"`. Implicit options:
`scheme="implicit", method="direct"` (gmres/bicgstab not implemented).

---

## 12. Snapshot and output customisation

```python
import numpy as np
import dustpy.constants as c

sim = Simulation()
sim.initialize()

# 5 snapshots per decade from 1 kyr to 1 Myr
sim.t.snapshots = np.hstack([
    sim.t,
    np.geomspace(1e3, 1e6, num=16) * c.year,   # [s]
])

sim.writer.datadir   = "runs/my_run"
sim.writer.overwrite = False                    # safe default

sim.run()
```

Reduce file size by dropping fields from the output (they are still
computed):
```python
sim.dust.kernel.save = False
sim.dust.v.rel.azi.save = False
```

Extend an existing run by appending snapshot times (see `output.md`):
```python
from dustpy import readdump
sim = readdump("data/frame.dmp")
sim.t.snapshots = np.hstack((sim.t.snapshots, [100001. * c.year]))
sim.run()
```

---

## 13. Dust backreaction (DustPyLib)

Dust backreaction (the dynamical feedback of drifting dust on the gas — pull
and push, coefficients A, B) is **off by default** in DustPy (`A = 1`,
`B = 0`). DustPyLib wires it in with one call (Garaté et al. 2019, 2020):

```python
from dustpylib.dynamics.backreaction import setup_backreaction

sim = Simulation()
sim.initialize()
setup_backreaction(sim)                # standard (vertically uniform d2g)
# or: setup_backreaction(sim, vertical_setup=True)  # accounts for dust settling
sim.run()
```

What it changes on the frame (verified in source):
- `sim.dust.backreaction.updater` → computes `A`, `B` per radius from the
  dust-to-gas and Stokes-number distribution (summed over mass bins); the
  outermost cell is pinned to (1, 0) for stability.
- `sim.dust.D.updater` → damps the diffusivity by `1/(1 + d2g_tot)`
  (important at high dust-to-gas ratios).
- With `vertical_setup=True`: adds `backreaction.A_dust_settling` /
  `B_dust_settling` (shape `Nr, Nm`) and replaces `sim.dust.v.rad.updater`
  with a per-species version.

The coefficients are consumed by DustPy's existing `gas.v.rad`
(`v_gas = A·v_visc + 2B·v_η`) and `dust.v.driftmax` — no further changes
needed. Cite Garaté et al. (2019), ApJ 871, 53, and Garaté et al. (2020),
A&A 635, A149 (the function prints the reminder).

Caveat: once backreaction is active (`A ≠ 1`, `B ≠ 0`), the "gas
hydrodynamics off" pattern in §4 (zeroing `nu`) is inconsistent — don't
combine them.

Full API and parameters: `references/dustpylib.md` §4.
