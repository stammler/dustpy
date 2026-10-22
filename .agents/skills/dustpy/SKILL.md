---
name: DustPy
description: >
  Set up, run, customise, and analyse DustPy (dustpy.Simulation) simulations of dust
  evolution in protoplanetary disks: gas and dust radial evolution, coagulation and
  fragmentation (Smoluchowski equation), radial drift, diffusion, ice lines, planetary
  gaps, and planetesimal formation. Also covers the companion library DustPyLib
  (dustpylib): backreaction, grid refinement, planetesimal formation, planetary gaps,
  RADMC-3D and slab radiative transfer. Use whenever a task involves DustPy code,
  writing or modifying a DustPy setup script, dustpy.std updaters, reading DustPy HDF5
  output (data/dataNNNN.hdf5), dustpy.plot, dustpy.utils.read_data, resuming a run from
  frame.dmp, or analysing dust-to-gas ratio, fragmentation barriers, or grain growth.
  Load this skill before writing or modifying any DustPy-related code.
---

# DustPy

DustPy (Stammler & Birnstiel 2022, ApJ 935, 35) simulates the coupled radial
evolution of gas and dust in protoplanetary disks: viscous gas evolution,
dust advection and diffusion, and grain growth via the Smoluchowski equation.
It is built on the `simframe` framework: a `Simulation` frame of `Field`
objects whose physics is set by `updater` functions, an `Integrator` that
steps the frame, and a `writer` that stores snapshots.

**All internal quantities are in CGS units** (cm, g, s, K). Unit confusion is
the single most common source of errors. Convert with `dustpy.constants as c`
(`c.au`, `c.M_sun`, `c.M_earth`, `c.M_jup`, `c.R_sun`, `c.m_p`, `c.year`, ...).

---

## Quick reference

| | |
|---|---|
| Install | `pip install dustpy` (Python 3 + Fortran compiler; CMake-free meson build) |
| Main class | `from dustpy import Simulation` |
| Constants | `from dustpy import constants as c` |
| Companion lib | `pip install dustpylib` — backreaction, grid refinement, planetesimal formation, gaps, RADMC-3D/slab RT; see `references/dustpylib.md` |
| Default model | 100 radial cells (1–1000 AU), 7 mass bins/decade (10⁻¹²–10⁵ g), α=10⁻³, M_disk=0.05 M☉, d2g=0.01, v_frag=100 cm/s; 21 snapshots over 1 kyr–100 kyr |
| Output | `<datadir>/data{N:04d}.hdf5` (snapshots) + `<datadir>/frame.dmp` (restart dump) |
| Analysis | `dustpy.utils.read_data(datadir)` (lightweight), `dustpy.hdf5writer` (full), `h5py` (direct) |
| Resume | `from dustpy import readdump; sim = readdump("data/frame.dmp")` |
| Plotting | `dustpy.plot.panel("data")` (interactive) or the headless scripts in `scripts/` |
| Citation | Always remind the user: Stammler & Birnstiel (2022), doi:10.3847/1538-4357/ac7d58 — plus the per-module reference(s) when DustPyLib is used (see `references/dustpylib.md`) |

---

## Workflow

Classify the task first, then follow the matching sequence:

**A. New setup and run**
1. Read the user's physical requirements; convert to CGS.
2. Write a setup script using the standard pattern below (or run
   `scripts/run_dustpy.py` for the default model with standard parameters).
3. Validate parameters (see iron rules: `Nmbpd ≥ 7`, `mmax` above the
   fragmentation barrier, `rmin` > 0, `t_end > t_start`).
4. Run, then sanity-check the output (see "Sanity checks").
5. Plot with `scripts/plot_dustpy.py` and report results.

**B. Analyse existing output**
1. Verify the output directory exists and list `data*.hdf5` files first.
2. Read a **single** snapshot before any batch processing; print shapes and
   ranges to confirm what you are looking at.
3. Sanity-check (below); report `[SANITY WARN]` items instead of ignoring them.
4. State the parameter provenance: read grid/star/gas parameters from the
   HDF5 file (`grid/r`, `grid/m`, `t`, ...) — never assume defaults.

**C. Customise physics**
1. Check `references/dustpylib.md` first — DustPyLib provides ready-made,
   well-tested versions of several patterns (backreaction, grid refinement,
   planetesimal formation, planetary gaps, radiative transfer). Prefer them
   over hand-rolled code when DustPyLib is installed.
2. Read `references/customisation.md` for the relevant pattern (updater
   replacement, custom grid, turning physics off, gaps, ice lines, sources,
   systoles/diastoles, boundaries, integrator).
3. Apply it with the exact timing rules in that file (before vs after
   `initialize()`).
4. After any updater change, call `sim.update()` before `sim.run()` so
   dependent fields are consistent.

**D. Resume / extend a run**
1. Use `readdump("<datadir>/frame.dmp")` (restores the full frame including
   updaters), append new snapshot times, then `sim.run()`. See `references/output.md`.

---

## Standard setup pattern

```python
from dustpy import Simulation
from dustpy import constants as c
import numpy as np

sim = Simulation()

# --- stellar parameters (CGS) ---
sim.ini.star.M = 1.0 * c.M_sun      # stellar mass [g]
sim.ini.star.R = 2.0 * c.R_sun      # stellar radius [cm]
sim.ini.star.T = 5772.              # effective temperature [K]

# --- grid ---
sim.ini.grid.Nr    = 100            # radial cells
sim.ini.grid.rmin  = 1.0  * c.au    # inner boundary [cm] (> 0)
sim.ini.grid.rmax  = 1000. * c.au   # outer boundary [cm]
sim.ini.grid.Nmbpd = 7              # mass bins per decade (>= 7 required)
sim.ini.grid.mmin  = 1.e-12         # [g]
sim.ini.grid.mmax  = 1.e5           # [g]; must exceed fragmentation barrier

# --- gas ---
sim.ini.gas.alpha    = 1.e-3        # turbulent alpha
sim.ini.gas.Mdisk    = 0.05 * c.M_sun  # [g]
sim.ini.gas.SigmaExp = -1.0         # power-law index of initial profile
sim.ini.gas.SigmaRc  = 60. * c.au   # cut-off radius [cm]

# --- dust ---
sim.ini.dust.d2gRatio   = 0.01      # initial dust-to-gas ratio
sim.ini.dust.aIniMax    = 1.e-4     # max initial grain size [cm]
sim.ini.dust.vFrag      = 100.      # fragmentation velocity [cm/s]
sim.ini.dust.rhoMonomer = 1.67      # monomer bulk density [g/cm^3]

sim.initialize()

# --- writer and snapshots (set AFTER initialize()) ---
sim.writer.datadir  = "data"
sim.writer.overwrite = False        # True only after explicit user confirmation
sim.t.snapshots = np.hstack([
    sim.t,
    np.geomspace(1e3, 1e5, num=21) * c.year   # [s], 10 per decade over 1 kyr-100 kyr
])

sim.run()
```

All other `sim.ini` parameters (with types, defaults, constraints) are
documented in `references/parameters.md`.

---

## Simulation object structure

```
sim.t          time [s]; sim.t.snapshots (output times), sim.t.suggest(dt)
sim.star       M [g], R [cm], T [K], L [erg/s]
sim.grid       Nr, r (cm), ri (cm, cell interfaces), A (cm^2),
               Nm, m (g, strictly log grid), OmegaK (1/s)
sim.gas        Sigma, SigmaFloor, T, cs, Hp, nu, mu, n, mfp, P, rho,
               alpha, eta, Fi (interface fluxes), v.visc, v.rad,
               S.ext, S.hyd, S.tot, torque.Lambda, torque.v
sim.dust       Sigma (Nr,Nm), SigmaFloor, a, St, H, rho, rhos, D, eps (Nr,),
               kernel (Nr,Nm,Nm), v.rad, v.frag (Nr,), v.driftmax,
               v.rel.{azi,brown,rad,turb,vert,tot}, delta.{rad,turb,vert},
               Fi.{adv,diff,tot}, p.{frag,stick}, S.{coag,hyd,ext,tot},
               backreaction.{A,B}, coagulation.{stick,stick_ind,A,eps,lf_ind,rm_ind,phi}
sim.integrator simframe.Integrator (instructions, preparator, finalizer)
sim.writer     simframe hdf5writer (datadir, overwrite, filename, extension)
sim.ini        parameter namespace: .star, .grid, .gas, .dust
```

- 1-D fields have shape `(Nr,)` (e.g. `gas.Sigma`, `dust.eps`).
- 2-D fields have shape `(Nr, Nm)` (e.g. `dust.Sigma`, `dust.St`).
- Collision-matrix fields have shape `(Nr, Nm, Nm)` (e.g. `dust.kernel`, `dust.p.frag`).
- Every field has an `.updater` (a callable `f(sim)` or `None`).
  Setting `.updater = None` freezes the field at its current value.
- Fields can be excluded from output with `sim.dust.kernel.save = False`.

---

## Reading output

HDF5 layout (one file per snapshot, CGS):

```
/t            scalar [s]
/grid/r       (Nr,) cm          /grid/m   (Nm,) g
/gas/Sigma    (Nr,) g/cm^2      /gas/T    (Nr,) K
/dust/Sigma   (Nr,Nm) g/cm^2    /dust/St  (Nr,Nm)
/dust/a       (Nr,Nm) cm        /dust/eps (Nr,) -
/dust/rho     (Nr,Nm) g/cm^3    (midplane volumetric, per mass bin)
```

```python
import h5py
import dustpy.constants as c

with h5py.File("data/data0000.hdf5", "r") as f:
    r     = f["grid/r"][:]        # (Nr,) cm
    m     = f["grid/m"][:]        # (Nm,) g
    t     = float(f["t"][()])     # s
    Sig_g = f["gas/Sigma"][:]     # (Nr,) g/cm^2
    Sig_d = f["dust/Sigma"][:]    # (Nr, Nm) g/cm^2
print(r.shape, m.shape, t / c.year)   # verify before batch processing
```

Prefer the helpers when available:

```python
from dustpy.utils import read_data          # lightweight: only plot-relevant fields,
data = read_data("data")                    # with derived quantities (masses, St limits)
# data.dust.sigma: bin-normalised dust density (Nt, Nr, Nm)
# data.dust.St_limits.frag / .drift: (Nt, Nr) growth-limit Stokes numbers

from dustpy import hdf5writer               # full reader, simframe-based
wrtr = hdf5writer(datadir="data")
one  = wrtr.read.output(0)                  # one snapshot, attribute access
seq  = wrtr.read.sequence("gas.Sigma")      # one field over all times (Nt, Nr)
```

`read_data` returns arrays with a leading time axis `(Nt, ...)`. Note
`read_data` exposes `dust.Sigma` (surface density per mass bin) and
`dust.sigma` (bin-normalised, ≈ Sigma / (Δm/m)); use the one you need
explicitly.

---

## Plotting

Headless (agent) plotting — no display required, wraps the same logic as
`dustpy.plot.panel`. Run from the skill directory (paths below are relative to
this SKILL.md):

```bash
# All three plot types for the last/first snapshots of a run directory:
python scripts/plot_dustpy.py --run_dir data/my_run --plot all --out plots/my_run

# Radial profiles at snapshot 0 and the last one only:
python scripts/plot_dustpy.py --run_dir data/my_run --plot radial --snaps 0 -1

# 6-panel overview at snapshot 5, PNG output:
python scripts/plot_dustpy.py --run_dir data/my_run --plot panel --snap 5 --fmt png
```

`--run_dir` accepts either a run directory containing a `data/` subdirectory or
a data directory with `data*.hdf5` files directly. Last stdout line on success:
`SUCCESS plot_dir=<path> files=<comma-separated paths>`; on failure:
`ERROR message=<description>`.

| Plot | Contents |
|---|---|
| `panel` | 6-panel overview: dust density map (r–m) with St=1/drift/frag limits, mass spectrum, gas+dust mass evolution, radial profile, Σ(r) + ε(r) |
| `radial` | Σ_gas/Σ_dust, ε(r), a_max(r), St_max(r) at chosen snapshots |
| `evolution` | space-time diagrams a_max(r,t) and Σ_dust(r,t) |

For interactive use in a notebook: `from dustpy import plot; plot.panel("data", it=10, ir=10, im=5)` or `plot.ipanel("data")`.

---

## Sanity checks (mandatory before reporting results)

Run on every snapshot you report on. Emit `[SANITY WARN]` lines for
out-of-range values; do not silently continue, and do not treat warnings as
failures without checking the setup.

```python
def sanity_check(snap, snap_index):
    """snap: dict from one snapshot with keys r, m, Sigma_gas, T_gas,
    Sigma_dust (Nr,Nm), St (Nr,Nm), a (Nr,Nm), eps (Nr,). All CGS."""
    import dustpy.constants as c
    warnings = []
    Sig_d_tot = snap["Sigma_dust"].sum(axis=1)
    if snap["Sigma_gas"].max() > 1e5:
        warnings.append(f"[SANITY WARN] snap {snap_index}: Sigma_gas max "
                        f"{snap['Sigma_gas'].max():.2e} g/cm^2 unusually high")
    if snap["T_gas"].max() > 3000 or snap["T_gas"].min() < 5:
        warnings.append(f"[SANITY WARN] snap {snap_index}: T range "
                        f"[{snap['T_gas'].min():.0f}, {snap['T_gas'].max():.0f}] K unusual")
    if Sig_d_tot.max() > 1e4:
        warnings.append(f"[SANITY WARN] snap {snap_index}: total dust Sigma "
                        f"max {Sig_d_tot.max():.2e} g/cm^2 — check sources/floor")
    if snap["St"].max() > 10:
        warnings.append(f"[SANITY WARN] snap {snap_index}: St_max "
                        f"{snap['St'].max():.2e} > 10 — drift-instability regime")
    if snap["eps"].max() > 1.0:
        warnings.append(f"[SANITY WARN] snap {snap_index}: eps_max "
                        f"{snap['eps'].max():.2f} > 1 — streaming instability likely")
    if snap["a"].max() > 1e3:   # > 10 m
        warnings.append(f"[SANITY WARN] snap {snap_index}: a_max "
                        f"{snap['a'].max():.2e} cm > 10 m — possible numerical runaway")
    return warnings
```

Typical ranges:

| Quantity | Typical | Notes |
|---|---|---|
| `gas.Sigma` | 1 – 10⁴ g cm⁻² | higher in the inner disk |
| `gas.T` | 10 – 2000 K | passively irradiated midplane |
| total dust `Sigma` | 0.01 – 100 g cm⁻² | mass-integrated |
| `dust.St` | 10⁻⁶ – 1 | warn above 1 |
| `dust.eps` | 10⁻⁴ – 0.1 (bulk) | warn above 1 |
| `dust.a` peak | 0.1 µm – 10 cm | fragmentation-limited |
| `dust.v.frag` | 100 – 1000 cm/s | default 100 cm/s |

---

## Customisation

Detailed, verified patterns live in `references/customisation.md`;
ready-made implementations of several of them live in **DustPyLib**
(`references/dustpylib.md`) — prefer DustPyLib when it is installed.
Quick index:

- Replace any physics function: `sim.<field>.updater = my_func` (before or
  after `initialize()`, then `sim.update()`); freeze with `= None`.
- Custom/refined radial grid: assign a plain `np.ndarray` to `sim.grid.ri`
  (cell interfaces) before `initialize()`. The mass grid must stay
  logarithmic — configure via `sim.ini.grid`, never set `sim.grid.m`.
  **DustPyLib:** `from dustpylib.grid.refinement import refine_radial_local`.
- Turn off fragmentation: `sim.dust.p.stick = 1.; sim.dust.p.frag = 0.;
  sim.dust.p.updater = None`.
- Turn off all coagulation: both probabilities to 0 + unset updaters of
  `v.rel`, `v.frag`, `kernel`.
- Turn off dust drift/diffusion: `sim.dust.v.rad = 0.; .updater = None`,
  `sim.dust.D = 0.; .updater = None` (docs recommend these over zeroing fluxes).
- Freeze gas: `sim.gas.nu[:] = 0.; sim.gas.nu.updater = None`.
- External sources/sinks (infall, winds, planetesimal formation):
  `sim.dust.S.ext.updater = f(sim)` / `sim.gas.S.ext.updater = f(sim)` [g cm⁻² s⁻¹].
- Planetary gaps: gap profile (Kanagawa 2017 / Duffell 2020) into
  `sim.gas.torque.Lambda` **or** the inverse gap profile onto `sim.gas.alpha`;
  refine the grid with `sim.grid.ri`.
  **DustPyLib:** `from dustpylib.substructures.gaps import kanagawa2017, duffell2020`.
- Ice lines: make `sim.dust.v.frag` depend on `sim.gas.T` (after
  `initialize()`), then `sim.update()`.
- Planetesimal formation sink: negative `sim.dust.S.ext` triggered by a
  midplane dust/pebble-to-gas threshold.
  **DustPyLib:** `from dustpylib.planetesimals.formation import drazkowska2016,
  schoonenberg2018, miller2021`.
- Dust backreaction (dust dynamical feedback on the gas, pull/push A, B).
  **DustPyLib:** `from dustpylib.dynamics.backreaction import setup_backreaction`
  (call after `initialize()`).
- Radiative transfer from a DustPy model (images/spectra).
  **DustPyLib:** `dustpylib.radtrans.radmc3d` (RADMC-3D inputs/output) and
  `dustpylib.radtrans.slab` (fast plane-parallel slab model).
- Track evolving planet masses: custom group via `sim.addgroup(...)` /
  `sim.<group>.addfield(...)` + updaters.
- Per-update hooks: systoles/diastoles via
  `sim.<field>.updater.systole = fn` / `sim.<field>.updater.diastole = fn`.
- Changing integrator: `sim.setdustintegrator(scheme="explicit", method="cash-karp")`
  (default: implicit 1st-order direct solver).

---

## Performance notes

- `Nmbpd` (mass bins per decade) is the main cost knob. `>= 7` is required for
  accuracy (Drążkowska et al. 2014); 10 roughly doubles the coagulation cost.
- Coagulation cost scales ~O(Nr × Nm²); gas evolution ~O(Nr).
- DustPy is single-process (no MPI); use job arrays for parameter surveys.
- Budget: the default model (100 cells, 100 kyr) runs in a few minutes on a
  modern single core; 1 Myr full-coagulation runs take ~30–60 min.
- The implicit dust integrator does a sparse LU solve each step; if the matrix
  is weakly conditioned and mass conservation degrades, switch to the explicit
  Cash-Karp scheme.

---

## Iron rules

1. **CGS only.** Never mix units silently; convert with `dustpy.constants`.
2. **Never overwrite existing output** (`sim.writer.overwrite = False`) without
   explicit user confirmation.
3. **Never assume parameters.** Read them from the setup script or the HDF5
   file (grid, mass range, t, star, gas fields are all stored).
4. **Read one snapshot first**, sanity-check, then batch-process.
5. `Nmbpd >= 7`, `rmin > 0`, `t_end > t_start`, `mmin < mmax`, `mmax` above the
   fragmentation barrier, `alpha` in (0, 0.1].
6. Mass grid is strictly logarithmic: set it via `sim.ini.grid.{mmin,mmax,Nmbpd}`,
   never by writing to `sim.grid.m`.
7. Set `sim.ini.*` **before** `sim.initialize()`; writer/snapshot settings
   **after** it (they need `sim.t` to exist).
8. After changing updaters or fields, call `sim.update()` before `sim.run()`.
9. Remind the user to cite Stammler & Birnstiel (2022) in any published use,
   plus the per-module reference(s) when DustPyLib modules are used.

---

## References

- `references/parameters.md` — full `sim.ini` parameter table: types, CGS
  defaults, constraints, cross-parameter interactions.
- `references/customisation.md` — verified customisation patterns with code.
- `references/dustpylib.md` — the DustPyLib companion library: grid
  refinement, planetary gaps, planetesimal formation, dust backreaction,
  RADMC-3D and slab radiative transfer (per-module API, usage, citations).
- `references/output.md` — output file layout, readers, resuming, units
  conventions for reporting.
- Docs: https://stammler.github.io/dustpy/ (basics, customisation, standard
  model, coagulation, examples: ice lines, planetary gaps, planetesimals)
- DustPyLib docs: https://dustpylib.readthedocs.io/ (repo:
  https://github.com/stammler/dustpylib)
- Simframe (framework): https://simframe.rtfd.io/
- Paper: https://doi.org/10.3847/1538-4357/ac7d58
