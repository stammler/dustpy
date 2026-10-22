# DustPyLib Reference

DustPyLib is the official companion library of DustPy: ready-made, well-tested
extensions (dust backreaction, grid refinement, planetesimal formation,
planetary gaps, radiative-transfer interfaces). All content here is verified
against the module source (v0.7.0, `stammler/dustpylib`).

- Install: `pip install dustpylib` (Python 3; pulls in `dustpy`, `dsharp_opac`,
  `astropy`, `matplotlib`, `numpy`, `scipy`, `simframe`)
- Docs: https://dustpylib.readthedocs.io/
- Repo: https://github.com/stammler/dustpylib
- **Citations:** every module requires citing its own publication(s) *in
  addition to* DustPy. The module sources below list the exact references;
  `setup_backreaction` and the opacity writers also print the citations to
  run. Never skip the per-module citation when reporting results.

Table of contents:
1. [Grid refinement — `dustpylib.grid.refinement`](#1-grid-refinement)
2. [Planetary gaps — `dustpylib.substructures.gaps`](#2-planetary-gaps)
3. [Planetesimal formation — `dustpylib.planetesimals.formation`](#3-planetesimal-formation)
4. [Dust backreaction — `dustpylib.dynamics.backreaction`](#4-dust-backreaction)
5. [RADMC-3D interface — `dustpylib.radtrans.radmc3d`](#5-radmc-3d-interface)
6. [Slab model — `dustpylib.radtrans.slab`](#6-slab-model)

All quantities CGS. All usage examples below are taken from the official
example notebooks in the DustPyLib repo.

---

## 1. Grid refinement

`dustpylib.grid.refinement.refine_radial_local(ri, r0, num=3)`

Refines the radial grid **locally** around `r0` by recursively splitting
cells (each level doubles the cells in a window of `num` cells around the
target; insertions are geometric means, so the grid stays log-spaced).

| Arg | Meaning |
|---|---|
| `ri` | radial grid cell interfaces, `(N,)` cm |
| `r0` | location to refine around, cm |
| `num` | refinement levels (default 3); each level doubles cells in the window |

Returns the refined interface array. **Use before `initialize()`** by
assigning the result to `sim.grid.ri` (see `customisation.md` §2):

```python
from dustpy import Simulation
import dustpy.constants as c
import numpy as np
from dustpylib.grid.refinement import refine_radial_local

s = Simulation()
# default grid cell interfaces (identical to what initialize() would build)
ri = np.geomspace(s.ini.grid.rmin, s.ini.grid.rmax, s.ini.grid.Nr)

ri = refine_radial_local(ri, 10. * c.au, num=3)   # refine around 10 AU
# several planets:
for a in [5.2 * c.au, 9.6 * c.au]:
    ri = refine_radial_local(ri, a, num=3)

s.grid.ri = ri
s.initialize()                                     # s.grid.Nr is now larger
```

Notes:
- Call it repeatedly to refine around several locations.
- The mass grid is unaffected — it stays strictly logarithmic.
- A hand-rolled equivalent is given in `customisation.md` §2 if DustPyLib
  is not installed.

---

## 2. Planetary gaps

`dustpylib.substructures.gaps` provides two analytical gap profiles
(fitness of 1D surface-density perturbation `f(r) = Σ(r)/Σ_0`, tending to 1
outside the gap — `kanagawa2017` has hard gap edges, `duffell2020` decays
smoothly/asymptotically to 1):

| Function | Signature | Reference |
|---|---|---|
| `kanagawa2017` | `(r, a, q, h, alpha0)` | Kanagawa et al. (2017), PASJ 69, 135, doi:10.1093/pasj/psx114 |
| `duffell2020` | `(r, a, q, h, alpha0)` | Duffell (2020), ApJ 897, 15, doi:10.3847/1538-4357/ab5b0f |

`r` may be the radial grid (or interfaces), `a` the planet semi-major axis,
`q = M_planet/M_star`, `h = Hp/r` at the planet, `alpha0` the unperturbed
turbulent viscosity.

### Recommended pattern: inverse gap on α

Because `ν·Σ = const` in steady state, imposing the *inverse* gap profile on
`alpha` produces the gap in the gas (and dust) surface density — and the gas
**still evolves** (unlike setting a static Σ profile). From the official
notebook:

```python
from dustpy import Simulation
import dustpy.constants as c
import numpy as np
from dustpylib.grid.refinement import refine_radial_local
from dustpylib.substructures.gaps import kanagawa2017
from scipy.interpolate import interp1d

planets = {
    "jupiter": {"a": 5.2038 * c.au, "M": 317.8 * c.M_earth},
    "saturn":  {"a": 9.5826 * c.au, "M": 95.159 * c.M_earth},
    "uranus":  {"a": 19.19126 * c.au, "M": 14.536 * c.M_earth},
    "neptune": {"a": 30.07 * c.au, "M": 17.147 * c.M_earth},
}

s = Simulation()
ri = np.geomspace(s.ini.grid.rmin, s.ini.grid.rmax, s.ini.grid.Nr)
for planet in planets.values():
    ri = refine_radial_local(ri, planet["a"], num=3)
s.grid.ri = ri
s.initialize()

s.addgroup("planets", description="Planets")
for name, planet in planets.items():
    s.planets.addgroup(name, description=name.title())
    s.planets.__dict__[name].addfield("M", planet["M"], description="Mass in g")
    s.planets.__dict__[name].addfield("a", planet["a"], description="Semi-major axis in cm")

alpha0 = s.gas.alpha.copy()          # unperturbed (constant) alpha

def alpha(s):
    a = alpha0.copy()
    for name, planet in s.planets.__dict__.items():
        if name.startswith("_"):
            continue
        q = planet.M / s.star.M
        f_h = interp1d(s.grid.r, s.gas.Hp / s.grid.r)
        h = f_h(planet.a)
        f_alp = interp1d(s.grid.r, alpha0)
        alp = f_alp(planet.a)
        a /= kanagawa2017(s.grid.r, planet.a, q, h, alp)   # inverse gap
    return a

s.gas.alpha.updater = alpha
s.update()

# Rescale the initial surface densities so nu*Sigma stays constant:
s.gas.Sigma[...]  /= s.gas.alpha / alpha0
s.dust.Sigma[...] /= (s.gas.alpha / alpha0)[:, None]
```

Notes:
- `duffell2020(...)` is a drop-in replacement for `kanagawa2017(...)` in the
  same place (note it assumes a *constant* `alpha0`).
- Alternative (DustPy-only, no DustPyLib): the torque-based pattern in
  `customisation.md` §6, which writes the gap gradient into
  `sim.gas.torque.Lambda`. Use it when DustPyLib is unavailable.
- Cite Kanagawa et al. (2017) or Duffell (2020) accordingly.

---

## 3. Planetesimal formation

`dustpylib.planetesimals.formation` provides three prescriptions as dust
*source terms* (negative = dust removed). All share the signature:

```python
S = func(OmegaK, rho_dust, rho_gas, Sigma_dust, St, **params)
# returns (Nr, Nm) array in g cm^-2 s^-1; boundary rows [0] and [-1] zeroed
```

| Function | Extra params (defaults) | Mechanism | Citation |
|---|---|---|---|
| `drazkowska2016` | `p2g_crit=1.`, `St_crit=0.01`, `zeta=0.01` | Only pebbles (St ≥ St_crit) count; if midplane pebble-to-gas ratio ≥ p2g_crit, a fraction ζ of pebbles is converted per orbit | Drążkowska et al. (2016), A&A 595, A45, doi:10.1051/0004-6361/201628983 |
| `schoonenberg2018` | `d2g_crit=1.`, `zeta=0.1` | All particles count; if midplane dust-to-gas ≥ d2g_crit, a fraction ζ is converted per settling time 1/(St·Ω_K) (small particles suppressed) | Schoonenberg et al. (2018), A&A 619, A119, doi:10.1051/0004-6361/201834047 |
| `miller2021` | `d2g_crit=1.`, `n=0.03`, `zeta=0.1` | Smooth tanh transition in midplane dust-to-gas (no hard trigger); conversion rate ∝ P·ζ·St·Ω_K | Miller et al. (2021), MNRAS 508, 1558, doi:10.1093/mnras/stab2935 |

### Wiring into a simulation (official pattern)

```python
from dustpy import Simulation
from dustpylib.planetesimals.formation import drazkowska2016
import numpy as np
from simframe import Instruction
from simframe import schemes

p2g_crit, St_crit, zeta = 1., 0.01, 0.01

s = Simulation()
s.initialize()

def S_ext(s):
    return drazkowska2016(
        s.grid.OmegaK,
        s.dust.rho,      # (Nr, Nm) g/cm^3 midplane
        s.gas.rho,       # (Nr,) g/cm^3
        s.dust.Sigma,    # (Nr, Nm) g/cm^2
        s.dust.St,       # (Nr, Nm)
        p2g_crit=p2g_crit, St_crit=St_crit, zeta=zeta)

s.dust.S.ext.updater = S_ext

# --- bookkeeping: track the converted planetesimal mass ---
s.addgroup("planetesimals", description="Planetesimal quantities")
s.planetesimals.addfield("Sigma", np.zeros_like(s.gas.Sigma),
                         description="Surface density of planetesimals [g/cm^2]")

def dSigma_planetesimals(s, t, Sigma_planetesimals):
    return -s.dust.S.ext.sum(-1)

s.planetesimals.Sigma.differentiator = dSigma_planetesimals
s.integrator.instructions.append(
    Instruction(schemes.expl_1_euler, s.planetesimals.Sigma,
                description="Planetesimals: explicit 1st-order Euler"))
s.update()
```

Notes:
- `schoonenberg2018` / `miller2021` are wired identically — only the
  function and its parameters change.
- The planetesimal mass field is 1-D (`Nr,`); the differentiator consumes
  the negative dust source summed over mass bins, so total dust+planetesimal
  mass is conserved.
- Cite the chosen prescription (and DustPy).
- The simple hand-rolled sink in `customisation.md` §8 is a minimal
  alternative when DustPyLib is unavailable.

---

## 4. Dust backreaction

`dustpylib.dynamics.backreaction.setup_backreaction(sim, vertical_setup=False)`

One-line setup of the dust dynamical feedback (pull/push coefficients A, B of
Garaté et al. 2019, 2020). **Call after `initialize()`, before `run()`:**

```python
from dustpy import Simulation
from dustpylib.dynamics.backreaction import setup_backreaction

sim = Simulation()
sim.ini.dust.d2gRatio = 0.01     # key parameters: d2g, vFrag, alpha, delta.*
sim.ini.dust.vFrag     = 1000.
sim.initialize()

setup_backreaction(sim)                  # standard (vertically uniform d2g)
# or
setup_backreaction(sim, vertical_setup=True)   # accounts for dust settling
sim.run()
```

What the setup does (verified in source):
- Standard mode: `sim.dust.backreaction.updater = BackreactionCoefficients` —
  computes per radius, summing over mass bins
  `X = Σ d2g_i/(1+St_i²)`, `Y = Σ St_i·d2g_i/(1+St_i²)`,
  `A = (X+1)/(X+1)² + Y²`, `B = Y/(X+1)² + Y²` (A→1, B→0 when dust-free).
- `vertical_setup=True`: uses a local vertical grid (Gaussian gas/dust
  structure) to compute X, Y per radius *and height*; stores per-species
  factors in new fields `sim.dust.backreaction.A_dust_settling` /
  `B_dust_settling` (`Nr, Nm`) and replaces `sim.dust.v.rad.updater` with
  `vrad_dust_BackreactionVerticalStructure`.
- Always: `sim.dust.D.updater = dustDiffusivity_Backreaction` — damps the
  dust diffusivity by `1/(1 + d2g_tot)` (important at high dust-to-gas).
- The outermost radial cell is pinned to (A=1, B=0) for stability.
- Effect on the gas (wired through DustPy's existing
  `sim.dust.backreaction.{A,B}` fields consumed by `gas.v.rad`):
  `v_gas = A·v_visc + 2B·v_η`, `Δv_φ,gas = −A·v_η + (B/2)·v_visc`.

Notes:
- Cite Garaté et al. (2019), ApJ 871, 53, and Garaté et al. (2020), A&A 635,
  A149. The function prints the reminder.
- With backreaction active, do **not** use the "gas hydrodynamics off"
  pattern of `customisation.md` §4 (zeroing `nu` while A≠1, B≠0 is
  inconsistent).
- Key controls: initial `d2gRatio`, `vFrag`, `gas.alpha`, and
  `sim.dust.delta.turb` (turbulent relative velocities).

---

## 5. RADMC-3D interface

`dustpylib.radtrans.radmc3d` — creates axisymmetric RADMC-3D input files from
a DustPy snapshot and reads RADMC-3D output back.

```python
from dustpylib.radtrans import radmc3d
# exported: Model, read_image, read_model, read_spectrum
```

### Writing an input model

```python
from dustpy import hdf5writer
from dustpylib.radtrans import radmc3d
import numpy as np

writer = hdf5writer()
writer.datadir = "data/my_run"
data = writer.read.output(21)        # one snapshot (namespace)

rt = radmc3d.Model(data)             # or a live Simulation
rt.phii_grid = np.array([0., 2.*np.pi])            # 1 phi cell (azimuthally uniform)
rt.ai_grid = np.geomspace(rt.a_dust_.min(), 1., 17)  # particle-size interfaces [cm]
rt.radmc3d_options["nphot"]        = 500_000
rt.radmc3d_options["nphot_scat"]   = 500_000
rt.radmc3d_options["nphot_spec"]   = 10_000
rt.radmc3d_options["mc_scat_maxtauabs"] = 5.
rt.radmc3d_options["dust_2daniso_nphi"]   = 60
rt.datadir = "radmc3d"
rt.write_files()                     # write_files(datadir=None, write_opacities=True,
                                     #              opacity=None, smoothing=False)
```

`Model(sim, ignore_last=True)`:
- `sim` may be a `dustpy.Simulation` **or** the namespace returned by
  `hdf5writer.read.output()` (this is the usual path — it works on saved
  data without a live run).
- `ignore_last=True` drops DustPy's outermost radial cell (its gas Σ sits on
  the floor and would distort the RT model).
- Imports `M_star_`, `R_star_`, `T_star_`, `T_gas_`, `a_dust_`, `H_dust_`,
  `Sigma_dust_` (trailing underscore = from DustPy).
- Grids: radial interfaces are auto-refined near the inner edge
  (`refine_radial_local(..., num=3)`); default wavelength grid is
  0.1 µm–1 mm (150 points); θ grid is 128 cells with fine resolution near the
  midplane; default `Nphi=16` and `Nspec=16` particle sizes
  (`geomspace(a_min, a_max)`).
- Settable grid attributes (interfaces only — the `*c_grid` centres are
  derived read-only properties that raise if set): `ri_grid`, `ai_grid`,
  `thetai_grid`, `phii_grid`.
- `opacity` attribute: `"birnstiel2018"` (default, DSHARP opacities) or
  `"ricci2010"`; opacities computed with the `dsharp_opac` package
  (opacities are per particle size, so the `ai_grid` controls the species).
- `write_files(..., smoothing=True)` averages opacities over several sizes
  (much slower); `write_opacities=False` skips the opacity files.
- Also writes `metadata.npz` (particle sizes, not needed by RADMC-3D itself
  but useful for post-processing).

Files written: `radmc3d.inp`, `stars.inp`, `wavelength_micron.inp`,
`amr_grid.inp`, `dust_density.inp`, `dust_temperature.dat`,
`dustopac.inp`, `dustkapscatmat.inp`, `metadata.npz`.

### Reading results back

```python
model = radmc3d.read_model(datadir="radmc3d")
# model.grid: r, ri, theta, thetai, phi, phii, a, ai
# model.rho:  (Nspecies, Nphi, Ntheta, Nr)  dust density
# model.T:    (Nspecies, Ntheta, Nphi, Nr)  (if dust_temperature.dat present;
#            after a run, copy dust_temperature_mctherm.dat over it first)

img = radmc3d.read_image("radmc3d/image_radio.out")
# img: dict with x, y [cm], lambda [cm], I/Q/U/V (if I-format 3)
# I has shape (Ny, Nx, Nlambda)

spec = radmc3d.read_spectrum("radmc3d/spectrum.out")
```

`read_model`/`read_image`/`read_spectrum` are only meant for models created
by this module; use `Radmc3dPy` for general RADMC-3D models.

Citations: RADMC-3D (Doll et al. 2009, A&A 507, 359, doi:10.1051/0004-6361/200810666)
plus whatever opacity model you used — Birnstiel et al. (2018), ApJL 864, L34
(doi:10.3847/2041-8213/aaf743) for the default, or Ricci et al. (2010), A&A 521, A321
(doi:10.1051/0004-6361/200913403) for `ricci2010`. The module prints the exact
references during the opacity computation.

---

## 6. Slab model

`dustpylib.radtrans.slab` — fast plane-parallel radiative-transfer estimate
of continuum emission from a DustPy model (intensity profiles, integrated
flux, effective size). Based on the analytical slab solution of Birnstiel et
al. (2018) (originally Miyake et al. 1993), with or without scattering.

```python
from dustpylib.radtrans import slab
# exported: Opacity, get_observables, get_all_observables
```

### Opacity

```python
opac = slab.Opacity()                      # default: "default_opacities_smooth.npz"
opac = slab.Opacity('default_opacities.npz')
opac = slab.Opacity('icefree_opacities_smooth.npz')
opac = slab.Opacity('ricci_compact.npz')
opac = slab.Opacity('dustkapabs.inp')      # or a RADMC-3D opacity .inp file
```

`Opacity(input=None, **kwargs)` — input can be a path to a `.npz`
(dsharp_opac format), a dsharp_opac data-file name (resolvable via
`dsharp_opac.get_datafile`; shipped names include `default_opacities`,
`default_opacities_smooth`, `default_opacities_extrapol`,
`icefree_opacities`, `icefree_opacities_smooth`, `ricci_compact`), or a
RADMC-3D `dustkap*.inp` file. Extra kwargs go to the log-log
`RegularGridInterpolator` (e.g. `bounds_error=True`).

Methods (all take `a` particle size [cm] and `lam` wavelength [cm]):
- `get_opacities(a, lam)` → `k_abs, k_sca` (cm²/g per cm)
- `get_g(a, lam)` → asymmetry parameter
- `get_k_ext_eff(a, lam)` → `k_abs + (1−g)·k_sca`
- `opac.rho_s` → material density [g/cm³]

Extrapolation outside the table is rejected unless `a > 100·λ/(2π)`
(optical-size regime), where it is safe.

### Computing observables

```python
lam = np.geomspace(0.05, 30., 20)           # cm — e.g. 5 mm to 30 cm radio window

obs = slab.get_all_observables(data, opac, lam)
```

- `data` for `get_all_observables` can be:
  - a **path string** to the DustPy data directory (internally
    `dustpy.hdf5writer().read.all()`), or
  - a `SimpleNamespace` from `hdf5writer.read.output()` (single snapshot) or
    `read.all()` (time series; fields must have time as first index).
- It loops over all times and returns arrays with a leading time axis:
  `obs.rf` [cm], `obs.flux_t` [Jy], `obs.tau`, `obs.I_nu` [Jy/arcsec²].
- Lower level: `slab.get_observables(r, sig_da, T, a, lam, opacity,
  scattering=True, inc=0.0, distance=140*pc, flux_fraction=0.68)` where
  `r` (Nr,) cm, `sig_da` (Nr, Na) g/cm², `T` (Nr,) K, `a` (Nr, Na) or
  (Na,) cm, `lam` (Nlam,) cm; `inc` is the inclination in radians
  (0 = face-on, increases path length).
- Module constants: `slab.pc` (cm), `slab.jy_sas` (Jy/arcsec² in
  erg s⁻¹ cm⁻² Hz⁻¹ sr⁻¹); `slab.bplanck(freq, temp)` for the Planck
  function [erg s⁻¹ cm⁻² Hz⁻¹ sr⁻¹].

Typical workflow (from the official notebook — comparing a radio image
against the slab model):

```python
from dustpy import hdf5writer
import dustpy.constants as c
from dustpylib.radtrans import slab

writer = hdf5writer()
writer.datadir = "data/my_run"
data = writer.read.output(21)

opac = slab.Opacity('default_opacities_smooth.npz')
lam = np.geomspace(0.05, 30., 20)          # cm — radio window 5 mm - 30 cm
obs = slab.get_all_observables(data, opac, lam)

# integrated flux [Jy] and effective radius per wavelength:
for i, lam_i in enumerate(lam):
    print(lam_i * 1.e4, "um:", obs.flux_t[0, i], "Jy, R_eff:", obs.rf[0, i] / c.au, "AU")
```

Citations: Birnstiel et al. (2018), ApJL 864, L34, doi:10.3847/2041-8213/aaf743
(slab solution + default opacities) and/or Miyake et al. (1993),
doi:10.1006/icar.1993.1156.

---

## Cross-cutting notes

- DustPyLib modules are **add-ons**: they never modify DustPy itself. The
  patterns in `customisation.md` (updaters, custom groups/instructions,
  `S.ext`) are the general mechanism these modules plug into.
- All DustPyLib setup happens **after** `sim.initialize()` (backreaction,
  planetesimal source, gap alpha-updater, rescaling of Σ) or **before**
  (`grid.ri` refinement). See the per-module sections above.
- If a module is missing from the environment (`ImportError`), the
  `customisation.md` fallbacks (hand-rolled refinement, Kanagawa torque
  pattern, plain eps-threshold sink) keep the task doable — but prefer the
  DustPyLib versions when available, and install with `pip install dustpylib`.
- For results using any DustPyLib module, cite DustPy **and** the module's
  reference(s) (listed per module above).
