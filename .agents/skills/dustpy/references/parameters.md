# DustPy Parameters Reference

All parameters in `sim.ini` must be set **before** `sim.initialize()`.
All quantities are in **CGS units** unless stated otherwise. Defaults are
from `dustpy.simulation.Simulation.__init__` (v1.0.9).

Table of contents:
1. [Stellar parameters — `sim.ini.star`](#1-stellar-parameters)
2. [Grid parameters — `sim.ini.grid`](#2-grid-parameters)
3. [Gas parameters — `sim.ini.gas`](#3-gas-parameters)
4. [Dust parameters — `sim.ini.dust`](#4-dust-parameters)
5. [Time / snapshot parameters — `sim.t`](#5-time--snapshot-parameters)
6. [Writer parameters — `sim.writer`](#6-writer-parameters)
7. [Cross-parameter constraints and interactions](#7-cross-parameter-constraints-and-interactions)

---

## 1. Stellar parameters

Namespace: `sim.ini.star`

| Parameter | Type | Default (CGS) | Default (physical) | Constraints | Notes |
|---|---|---|---|---|---|
| `M` | float | `1.9885e+33` g | 1 M☉ | > 0 | Sets Keplerian frequency and dynamical timescales throughout the disk |
| `R` | float | `2 × 6.96e+10` cm | 2 R☉ | > 0 | Used with `T` to compute luminosity `L = 4πR²σT⁴` for the default passively irradiated temperature profile; has no effect if `gas.T.updater` is replaced |
| `T` | float | `5772.0` K | solar | > 0 | Effective temperature; enters the irradiation temperature profile together with `R` |

Derived (computed by `initialize()`):

| Field | Unit | Description |
|---|---|---|
| `sim.star.L` | erg/s | Bolometric luminosity, `L = 4πR²σ_SB T⁴` (updater: `std.star.luminosity`) |

---

## 2. Grid parameters

Namespace: `sim.ini.grid`

| Parameter | Type | Default (CGS) | Default (physical) | Constraints | Notes |
|---|---|---|---|---|---|
| `Nr` | int | `100` | — | ≥ 5 (practical) | Number of radial grid cells. Runtime scales ~O(Nr) for gas, ~O(Nr·Nm²) for coagulation. Production runs: 100–300 |
| `rmin` | float | `1.496e+13` cm | 1 AU | > 0, < `rmax` | Inner radial boundary. Must not be 0 (Keplerian singularity) |
| `rmax` | float | `1.496e+16` cm | 1000 AU | > `rmin` | Outer radial boundary |
| `Nmbpd` | int | `7` | — | **≥ 7 (hard)** | Mass bins per decade. Lower values make the coagulation mass-redistribution inaccurate (Drążkowska et al. 2014). 10 roughly doubles coagulation cost |
| `mmin` | float | `1e-12` g | ~1 nm grain | > 0, < `mmax` | Smallest particle mass on the mass grid |
| `mmax` | float | `1e+5` g | ~46 cm boulder | > `mmin` | Largest particle mass. Must exceed the expected fragmentation barrier, otherwise the top bin acts as an artificial sink |

Derived (computed by `makegrids()`; read-only afterwards):

| Field | Shape | Unit | Description |
|---|---|---|---|
| `sim.grid.r` | `(Nr,)` | cm | Radial cell centres (log-spaced by default) |
| `sim.grid.ri` | `(Nr+1,)` | cm | Radial cell interfaces; **can be assigned a plain `np.ndarray` before `initialize()` for a custom grid** |
| `sim.grid.Nm` | scalar | — | Number of mass bins = `ceil(log10(mmax/mmin)) × Nmbpd + 1` |
| `sim.grid.m` | `(Nm,)` | g | Mass grid (strictly log-spaced; never modify directly) |
| `sim.grid.A` | `(Nr,)` | cm² | Annulus area per radial cell, `π(ri[i+1]² − ri[i]²)` |
| `sim.grid.OmegaK` | `(Nr,)` | 1/s | Keplerian angular velocity (updater: `std.grid.OmegaK`) |

**Custom radial grid:** assign `sim.grid.ri` (cell interfaces, shape `Nr+1`)
before `sim.initialize()`; centres are the geometric means of adjacent
interfaces. For a local-refinement helper use
`dustpylib.grid.refinement.refine_radial_local` (see `dustpylib.md` §1) or
the hand-rolled version in `customisation.md` §2.

---

## 3. Gas parameters

Namespace: `sim.ini.gas`

| Parameter | Type | Default (CGS) | Default (physical) | Constraints | Notes |
|---|---|---|---|---|---|
| `alpha` | float | `1e-3` | — | 0 < α ≤ 0.1 | Shakura–Sunyaev viscosity: `ν = α c_s H_p`. Also the default for `dust.delta.{rad,turb,vert}` unless those are set independently |
| `Mdisk` | float | `0.05 M_sun` ≈ `9.94e+31` g | 0.05 M☉ | > 0 | Initial gas disk mass; normalises the Lynden-Bell & Pringle (1974) surface density profile |
| `mu` | float | `2.3 m_p` ≈ `3.85e-24` g | — | > 0 | Mean molecular weight (H₂/He mixture); affects sound speed and pressure scale height |
| `SigmaExp` | float | `-1.0` | — | any | Power-law index of the initial Σ_gas profile (LB&P 1974). Typical: −0.5 to −1.5 |
| `SigmaRc` | float | `60 × c.au` ≈ `8.98e+15` cm | 60 AU | > 0 | Cut-off radius of the self-similar profile; controls where the exponential taper starts |

Initial surface density profile (updater: `std.gas.lyndenbellpringle1974`):

```
Sigma(r) = (2+p) M_disk / (2π rc²) (r/rc)^p exp(-(r/rc)^(2+p))
```

Derived gas fields (created by `initialize()`):

| Field | Shape | Unit | Default updater | Description |
|---|---|---|---|---|
| `sim.gas.Sigma` | `(Nr,)` | g/cm² | — (evolved) | Gas surface density (integrated variable) |
| `sim.gas.SigmaFloor` | `(Nr,)` | g/cm² | constant 1e-100 | Floor preventing Σ → 0 |
| `sim.gas.T` | `(Nr,)` | K | `std.gas.T_passive` | Midplane temperature (passive irradiation, θ_irr = 0.05) |
| `sim.gas.cs` | `(Nr,)` | cm/s | `std.gas.cs_isothermal` | Isothermal sound speed |
| `sim.gas.Hp` | `(Nr,)` | cm | `std.gas.Hp` | Pressure scale height |
| `sim.gas.nu` | `(Nr,)` | cm²/s | `std.gas.nu` | Kinematic viscosity `α c_s H_p` |
| `sim.gas.eta` | `(Nr,)` | — | `std.gas.eta_midplane` | Pressure gradient parameter (sub-Keplerian measure) |
| `sim.gas.rho` | `(Nr,)` | g/cm³ | `std.gas.rho_midplane` | Midplane volume density |
| `sim.gas.n` | `(Nr,)` | 1/cm³ | `std.gas.n_midplane` | Midplane number density |
| `sim.gas.mfp` | `(Nr,)` | cm | `std.gas.mfp_midplane` | Mean free path (selects Epstein/Stokes drag) |
| `sim.gas.P` | `(Nr,)` | dyn/cm² | `std.gas.P_midplane` | Midplane pressure |
| `sim.gas.Fi` | `(Nr+1,)` | g/cm/s | `std.gas.Fi` | Gas fluxes at cell interfaces |
| `sim.gas.v.visc` | `(Nr,)` | cm/s | `std.gas.vvisc` | Viscous accretion velocity |
| `sim.gas.v.rad` | `(Nr,)` | cm/s | `std.gas.vrad` | Radial velocity (incl. backreaction + torque) |
| `sim.gas.S.ext` | `(Nr,)` | g/cm²/s | none (0) | External sources (infall, wind, photoevaporation); boundary cells overwritten by BCs |
| `sim.gas.S.hyd` | `(Nr,)` | g/cm²/s | `std.gas.S_hyd` | Hydrodynamic sources |
| `sim.gas.S.tot` | `(Nr,)` | g/cm²/s | `std.gas.S_tot` | Total sources |
| `sim.gas.torque.Lambda` | `(Nr,)` | cm²/s² | none (0) | Specific angular momentum injection rate (planetary gaps) |
| `sim.gas.torque.v` | `(Nr,)` | cm/s | `std.gas.vtorque` | Effective velocity from torque, `2Λ/(Ω_K r)` |

---

## 4. Dust parameters

Namespace: `sim.ini.dust`

### 4a. Initial conditions

| Parameter | Type | Default (CGS) | Default (physical) | Constraints | Notes |
|---|---|---|---|---|---|
| `d2gRatio` | float | `0.01` | — | > 0 | Initial dust-to-gas mass ratio, applied uniformly to set initial `dust.Sigma` |
| `aIniMax` | float | `1e-4` cm | 1 µm | > 0 | Maximum grain size filled at t=0 (subject to `allowDriftingParticles`) |
| `distExp` | float | `-3.5` | −7/2 (MRN) | any | Power-law index of the initial size distribution `n(a)da ∝ a^distExp da`; default is the MRN ISM value |
| `allowDriftingParticles` | bool | `False` | — | — | If `False` (default), particles in the outer disk that would drift inward faster than one orbit are removed from the initial condition, preventing an artificial initial drift wave. Set `True` only when studying that wave |

### 4b. Collision physics (coagulation model)

The model is the sticking–erosion–fragmentation scheme of Brauer et al.
(2008) / Blum & Hahn (2007). Coagulation matrices are pre-computed at
`initialize()` (updater: `std.dust.coagulation_parameters`); parameters below
therefore must be set **before** `initialize()`.

| Parameter | Type | Default | Constraints | Notes |
|---|---|---|---|---|
| `vFrag` | float | `100.0` cm/s (1 m/s) | > 0 | Fragmentation velocity threshold; collisions above it fragment. Typical: 100 cm/s (silicates) to 1000 cm/s (icy grains). Can be made spatially varying via `sim.dust.v.frag` (see ice lines) |
| `rhoMonomer` | float | `1.67` g/cm³ | > 0 | Solid (bulk) density of a monomer; silicate default, use ~1.0 for icy aggregates. Enters the Stokes number via the mass–size relation |
| `erosionMassRatio` | float | `10.0` | > 1 | Mass ratio m_large/m_small above which only the small particle fragments (erosion); below it both fragment fully |
| `excavatedMass` | float | `1.0` | > 0 | Mass chipped from the larger particle in an erosive collision, in units of the smaller particle's mass |
| `fragmentDistribution` | float | `-11/6 ≈ -1.833` | < −1 | Power-law index of the fragment mass distribution `n(m)dm ∝ m^γ dm` (Dohnanyi 1969); more negative = more small fragments |

### 4c. Key derived dust fields (read-only after `initialize()`)

| Field | Shape | Unit | Description |
|---|---|---|---|
| `sim.dust.Sigma` | `(Nr, Nm)` | g/cm² | Dust surface density per mass bin — the primary evolved variable |
| `sim.dust.SigmaFloor` | `(Nr, Nm)` | g/cm² | Floor: one particle per annulus per mass bin, times 0.1 |
| `sim.dust.a` | `(Nr, Nm)` | cm | Particle size per mass bin (updater: `std.dust.a`) |
| `sim.dust.rho` | `(Nr, Nm)` | g/cm³ | Midplane volume mass density per mass bin |
| `sim.dust.rhos` | `(Nr, Nm)` | g/cm³ | Solid-state (bulk) density per mass bin (initialized to `rhoMonomer`) |
| `sim.dust.St` | `(Nr, Nm)` | — | Stokes number (Epstein + Stokes I drag, `std.dust.St_Epstein_StokesI`) |
| `sim.dust.eps` | `(Nr,)` | — | Dust-to-gas ratio, `Σ_dust,Σgas` (mass-integrated) |
| `sim.dust.H` | `(Nr, Nm)` | cm | Dust scale height (Dubrulle et al. 1995, `std.dust.H`) |
| `sim.dust.D` | `(Nr, Nm)` | cm²/s | Dust diffusivity (`std.dust.D`; zero in the 2 outermost/innermost cells) |
| `sim.dust.v.rad` | `(Nr, Nm)` | cm/s | Radial drift velocity per mass bin |
| `sim.dust.v.frag` | `(Nr,)` | cm/s | Fragmentation velocity (make spatially varying for ice lines) |
| `sim.dust.v.driftmax` | `(Nr,)` | cm/s | Maximum drift velocity incl. backreaction |
| `sim.dust.v.rel.{azi,brown,rad,turb,vert,tot}` | `(Nr, Nm, Nm)` | cm/s | Relative-velocity contributions and total |
| `sim.dust.Fi.{adv,diff,tot}` | `(Nr+1, Nm)` | g/cm/s | Fluxes at interfaces |
| `sim.dust.kernel` | `(Nr, Nm, Nm)` | cm²/s | Vertically integrated collision kernel |
| `sim.dust.p.{stick,frag}` | `(Nr, Nm, Nm)` | — | Probability matrices (linear transition at `v.frag`) |
| `sim.dust.S.{coag,hyd,ext,tot}` | `(Nr, Nm)` | g/cm²/s | Source terms; `S.ext` has no updater by default |
| `sim.dust.backreaction.{A,B}` | `(Nr,)` | — | Pull/push backreaction coefficients (A=1, B=0 by default, i.e. off). Enable with DustPyLib `setup_backreaction` — see `dustpylib.md` §4 |

### 4d. Mixing / diffusion parameters

Set directly on the frame (default: all equal to `gas.alpha`); set before or
after `initialize()` — set them once and freeze, or give an updater:

| Field | Default | Notes |
|---|---|---|
| `sim.dust.delta.rad` | `alpha` | Radial diffusion strength: `D = δ_rad α-like · v²` enters the diffusive flux |
| `sim.dust.delta.turb` | `alpha` | Turbulent relative-velocity strength (Ormel & Cuzzi 2007); sets the coagulation/fragmentation boundary |
| `sim.dust.delta.vert` | `alpha` | Vertical settling parameter: `H_dust = H_gas / sqrt(1 + St/δ_vert)` |

---

## 5. Time / snapshot parameters

`sim.t.snapshots`: 1-D array of absolute output times in **seconds**. The
first element must be the current time `sim.t` (0 before `initialize()`).
Set after `initialize()`.

```python
# 10 snapshots per decade from 1 kyr to 1 Myr
sim.t.snapshots = np.hstack([
    sim.t,
    np.geomspace(1e3, 1e6, num=31) * c.year   # 31 points = 10/decade over 3 decades
])
```

Computing `num` from `snaps_per_decade`:

```python
num = int(snaps_per_decade * np.log10(t_end_yr / t_start_yr)) + 1
```

Guidelines: start snapshots at ≥ 1 orbit at `rmin` (avoids initialisation
transients); 5–10 per decade is usually sufficient.

---

## 6. Writer parameters

Set on `sim.writer` (a simframe `hdf5writer` instance created during
`initialize()`):

| Attribute | Type | Default | Notes |
|---|---|---|---|
| `sim.writer.datadir` | str | `"data"` | Output directory; created if absent |
| `sim.writer.overwrite` | bool | `False` | If `False`, raises when the directory already contains output — the safe default. Set `True` only after explicit user confirmation |
| `sim.writer.filename` | str | `"data"` | File stem: `data0000.hdf5`, `data0001.hdf5`, ... |

Output written per snapshot: all (non-hidden) frame fields, e.g.
`/t`, `/star/M`, `/grid/r`, `/grid/m`, `/gas/Sigma`, `/gas/T`,
`/dust/Sigma`, `/dust/St`, `/dust/a`, `/dust/eps`, plus the binary dump
`frame.dmp` (restored with `dustpy.readdump`).

Exclude individual fields from the output with `sim.<field>.save = False`
(e.g. `sim.dust.kernel.save = False` to save file size — the kernel is
recomputed anyway).

---

## 7. Cross-parameter constraints and interactions

**`Nmbpd ≥ 7` is a hard requirement.** With fewer bins per decade the mass
redistribution in sticking events (into the two neighbouring bins) becomes
inaccurate; results can differ by orders of magnitude. Ref: Drążkowska et
al. (2014).

**`mmax` must exceed the expected fragmentation barrier.** A rough
fragmentation-limited size estimate:

```
a_frag ≈ (2/3π) (Σ_gas / (ρ_s π)) (v_frag / c_s)² (1/α)
```

If `mmax` is too small the top mass bin acts as a sink. Default `1e5` g is
safe for typical disks.

**`alpha` is shared by gas viscosity and dust mixing by default.** All three
`dust.delta.*` default to `gas.alpha`. Decouple explicitly if needed, e.g.

```python
sim.initialize()
sim.dust.delta.rad[:]  = 1e-4
sim.dust.delta.rad.updater = None
```

**`vFrag` and `rhoMonomer` set the fragmentation barrier location.** Higher
`vFrag` or lower `rhoMonomer` push the barrier to larger grains;
`St_frag ∝ v_frag² / (α c_s²)`.

**Erosion parameters are baked in at `initialize()`.** `erosionMassRatio`,
`excavatedMass`, `fragmentDistribution` enter the pre-computed coagulation
matrices; changing them after `initialize()` has no effect until re-init.

**`SigmaRc` interacts with `rmin`/`rmax`.** If `SigmaRc` ≪ `rmax`, the
exponential taper empties the outer disk; rule of thumb
`rmin ≲ 0.1·SigmaRc ≲ rmax`.

**`allowDriftingParticles = False` (default) trims the initial outer-dust
condition** so that no particle drifts inward on a time shorter than one
orbit initially. Set `True` only when studying that initial drift wave.

**Backreaction is off by default** (`backreaction.A = 1`, `B = 0`).
Enable with the companion library: `from dustpylib.dynamics.backreaction import
setup_backreaction; setup_backreaction(sim)` after `initialize()` (optionally
`vertical_setup=True`). This wires up the A, B coefficients, damps the dust
diffusivity, and `gas.v.rad` (which already consumes A, B) and `dust.v.driftmax`
pick them up automatically. See `dustpylib.md` §4 for details.
