# DustPy Output Reference

Everything about the files a DustPy run produces and how to read them back.
All stored values are in **CGS** (time in seconds, length cm, mass g).

---

## 1. Files

A run with `sim.writer.datadir = "data"` writes:

| File | Content |
|---|---|
| `data/data0000.hdf5`, `data0001.hdf5`, ... | One HDF5 file per snapshot, containing all non-hidden frame fields grouped by frame group (`/t`, `/star/...`, `/grid/...`, `/gas/...`, `/dust/...`) |
| `data/frame.dmp` | Binary simframe dump of the full frame (values **and** updater wiring); updated at each snapshot; used for resuming |

File names use the writer's `filename` stem (`data` by default) with a
4-digit zero-padded index and the `extension` (`hdf5`).

---

## 2. HDF5 layout (per snapshot)

| Path | Shape | Unit | Description |
|---|---|---|---|
| `t` | scalar | s | Snapshot time |
| `star/M` | scalar | g | Stellar mass |
| `star/R` | scalar | cm | Stellar radius |
| `star/T` | scalar | K | Effective temperature |
| `star/L` | scalar | erg/s | Luminosity |
| `grid/Nr`, `grid/Nm` | scalar | — | Grid sizes |
| `grid/r` | (Nr,) | cm | Radial cell centres |
| `grid/ri` | (Nr+1,) | cm | Radial cell interfaces |
| `grid/m` | (Nm,) | g | Mass grid (log) |
| `grid/OmegaK` | (Nr,) | 1/s | Keplerian frequency |
| `gas/Sigma` | (Nr,) | g/cm² | Gas surface density |
| `gas/SigmaFloor` | (Nr,) | g/cm² | Gas floor |
| `gas/T` | (Nr,) | K | Midplane temperature |
| `gas/cs` | (Nr,) | cm/s | Sound speed |
| `gas/Hp` | (Nr,) | cm | Pressure scale height |
| `gas/nu` | (Nr,) | cm²/s | Kinematic viscosity |
| `gas/rho` | (Nr,) | g/cm³ | Midplane volume density |
| `gas/eta` | (Nr,) | — | Pressure gradient parameter |
| `gas/Fi` | (Nr+1,) | g/cm/s | Interface gas fluxes |
| `dust/Sigma` | (Nr, Nm) | g/cm² | Dust surface density per mass bin |
| `dust/SigmaFloor` | (Nr, Nm) | g/cm² | Dust floor |
| `dust/a` | (Nr, Nm) | cm | Particle size per bin |
| `dust/St` | (Nr, Nm) | — | Stokes number per bin |
| `dust/eps` | (Nr,) | — | Dust-to-gas ratio (mass-integrated) |
| `dust/rho` | (Nr, Nm) | g/cm³ | Midplane volume density per bin |
| `dust/D` | (Nr, Nm) | cm²/s | Diffusivity |
| `dust/H` | (Nr, Nm) | cm | Dust scale height |
| `dust/v/rad` | (Nr, Nm) | cm/s | Drift velocity |
| `dust/v/frag` | (Nr,) | cm/s | Fragmentation velocity |
| `dust/p/stick`, `dust/p/frag` | (Nr, Nm, Nm) | — | Probability matrices |
| `dust/kernel` | (Nr, Nm, Nm) | cm²/s | Collision kernel (unless `save=False`) |
| `dust/Fi/adv`, `dust/Fi/diff`, `dust/Fi/tot` | (Nr+1, Nm) | g/cm/s | Fluxes at interfaces |

Only fields with `save != False` are written; custom fields/groups added via
`addgroup`/`addfield` are stored under their own path (e.g. `/planets/jupiter/M`).

---

## 3. Reading

### 3a. Direct (h5py) — for a single file or one field

```python
import h5py
import dustpy.constants as c

with h5py.File("data/data0017.hdf5", "r") as f:
    r     = f["grid/r"][:]            # (Nr,) cm
    m     = f["grid/m"][:]            # (Nm,) g
    t     = float(f["t"][()])         # s
    Sig_g = f["gas/Sigma"][:]         # (Nr,) g/cm^2
    Sig_d = f["dust/Sigma"][:]        # (Nr, Nm) g/cm^2
    St    = f["dust/St"][:]
    a     = f["dust/a"][:]
    eps   = f["dust/eps"][:]
```

List all snapshots first:

```python
import glob
snaps = sorted(glob.glob("data/data*.hdf5"))   # data0000.hdf5, ...
```

### 3b. simframe writer — for frame-like access

```python
from dustpy import hdf5writer

wrtr = hdf5writer(datadir="data")
one  = wrtr.read.output(0)        # snapshot 0, attribute access: one.gas.Sigma
all  = wrtr.read.all()            # everything; adds a leading time axis
seq  = wrtr.read.sequence("gas.Sigma")   # one field, shape (Nt, Nr)
```

### 3c. `dustpy.utils.read_data` — lightweight, plot-ready

```python
from dustpy.utils import read_data

data = read_data("data")     # or read_data(sim) with a live Simulation
# time:      data.t            (Nt,) s ; data.Nt
# grid:      data.grid.r (Nt, Nr) cm, data.grid.m (Nt, Nm) g, data.grid.ri,
#            data.grid.Nr, data.grid.Nm, data.grid.OmegaK
# gas:       data.gas.Sigma (Nt, Nr), data.gas.cs (Nt, Nr), data.gas.M (Nt,) (total mass, M_sun-free)
# dust:      data.dust.Sigma (Nt, Nr, Nm) g/cm^2 per bin
#            data.dust.sigma (Nt, Nr, Nm) bin-normalised density (≈ Sigma / (Δm/m))
#            data.dust.St (Nt, Nr, Nm), data.dust.eps (Nt, Nr), data.dust.M (Nt,)
#            data.dust.St_limits.frag (Nt, Nr), .drift (Nt, Nr) — growth-limit Stokes numbers
#            data.dust.v.frag (Nt, Nr)
```

This is what `dustpy.plot.panel` uses. Note the difference between
`dust.Sigma` (per mass bin) and `dust.sigma` (bin-normalised); use the
former for surface densities, the latter for density maps.

### 3d. Resuming a run

`frame.dmp` stores the complete frame (values + updaters + writer config):

```python
from dustpy import readdump
import dustpy.constants as c
import numpy as np

sim = readdump("data/frame.dmp")          # full Simulation, identical to the original
print(sim.t / c.year)                      # current time

# extend by snapshots and continue
sim.t.snapshots = np.hstack((sim.t.snapshots,
                             [1.e5, 5.e5, 1.e6] * c.year))
sim.run()
```

For a script-driven run, the equivalent "fresh Simulation" path works too:
build the same setup, `sim.initialize()`, then `sim.load("data/frame.dmp")` —
but `readdump` is the documented, safest route because it restores the exact
updater wiring of the original run.

---

## 4. Derived quantities (common)

```python
# Total dust surface density per radius [g/cm^2]
Sigma_dust_tot = Sig_d.sum(axis=1)

# Total disk masses [g]
A = np.pi * (r[1:]**2 if r.ndim == 1 else np.pi * (ri[1:]**2 - ri[:-1]**2))
M_gas = (A * Sig_g).sum()
M_dust = (A * Sig_dust_tot).sum()

# Peak particle size per radius (at the bin with max Sigma) [cm]
i_peak = np.argmax(Sig_d, axis=1)
a_max = a[np.arange(len(r)), i_peak]

# Largest Stokes among bins with significant dust
St_max = np.where(Sig_d > 1e-20, St, 0.).max(axis=1)

# Fragmentation/drift-limited Stokes numbers: read_data(...).dust.St_limits
```

To convert `Sigma_d[i, j]` (per bin) into a density per logarithmic mass
interval, divide by `B = 2 (a−1)/(a+1)` with `a = mean(m[1:]/m[:-1])` — this is
what `read_data` already does in `dust.sigma`.

---

## 5. Reporting conventions

When summarising results for the user:

- Report times in **years** (`t / c.year`); raw data are seconds.
- Report radii in **AU** (`r / c.au`); masses in **g** or **M_earth**
  (`M_earth ≈ 5.97e+27 g`).
- State the parameter provenance explicitly (from the setup script or read
  from the HDF5: grid extent, Nm, star mass, alpha, Mdisk, d2g, vFrag).
- Distinguish:
  - `dust/rho` — midplane **volumetric** density per mass bin, g/cm³;
  - `dust/Sigma` — surface density per mass bin, g/cm²;
  - `dust/eps` — vertically integrated dust-to-gas ratio, dimensionless.
- Quote the citation with every published use: Stammler & Birnstiel (2022),
  ApJ 935, 35, doi:10.3847/1538-4357/ac7d58.
