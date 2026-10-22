#!/usr/bin/env python3
"""
plot_dustpy.py — Headless diagnostic plots for DustPy HDF5 output.

Produces three plot types (selectable via --plot):
  panel     : 6-panel overview (density map, mass spectrum, gas/dust mass
               evolution, radial profile, Sigma(r) + eps(r)) — mirrors
               dustpy.plot.panel but saves to file without a display
  radial    : radial profiles of Sigma_gas, Sigma_dust, eps, a_max, St_max
               at one or more snapshots
  evolution : space-time diagrams of a_max(r, t) and Sigma_dust(r, t)

Usage
-----
  python plot_dustpy.py --run_dir data/dust/my_run [options]

  # All three plot types:
  python plot_dustpy.py --run_dir runs/my_run --plot all --out plots/my_run

  # Radial profiles at snapshots 0 and -1 (last):
  python plot_dustpy.py --run_dir runs/my_run --plot radial --snaps 0 -1

  # 6-panel overview at snapshot 5 as PNG:
  python plot_dustpy.py --run_dir runs/my_run --plot panel --snap 5 --fmt png

--run_dir may be either a run directory containing a `data/` subdirectory
with the HDF5 files, or a data directory containing data*.hdf5 directly.

Output (stdout, last line):
  SUCCESS plot_dir=<path>  files=<comma-separated list>
  ERROR   message=<description>
"""

import argparse
import glob
import os
import sys
import traceback
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

matplotlib.rcParams.update(
    {
        "figure.dpi": 300,
        "font.size": 11,
        "axes.labelsize": 11,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
    }
)


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------


def resolve_datadir(run_dir: str) -> str:
    """Accept a run dir (with data/ subdir) or a data dir directly."""
    sub = os.path.join(run_dir, "data")
    if os.path.isdir(sub):
        return sub
    if glob.glob(os.path.join(run_dir, "data*.hdf5")):
        return run_dir
    raise FileNotFoundError(
        f"No data/ subdirectory with data*.hdf5 files found in {run_dir}")


def _hdf5_files(datadir: str):
    files = sorted(glob.glob(os.path.join(datadir, "data*.hdf5")))
    if not files:
        raise FileNotFoundError(f"No data*.hdf5 files in {datadir}")
    return files


def _read_snap(path: str) -> dict:
    """Read one DustPy HDF5 snapshot. Returns dict of CGS arrays."""
    import h5py
    with h5py.File(path, "r") as f:
        r_cm = f["grid/r"][:]            # (Nr,)  cm
        m_g = f["grid/m"][:]             # (Nm,)  g
        t_s = float(f["t"][()])          # s
        Sigma_gas = f["gas/Sigma"][:]    # (Nr,)  g/cm^2
        Sigma_dust = f["dust/Sigma"][:]  # (Nr, Nm) g/cm^2
        St = f["dust/St"][:]             # (Nr, Nm)
        a_cm = f["dust/a"][:]            # (Nr, Nm) cm
        eps = f["dust/eps"][:]           # (Nr,)
    Sigma_dust_tot = Sigma_dust.sum(axis=1)  # (Nr,)
    # a_max: size at the bin with peak surface density per radial cell
    a_max = a_cm[np.arange(len(r_cm)), np.argmax(Sigma_dust, axis=1)]
    return dict(
        r_au=r_cm / c_au(),
        m_g=m_g,
        t_yr=t_s / c_year(),
        Sigma_gas=Sigma_gas,
        Sigma_dust=Sigma_dust,
        Sigma_dust_tot=Sigma_dust_tot,
        St=St,
        a_cm=a_cm,
        a_max_cm=a_max,
        eps=eps,
    )


def c_au():
    from dustpy import constants as c
    return c.au


def c_year():
    from dustpy import constants as c
    return c.year


# ---------------------------------------------------------------------------
# Plot: 6-panel overview (mirrors dustpy.plot.panel, saved to file)
# ---------------------------------------------------------------------------


def plot_panel(datadir: str, snap_index: int, out_path: str) -> str:
    from dustpy.utils import read_data

    data = read_data(datadir)
    Nt = data.Nt
    it = snap_index % Nt

    sd_max = np.ceil(np.log10(data.dust.sigma.max()))
    Mmax = np.ceil(np.log10(data.gas.M.max() / c_msun())) + 1
    levels = np.linspace(sd_max - 6, sd_max, 7)

    width = 3.5
    fig, axes = plt.subplots(2, 3, figsize=(3 * width, 2 * width / 1.618))
    ax00, ax01, ax02 = axes[0]
    ax10, ax11, _ = axes[1]
    ax11r = ax11.twinx()
    axes[1, 2].set_visible(False)

    # Panel 0,0 — dust density map (r-m) with growth limits
    cf = ax00.contourf(
        data.grid.r[it] / c_au(),
        data.grid.m[it],
        np.log10(data.dust.sigma[it].T),
        levels=levels,
        cmap="magma",
        extend="both",
    )
    ax00.contour(
        data.grid.r[it] / c_au(),
        data.grid.m[it],
        data.dust.St[it].T,
        levels=[1.0],
        colors="white",
        linewidths=1.5,
    )
    ax00.contour(
        data.grid.r[it] / c_au(),
        data.grid.m[it],
        (data.dust.St - data.dust.St_limits.drift[..., None])[it].T,
        levels=[0.0],
        colors="C2",
        linewidths=1,
    )
    ax00.contour(
        data.grid.r[it] / c_au(),
        data.grid.m[it],
        (data.dust.St - data.dust.St_limits.frag[..., None])[it].T,
        levels=[0.0],
        colors="C0",
        linewidths=1,
    )
    cbar = fig.colorbar(cf, ax=ax00)
    cbar.ax.set_ylabel(r"$\log_{10}\,\sigma_\mathrm{d}$ [g cm$^{-2}$]")
    ax00.set_xscale("log")
    ax00.set_yscale("log")
    ax00.set_xlabel("Distance [AU]")
    ax00.set_ylabel("Particle mass [g]")
    ax00.set_title(f"t = {data.t[it] / c_year():.2e} yr")

    # Panel 0,1 — mass spectrum at mid-disk radius
    ir_mid = int(data.grid.Nr) // 2
    ax01.loglog(data.grid.m[it], data.dust.sigma[it, ir_mid], c="C3")
    ax01.set_xlabel("Particle mass [g]")
    ax01.set_ylabel(r"$\sigma_\mathrm{d}$ [g cm$^{-2}$]")
    ax01.set_title(f"r = {data.grid.r[it, ir_mid] / c_au():.1f} AU")

    # Panel 0,2 — gas and dust mass vs time
    if data.Nt >= 3:
        ax02.loglog(data.t / c_year(), data.gas.M / c_msun(), label="Gas")
        ax02.loglog(data.t / c_year(), data.dust.M / c_msun(), label="Dust")
        ax02.axvline(data.t[it] / c_year(), color="#AAAAAA", lw=1, ls="--")
        ax02.set_xlim(data.t[1] / c_year(), data.t[-1] / c_year())
        ax02.set_ylim(10 ** (Mmax - 6), 10 ** Mmax)
        ax02.legend()
    ax02.set_xlabel("Time [yr]")
    ax02.set_ylabel(r"Mass [$M_\odot$]")

    # Panel 1,0 — radial dust profile at a mid-range mass bin
    im_mid = int(data.grid.Nm) // 2
    ax10.loglog(data.grid.r[it] / c_au(), data.dust.sigma[it, :, im_mid], c="C3")
    ax10.set_xlabel("Distance [AU]")
    ax10.set_ylabel(r"$\sigma_\mathrm{d}$ [g cm$^{-2}$]")
    ax10.set_title(f"m = {data.grid.m[it, im_mid]:.2e} g")

    # Panel 1,1 — Sigma and dust-to-gas ratio
    ax11.loglog(data.grid.r[it] / c_au(), data.gas.Sigma[it], label="Gas")
    ax11.loglog(data.grid.r[it] / c_au(), data.dust.Sigma[it].sum(-1), label="Dust")
    ax11.set_xlabel("Distance [AU]")
    ax11.set_ylabel(r"$\Sigma$ [g cm$^{-2}$]")
    ax11.legend()
    ax11r.loglog(data.grid.r[it] / c_au(), data.dust.eps[it], color="C7", lw=1, ls="--")
    ax11r.set_ylim(1e-5, 1e1)
    ax11r.set_ylabel("Dust-to-gas ratio")

    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


def c_msun():
    from dustpy import constants as c
    return c.M_sun


# ---------------------------------------------------------------------------
# Plot: radial profiles
# ---------------------------------------------------------------------------


def plot_radial(datadir: str, snap_indices, out_path: str) -> str:
    """Radial profiles of Sigma_gas, Sigma_dust, eps, a_max, St_max."""
    files = _hdf5_files(datadir)
    n = len(files)
    snaps = [i % n for i in snap_indices]

    fig, axes = plt.subplots(2, 2, figsize=(9, 7))
    ax_sig, ax_eps, ax_amax, ax_st = axes.flat

    cmap = plt.get_cmap("viridis")
    colors = [cmap(i / max(len(snaps) - 1, 1)) for i in range(len(snaps))]

    for color, idx in zip(colors, snaps):
        s = _read_snap(files[idx])
        lbl = f"t = {s['t_yr']:.2e} yr"
        r = s["r_au"]

        ax_sig.loglog(r, s["Sigma_gas"], color=color, ls="-", label=lbl)
        ax_sig.loglog(r, s["Sigma_dust_tot"], color=color, ls="--")

        ax_eps.loglog(r, s["eps"], color=color, label=lbl)

        ax_amax.loglog(r, s["a_max_cm"] * 10, color=color, label=lbl)  # cm -> mm

        # Max Stokes among bins with significant dust
        sig_mask = s["Sigma_dust"] > 1e-20
        St_max_r = np.where(sig_mask, s["St"], 0).max(axis=1)
        ax_st.loglog(r, St_max_r, color=color, label=lbl)

    ax_sig.set_xlabel("Distance [AU]")
    ax_sig.set_ylabel(r"$\Sigma$ [g cm$^{-2}$]")
    ax_sig.set_title("Surface density (solid=gas, dashed=dust)")
    ax_sig.legend(fontsize=7)

    ax_eps.set_xlabel("Distance [AU]")
    ax_eps.set_ylabel("Dust-to-gas ratio")
    ax_eps.axhline(0.5, color="r", ls=":", lw=1, label="SI threshold")
    ax_eps.legend(fontsize=7)

    ax_amax.set_xlabel("Distance [AU]")
    ax_amax.set_ylabel(r"$a_\mathrm{max}$ [mm]")
    ax_amax.set_title(r"Maximum grain size (at peak $\Sigma_\mathrm{dust}$)")
    ax_amax.legend(fontsize=7)

    ax_st.set_xlabel("Distance [AU]")
    ax_st.set_ylabel(r"$\mathrm{St}_\mathrm{max}$")
    ax_st.axhline(1.0, color="r", ls=":", lw=1, label="St = 1")
    ax_st.set_title("Max Stokes number (significant bins only)")
    ax_st.legend(fontsize=7)

    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


# ---------------------------------------------------------------------------
# Plot: space-time evolution
# ---------------------------------------------------------------------------


def plot_evolution(datadir: str, out_path: str) -> str:
    """Space-time diagrams of a_max(r,t) and Sigma_dust(r,t)."""
    files = _hdf5_files(datadir)

    snaps = [_read_snap(f) for f in files]
    r_au = snaps[-1]["r_au"]
    t_yr = np.array([s["t_yr"] for s in snaps])

    a_max_grid = np.array([s["a_max_cm"] * 10 for s in snaps])   # mm, (Nt, Nr)
    Sigma_d_grid = np.array([s["Sigma_dust_tot"] for s in snaps])  # (Nt, Nr)

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(12, 5))

    def _pcolormesh_log(ax, x, y, z, cmap, label):
        z_safe = np.where(z > 0, z, np.nan)
        pm = ax.pcolormesh(x, y, np.log10(z_safe), cmap=cmap, shading="auto")
        cb = fig.colorbar(pm, ax=ax)
        cb.set_label(label)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Distance [AU]")
        ax.set_ylabel("Time [yr]")

    _pcolormesh_log(ax0, r_au, t_yr, a_max_grid, "magma",
                    r"$\log_{10}\,a_\mathrm{max}$ [mm]")
    ax0.set_title("Maximum grain size evolution")

    _pcolormesh_log(ax1, r_au, t_yr, Sigma_d_grid, "viridis",
                    r"$\log_{10}\,\Sigma_\mathrm{dust}$ [g cm$^{-2}$]")
    ax1.set_title("Dust surface density evolution")

    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--run_dir", required=True,
                   help="Run directory (with data/ sub-directory) or data "
                        "directory containing data*.hdf5 files")
    p.add_argument("--plot", default="all",
                   choices=["panel", "radial", "evolution", "all"],
                   help="Which plot(s) to produce (default: all)")
    p.add_argument("--snaps", type=int, nargs="+", default=[0, -1],
                   help="Snapshot indices for --plot radial (default: 0 -1; negative = from end)")
    p.add_argument("--snap", type=int, default=-1,
                   help="Snapshot index for --plot panel (default: -1 = last)")
    p.add_argument("--out", default=None,
                   help="Output directory for plots (default: plots/<run_name>/)")
    p.add_argument("--fmt", default="pdf", choices=["pdf", "png", "both"],
                   help="Output format (default: pdf)")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    run_dir = os.path.abspath(args.run_dir)
    try:
        datadir = resolve_datadir(run_dir)
    except FileNotFoundError as e:
        print(f"ERROR message={e}", flush=True)
        sys.exit(1)

    run_name = Path(datadir).name
    out_dir = args.out or os.path.join("plots", run_name)
    os.makedirs(out_dir, exist_ok=True)

    fmts = ["pdf", "png"] if args.fmt == "both" else [args.fmt]
    made = []

    try:
        plots_requested = (["panel", "radial", "evolution"] if args.plot == "all"
                            else [args.plot])

        for fmt in fmts:
            if "panel" in plots_requested:
                out = os.path.join(out_dir, f"{run_name}_panel.{fmt}")
                plot_panel(datadir, args.snap, out)
                made.append(out)

            if "radial" in plots_requested:
                out = os.path.join(out_dir, f"{run_name}_radial.{fmt}")
                plot_radial(datadir, args.snaps, out)
                made.append(out)

            if "evolution" in plots_requested:
                out = os.path.join(out_dir, f"{run_name}_evolution.{fmt}")
                plot_evolution(datadir, out)
                made.append(out)

        print(f"SUCCESS plot_dir={out_dir}  files={','.join(made)}", flush=True)

    except Exception:
        print(f"ERROR message={traceback.format_exc()}", flush=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
