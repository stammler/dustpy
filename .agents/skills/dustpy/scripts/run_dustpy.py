#!/usr/bin/env python3
"""
run_dustpy.py — Launch script for DustPy dust-evolution simulations.

Usage (local):
    python run_dustpy.py --run_dir runs/my_run [options]
    python run_dustpy.py --run_dir runs/my_run --json '{"alpha": 1e-3, ...}'

Usage (dry run — validates parameters, writes setup script only):
    python run_dustpy.py --run_dir runs/my_run --dry_run

Usage (resume a run from its frame.dmp dump, extending it to t_end_yr):
    python run_dustpy.py --run_dir runs/my_run --resume --t_end_yr 1e6

Usage (HPC / SLURM):
    python run_dustpy.py --run_dir runs/my_run --hpc --dry_run   # write submit.sh
    python run_dustpy.py --run_dir runs/my_run --hpc             # write + sbatch

Output (stdout, last line):
    SUCCESS run_dir=<path> n_snapshots=<N> t_end_yr=<float> wall_clock_s=<float>
    ERROR message=<description>

All physical parameters must be in CGS unless the argument name ends in a
unit suffix (_msun, _rsun, _au, _yr, _cms) — those are converted internally.
"""

import argparse
import glob
import json
import os
import sys
import time
import traceback

import numpy as np

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------


def parse_args():
    p = argparse.ArgumentParser(description="Launch a DustPy dust-evolution simulation.")

    # --- run control ---
    p.add_argument("--run_dir", required=True,
                   help="Output directory for this run (created if absent); HDF5 output goes to <run_dir>/data.")
    p.add_argument("--json", default=None,
                   help="JSON string of parameter overrides (see known keys in the source).")
    p.add_argument("--overwrite", action="store_true",
                   help="Overwrite existing data directory. Destructive — use deliberately.")
    p.add_argument("--resume", action="store_true",
                   help="Resume from <run_dir>/data/frame.dmp and extend to t_end_yr.")
    p.add_argument("--dry_run", action="store_true",
                   help="Validate parameters and write the setup script/manifest only; do not run.")

    # --- HPC / SLURM ---
    p.add_argument("--hpc", action="store_true",
                   help="Write a SLURM jobscript (submit.sh) instead of running locally.")
    p.add_argument("--partition", default="serial", help="SLURM partition name.")
    p.add_argument("--walltime", default="04:00:00", help="SLURM wall time (HH:MM:SS).")
    p.add_argument("--memory_gb", type=float, default=8.0, help="SLURM memory per node in GB.")

    # --- stellar (physical units converted internally) ---
    p.add_argument("--star_mass_msun", type=float, default=1.0)
    p.add_argument("--star_radius_rsun", type=float, default=2.0)
    p.add_argument("--star_teff", type=float, default=5772.0)

    # --- grid ---
    p.add_argument("--nr", type=int, default=100, help="Radial grid cells.")
    p.add_argument("--rmin_au", type=float, default=1.0)
    p.add_argument("--rmax_au", type=float, default=1000.0)
    p.add_argument("--nmbpd", type=int, default=7, help="Mass bins per decade; minimum 7.")
    p.add_argument("--mmin_g", type=float, default=1e-12)
    p.add_argument("--mmax_g", type=float, default=1e5)

    # --- gas ---
    p.add_argument("--alpha", type=float, default=1e-3)
    p.add_argument("--mdisk_msun", type=float, default=0.05)
    p.add_argument("--sigma_exp", type=float, default=-1.0)
    p.add_argument("--sigma_rc_au", type=float, default=60.0)

    # --- dust ---
    p.add_argument("--d2g", type=float, default=0.01)
    p.add_argument("--vfrag_cms", type=float, default=100.0,
                   help="Fragmentation velocity in cm/s (default 100 = 1 m/s).")
    p.add_argument("--rho_monomer", type=float, default=1.67, help="Monomer bulk density [g/cm^3].")
    p.add_argument("--a_ini_max_cm", type=float, default=1e-4,
                   help="Max initial grain size [cm] (default 1 um).")

    # --- time ---
    p.add_argument("--t_start_yr", type=float, default=1e3)
    p.add_argument("--t_end_yr", type=float, default=1e5)
    p.add_argument("--snaps_per_decade", type=int, default=10)

    return p.parse_args()


# ---------------------------------------------------------------------------
# Parameter override via --json
# ---------------------------------------------------------------------------

KNOWN_JSON_KEYS = {
    "star_mass_msun", "star_radius_rsun", "star_teff",
    "nr", "rmin_au", "rmax_au", "nmbpd", "mmin_g", "mmax_g",
    "alpha", "mdisk_msun", "sigma_exp", "sigma_rc_au",
    "d2g", "vfrag_cms", "rho_monomer", "a_ini_max_cm",
    "t_start_yr", "t_end_yr", "snaps_per_decade",
}


def apply_json_overrides(args, json_str):
    try:
        overrides = json.loads(json_str)
    except json.JSONDecodeError as e:
        print(f"ERROR message=Invalid --json: {e}", flush=True)
        sys.exit(1)
    unknown = set(overrides) - KNOWN_JSON_KEYS
    if unknown:
        print(f"ERROR message=Unknown JSON keys: {sorted(unknown)}. "
              f"Known keys: {sorted(KNOWN_JSON_KEYS)}", flush=True)
        sys.exit(1)
    for key, val in overrides.items():
        setattr(args, key, val)
    return args


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def validate_params(args):
    errors = []
    if args.nmbpd < 7:
        errors.append(f"nmbpd={args.nmbpd} < 7 — required for accurate "
                      f"coagulation (Drążkowska et al. 2014)")
    if args.rmin_au <= 0 or args.rmax_au <= args.rmin_au:
        errors.append(f"Invalid grid: rmin={args.rmin_au} AU, rmax={args.rmax_au} AU")
    if args.t_end_yr <= args.t_start_yr:
        errors.append(f"t_end_yr={args.t_end_yr} must be > t_start_yr={args.t_start_yr}")
    if args.mmin_g <= 0 or args.mmax_g <= args.mmin_g:
        errors.append(f"Invalid mass grid: mmin={args.mmin_g} g, mmax={args.mmax_g} g")
    if args.alpha <= 0 or args.alpha > 0.1:
        errors.append(f"alpha={args.alpha} out of physical range (0, 0.1]")
    if args.vfrag_cms <= 0:
        errors.append(f"vfrag_cms={args.vfrag_cms} must be positive")
    if args.nr < 5:
        errors.append(f"nr={args.nr} < 5 — too few radial cells")
    if errors:
        for e in errors:
            print(f"ERROR message={e}", flush=True)
        sys.exit(1)


# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------

PARAM_KEYS = [
    "star_mass_msun", "star_radius_rsun", "star_teff",
    "nr", "rmin_au", "rmax_au", "nmbpd", "mmin_g", "mmax_g",
    "alpha", "mdisk_msun", "sigma_exp", "sigma_rc_au",
    "d2g", "vfrag_cms", "rho_monomer", "a_ini_max_cm",
    "t_start_yr", "t_end_yr", "snaps_per_decade",
]


def write_manifest(run_dir, datadir, args, mode):
    manifest = {
        "code": "DustPy",
        "mode": mode,
        "run_dir": run_dir,
        "data_dir": datadir,
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "parameters": {k: getattr(args, k) for k in PARAM_KEYS},
    }
    with open(os.path.join(run_dir, "run_manifest.json"), "w") as mf:
        json.dump(manifest, mf, indent=2)


# ---------------------------------------------------------------------------
# SLURM jobscript writer
# ---------------------------------------------------------------------------

SLURM_TEMPLATE = """\
#!/bin/bash
#SBATCH --job-name=dustpy_{run_name}
#SBATCH --partition={partition}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem={memory_mb}M
#SBATCH --time={walltime}
#SBATCH --output={run_dir}/slurm_%j.out
#SBATCH --error={run_dir}/slurm_%j.err

# Activate your environment if needed:
# source /path/to/venv/bin/activate

python {script_path} \\
    --run_dir {run_dir} \\
    --json '{params_json}'
"""


def write_slurm_script(args, params_json, run_dir):
    script_path = os.path.abspath(__file__)
    run_name = os.path.basename(run_dir.rstrip("/"))
    slurm_path = os.path.join(run_dir, "submit.sh")
    os.makedirs(run_dir, exist_ok=True)
    content = SLURM_TEMPLATE.format(
        run_name=run_name,
        partition=args.partition,
        memory_mb=int(args.memory_gb * 1024),
        walltime=args.walltime,
        run_dir=os.path.abspath(run_dir),
        script_path=script_path,
        params_json=params_json.replace("'", '"'),
    )
    with open(slurm_path, "w") as f:
        f.write(content)
    return slurm_path


# ---------------------------------------------------------------------------
# Snapshot helper
# ---------------------------------------------------------------------------


def snapshot_array(t_start_yr, t_end_yr, snaps_per_decade):
    """Log-spaced snapshot times in seconds (geometric, inclusive)."""
    n_dec = np.log10(t_end_yr / t_start_yr)
    num = max(int(snaps_per_decade * n_dec) + 1, 2)
    return np.geomspace(t_start_yr, t_end_yr, num=num) * c_year()


def c_year():
    from dustpy import constants as c
    return c.year


# ---------------------------------------------------------------------------
# Simulation setup and run
# ---------------------------------------------------------------------------


def build_and_run(args):
    from dustpy import Simulation
    from dustpy import constants as c

    run_dir = os.path.abspath(args.run_dir)
    datadir = os.path.join(run_dir, "data")
    os.makedirs(run_dir, exist_ok=True)

    if args.resume:
        # ------------------------------------------------------------------
        # Resume: restore the full frame from the dump
        # ------------------------------------------------------------------
        from dustpy import readdump
        dump_path = os.path.join(datadir, "frame.dmp")
        if not os.path.exists(dump_path):
            print(f"ERROR message=Resume requested but dump file not found: {dump_path}",
                  flush=True)
            sys.exit(1)
        sim = readdump(dump_path)
        t_now = float(sim.t)
        write_manifest(run_dir, datadir, args, "resume")
        if args.dry_run:
            print(f"SUCCESS run_dir={run_dir} dry_run=True resume=True "
                  f"t_now_yr={t_now / c.year:.4g}", flush=True)
            return

        # Extend snapshots if t_end is beyond the current time
        if args.t_end_yr * c.year > t_now:
            t0 = max(t_now, args.t_start_yr * c.year)
            new_snaps = snapshot_array(t0 / c.year, args.t_end_yr, args.snaps_per_decade)
            sim.t.snapshots = np.unique(
                np.concatenate([np.asarray(sim.t.snapshots), new_snaps]))
        sim.writer.datadir = datadir
        sim.writer.overwrite = args.overwrite
    else:
        # ------------------------------------------------------------------
        # Fresh run
        # ------------------------------------------------------------------
        sim = Simulation()

        # stellar
        sim.ini.star.M = args.star_mass_msun * c.M_sun
        sim.ini.star.R = args.star_radius_rsun * c.R_sun
        sim.ini.star.T = args.star_teff

        # grid
        sim.ini.grid.Nr = args.nr
        sim.ini.grid.rmin = args.rmin_au * c.au
        sim.ini.grid.rmax = args.rmax_au * c.au
        sim.ini.grid.Nmbpd = args.nmbpd
        sim.ini.grid.mmin = args.mmin_g
        sim.ini.grid.mmax = args.mmax_g

        # gas
        sim.ini.gas.alpha = args.alpha
        sim.ini.gas.Mdisk = args.mdisk_msun * c.M_sun
        sim.ini.gas.SigmaExp = args.sigma_exp
        sim.ini.gas.SigmaRc = args.sigma_rc_au * c.au

        # dust
        sim.ini.dust.d2gRatio = args.d2g
        sim.ini.dust.vFrag = args.vfrag_cms
        sim.ini.dust.rhoMonomer = args.rho_monomer
        sim.ini.dust.aIniMax = args.a_ini_max_cm

        write_manifest(run_dir, datadir, args, "fresh")
        if args.dry_run:
            n_snaps_dry = len(snapshot_array(args.t_start_yr, args.t_end_yr,
                                             args.snaps_per_decade)) + 1
            print(f"SUCCESS run_dir={run_dir} dry_run=True "
                  f"n_snapshots={n_snaps_dry} t_end_yr={args.t_end_yr:.4g}", flush=True)
            return

        sim.initialize()

        # writer and snapshots — set after initialize()
        sim.writer.datadir = datadir
        sim.writer.overwrite = args.overwrite
        sim.t.snapshots = np.hstack(
            [sim.t,
             snapshot_array(args.t_start_yr, args.t_end_yr, args.snaps_per_decade)])

    # ------------------------------------------------------------------
    # Run
    # ------------------------------------------------------------------
    t0 = time.perf_counter()
    sim.run()
    wall_clock_s = time.perf_counter() - t0

    hdf5_files = sorted(glob.glob(os.path.join(datadir, "data*.hdf5")))
    n_snaps = len(hdf5_files)
    if n_snaps == 0:
        print(f"ERROR message=Run completed but no HDF5 files found in {datadir}",
              flush=True)
        sys.exit(1)

    # Verify last snapshot is readable
    import h5py as _h5
    with _h5.File(hdf5_files[-1], "r") as f:
        t_end_s = float(f["t"][()])
    t_end_yr = t_end_s / c.year

    print(f"SUCCESS run_dir={run_dir} n_snapshots={n_snaps} "
          f"t_end_yr={t_end_yr:.4g} wall_clock_s={wall_clock_s:.1f}", flush=True)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main():
    args = parse_args()

    if args.json:
        args = apply_json_overrides(args, args.json)

    validate_params(args)

    run_dir = os.path.abspath(args.run_dir)
    os.makedirs(run_dir, exist_ok=True)
    params_json = json.dumps({k: getattr(args, k) for k in PARAM_KEYS})

    if args.hpc:
        if args.resume:
            print("ERROR message=--resume is not supported together with --hpc "
                  "(SLURM jobscripts are for fresh runs)", flush=True)
            sys.exit(1)
        slurm_path = write_slurm_script(args, params_json, args.run_dir)
        if args.dry_run:
            print(f"SUCCESS run_dir={run_dir} slurm_script={slurm_path} dry_run=True",
                  flush=True)
        else:
            import subprocess
            result = subprocess.run(["sbatch", slurm_path],
                                    capture_output=True, text=True)
            if result.returncode != 0:
                print(f"ERROR message=sbatch failed: {result.stderr.strip()}", flush=True)
                sys.exit(1)
            job_id = result.stdout.strip().split()[-1]
            print(f"SUCCESS run_dir={run_dir} slurm_job_id={job_id} slurm_script={slurm_path}",
                  flush=True)
        return

    try:
        build_and_run(args)
    except SystemExit:
        raise
    except Exception:
        tb = traceback.format_exc()
        print(f"ERROR message={tb[-1800:]}", flush=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
