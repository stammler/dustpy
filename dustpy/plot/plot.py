import dustpy.constants as c
from dustpy.simulation import Simulation

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from simframe.io.writers import hdf5writer as writer
import os


def panel(data, im=0, ir=0, it=0, filename="data", extension="hdf5"):
    '''Simple plotting script for data files or simulation objects.

    Parameters
    ----------
    data : ``dustpy.Simulation`` or string
        Either instance of ``dustpy.Simulation`` or path to data directory to be plotted
    im : int, optional, default : 0
        Number of mass bin along which density distribution is plotted
    ir : int, optional, default : 0
        Number of radial grid index along density distribution is plotted
    it : int, optional, default : 0
        Index of snapshot to be plotted
    filename : string, optional, default : "data"
    extension : string, optional, default : "hdf5"
        Plotting script is looking for files with pattern ``<data>/<filename>*.<extension>``'''

    # Reading data
    single = False
    isframe = False

    # If input is simulation frame
    if isinstance(data, Simulation):
        Nm = data.grid.Nm
        Nr = data.grid.Nr
        m = data.grid.m
        r = data.grid.r
        ri = data.grid.ri
        t = data.t

        SigmaD = data.dust.Sigma
        SigmaG = data.gas.Sigma
        eps = data.dust.eps

        cs = data.gas.cs
        delta = data.dust.delta.turb
        OmegaK = data.grid.OmegaK
        St = data.dust.St
        vK = OmegaK * r
        vFrag = data.dust.v.frag

        single = True
        isframe = True
    # If input is data dir
    elif os.path.isdir(data):
        writer.filename = filename
        writer.extension = extension
        writer.datadir = data
        data = writer.read.all()

        Nm = data.grid.Nm[0, ...]
        Nr = data.grid.Nr[0, ...]
        m = data.grid.m[0, ...]
        r = data.grid.r[0, ...]
        ri = data.grid.ri[0, ...]
        t = data.t

        SigmaD = data.dust.Sigma
        SigmaG = data.gas.Sigma
        eps = data.dust.eps

        cs = data.gas.cs
        delta = data.dust.delta.turb
        OmegaK = data.grid.OmegaK
        St = data.dust.St
        vK = OmegaK * r
        vFrag = data.dust.v.frag

        Mgas = (np.pi * (ri[1:]**2 - ri[:-1]**2) * SigmaG[...]).sum(-1)
        Mdust = (np.pi * (ri[1:]**2 - ri[:-1]**2)
                 * SigmaD[...].sum(-1)).sum(-1)

        if len(data.t) < 3:
            single = True

    # Fix indices if necessary
    im = np.maximum(0, im)
    im = np.minimum(im, Nm-1)
    ir = np.maximum(0, ir)
    ir = np.minimum(ir, Nr-1)
    it = np.maximum(0, it)
    it = np.minimum(it, len(t)-1)

    # Transformation of distribution
    a = np.mean(m[..., 1:] / m[..., :-1], axis=-1)
    dm = 2. * (a - 1.) / (a + 1.)
    sigmaD = SigmaD[..., :] / dm

    # Calculating limits

    # Fragmentation limit
    b = vFrag**2 / (delta * cs**2)
    with np.warnings.catch_warnings():
        np.warnings.filterwarnings(
            'ignore',
            r'invalid value encountered in sqrt')
        St_fr = 1 / (2 * b) * (3 - np.sqrt(9 - 4 * b**2))

    # Drift limit
    p = SigmaG * OmegaK * cs

    _f = interp1d(np.log10(r), np.log10(p), fill_value='extrapolate')
    pi = 10.**_f(np.log10(ri))
    gamma = np.abs(r / p * np.diff(pi) / np.diff(ri))
    St_dr = eps / gamma * (vK / cs)**2

    # Get limits
    sd_max = np.ceil(np.log10(sigmaD.max()))
    sg_max = np.ceil(np.log10(SigmaG.max()))
    if not single:
        Mmax = np.ceil(np.log10(Mgas.max()/c.M_sun)) + 1

    width = 4.
    fig = plt.figure(figsize=(3.*width, 2.*width/1.618), dpi=150)
    axcmap = fig.add_subplot(231)
    axsigm = fig.add_subplot(232)
    axMass = fig.add_subplot(233)
    axsigr = fig.add_subplot(234)
    axSigma = fig.add_subplot(235)
    axd2g = axSigma.twinx()

    if isframe:
        pltcmap = axcmap.contourf(r/c.au,
                                  m,
                                  np.log10(sigmaD.T),
                                  levels=np.linspace(sd_max-6, sd_max, 7),
                                  cmap="magma",
                                  extend="both"
                                  )
        axcmap.contour(r/c.au,
                       m,
                       St.T,
                       levels=[1.],
                       colors="white",
                       linewidths=2
                       )
        axcmap.contour(r/c.au,
                       m,
                       (St - St_dr[..., None]).T,
                       levels=[0.],
                       colors="C2",
                       linewidths=1
                       )
        axcmap.contour(r/c.au,
                       m,
                       (St - St_fr[..., None]).T,
                       levels=[0.],
                       colors="C0",
                       linewidths=1
                       )
    else:
        pltcmap = axcmap.contourf(r/c.au,
                                  m,
                                  np.log10(sigmaD[it, ...].T),
                                  levels=np.linspace(sd_max-6, sd_max, 7),
                                  cmap="magma",
                                  extend="both"
                                  )
        axcmap.contour(r/c.au,
                       m,
                       St[it, ...].T,
                       levels=[1.],
                       colors="white",
                       linewidths=2
                       )
        axcmap.contour(r/c.au,
                       m,
                       (St - St_dr[..., None])[it, ...].T,
                       levels=[0.],
                       colors="C2",
                       linewidths=1
                       )
        axcmap.contour(r/c.au,
                       m,
                       (St - St_fr[..., None])[it, ...].T,
                       levels=[0.],
                       colors="C0",
                       linewidths=1
                       )

    axcmap.axhline(m[im], color="#AAAAAA", lw=1, ls="--")
    axcmap.axvline(r[ir]/c.au, color="#AAAAAA", lw=1, ls="--")
    cbarcmap = plt.colorbar(pltcmap, ax=axcmap)
    cbarcmap.ax.set_ylabel("$\log\ \sigma$ [g/cm²]")
    axcmap.set_xlim(r[0]/c.au, r[-1]/c.au)
    axcmap.set_ylim(m[0], m[-1])
    axcmap.set_xscale("log")
    axcmap.set_yscale("log")
    axcmap.set_xlabel("Distance from star [au]")
    axcmap.set_ylabel("Particle mass [g]")

    if isframe:
        axsigm.loglog(m, sigmaD[ir, :][0], c="C3")
    else:
        axsigm.loglog(m, sigmaD[it, ir, :], c="C3")
    axsigm.set_xlim(m[0], m[-1])
    axsigm.set_ylim(10.**(sd_max-6.), 10.**sd_max)
    axsigm.set_xlabel("Particle mass [g]")
    axsigm.set_ylabel("$\sigma$ [g/cm²]")

    axMass.set_xlabel("Time [yrs]")
    axMass.set_ylabel("Mass [$M_\odot$]")
    if single:
        axMass.set_xticks([0., 1.])
        axMass.set_yticks([0., 1.])
        axMass.text(0.5,
                    0.5,
                    "Not enough data points.",
                    verticalalignment="center",
                    horizontalalignment="center",
                    size="large")
    else:
        axMass.loglog(t/c.year, Mgas/c.M_sun, c="C0", label="Gas")
        axMass.loglog(t/c.year, Mdust/c.M_sun, c="C1", label="Dust")
        axMass.axvline(t[it]/c.year, c="gray", lw=1)
        axMass.set_xlim(t[1]/c.year, t[-1]/c.year)
        axMass.set_ylim(10.**(Mmax-6.), 10.**Mmax)
        axMass.legend()

    if isframe:
        axsigr.loglog(r/c.au, sigmaD[..., :, im], c="C3")
    else:
        axsigr.loglog(r/c.au, sigmaD[it, :, im], c="C3")
    axsigr.set_xlim(r[0]/c.au, r[-1]/c.au)
    axsigr.set_ylim(10.**(sd_max-6.), 10.**sd_max)
    axsigr.set_xlabel("Distance from star [au]")
    axsigr.set_ylabel("$\sigma$ [g/cm²]")

    if isframe:
        axSigma.loglog(r/c.au, SigmaG, label="Gas")
        axSigma.loglog(r/c.au, SigmaD.sum(-1), label="Dust")
    else:
        axSigma.loglog(r/c.au, SigmaG[it, ...], label="Gas")
        axSigma.loglog(r/c.au, SigmaD[it, ...].sum(-1), label="Dust")
    axSigma.set_xlim(r[0]/c.au, r[-1]/c.au)
    axSigma.set_ylim(10.**(sg_max-6), 10.**sg_max)
    axSigma.set_xlabel("Distance from star [au]")
    axSigma.set_ylabel("$\Sigma$ [g/cm²]")
    axSigma.legend()
    if isframe:
        axd2g.loglog(r/c.au, eps, color="C7", lw=1)
    else:
        axd2g.loglog(r/c.au, eps[it, ...], color="C7", lw=1)
    axd2g.set_ylim(1.e-5, 1.e1)
    axd2g.set_ylabel("Dust-to-gas ratio")

    fig.tight_layout()

    plt.show()
