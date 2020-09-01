import dustpy.constants as c
from dustpy.std import dust_f

import numpy as np
from scipy.interpolate import interp1d


def boundary(sim):
    """Function sets the boundary condition of dust surface density.
    Not implemented, yet.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame"""
    pass


def enforce_floor_value(sim):
    """Function enforces floor value onto dust surface density.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame"""
    sim.dust.Sigma = np.maximum(sim.dust.Sigma, sim.dust.SigmaFloor)


def finalize(sim):
    """Function finalizes integration step.

    Parameters
    ----------
    sim : Frame
        Parent integration frame"""
    boundary(sim)
    enforce_floor_value(sim)


def a(sim):
    """Function calculates the particle size from the solid density and the filling factor.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    a : Field
        Particle sizes"""
    rho = sim.dust.fill * sim.dust.rhos
    return dust_f.a(sim.grid.m, rho)


def D(sim):
    """Function calculates the dust diffusivity.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    D : Field
        Dust diffusivity"""
    v2 = sim.dust.delta.rad * sim.gas.cs**2
    return dust_f.d(v2, sim.grid.OmegaK, sim.dust.St)


def eps(sim):
    """Function returns the vertically integrated dust-to-gas ratio.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    eps : Field
        vertically integrated dust-to-gas ratio"""
    return np.sum(sim.dust.Sigma, axis=-1) / sim.gas.Sigma


def F_adv(sim, Sigma=None):
    """Function calculates the advective flux at the cell interfaces. It is linearly interpolating
    the velocities onto the grid cel interfaces and is assuming
    vi(0, :) = vi(1, :) and vi(Nr, :) = vi(Nr-1, :).

    Parameters
    ----------
    sim : Frame
        Parent simulation frame
    Sigma : Field, optional, default : None
        Surface density to be used if not None

    Returns
    -------
    Fi : Field
        Advective mass fluxes through the grid cell interfaces"""
    Sigma = Sigma if Sigma is not None else sim.dust.Sigma
    return dust_f.fi_adv(sim.dust.Sigma, sim.dust.v.rad, sim.grid.r, sim.grid.ri)


def F_diff(sim, Sigma=None):
    '''Function calculates the diffusive flux at the cell interfaces'''
    if Sigma is None:
        Sigma = sim.dust.Sigma
    return dust_f.fi_diff(sim.dust.D,
                          sim.dust.eps,
                          sim.dust.H,
                          sim.gas.Hp,
                          sim.grid.OmegaK,
                          sim.dust.rho,
                          sim.gas.rho,
                          sim.grid.r,
                          sim.grid.ri,
                          Sigma,
                          sim.gas.Sigma,
                          sim.dust.St,
                          sim.dust.delta.rad*sim.gas.cs**2)


def f(St):
    """This is a helper function for the calculation of the diffusive fluxes."""
    return 2*St / (1 + 2*St**3)


def F_tot(sim, Sigma=None):
    """Function calculates the total mass fluxes through grid cell interfaces.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame
    Sigma : Field, optional, default : None
        Surface density to be used if not None

    Returns
    -------
    Ftot : Field
        Total mass flux through interfaces"""
    Fi = np.zeros((np.int(sim.grid.Nr+1), np.int(sim.grid.Nm)))
    if Sigma is None:
        Sigma = sim.dust.Sigma
        Fdiff = sim.dust.Fi.diff
        Fadv = sim.dust.Fi.adv
    else:
        Fdiff = sim.dust.Fi.diff.updater.beat(sim, Sigma=Sigma)
        Fadv = sim.dust.Fi.adv.updater.beat(sim, Sigma=Sigma)
    if Fdiff is not None:
        Fi += Fdiff
    if Fadv is not None:
        Fi += Fadv
    return Fi


def H(sim):
    """Function calculates the dust scale height according Dubrulle et al. (1995).

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    H : Field
        Dust scale heights"""
    return dust_f.h_dubrulle1995(sim.gas.Hp, sim.dust.St, sim.dust.delta.vert)


def MRN_distribution(sim):
    """Function calculates the initial particle mass distribution. The parameters are taken from the
    ``Simulation.ini`` object.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    Sigma : Field
        Initial dust surface density

    Notes
    -----
    ``sim.ini.dust.aIniMax`` : maximum initial particle size
    ``sim.ini.dust.d2gRatio`` : initial dust-to-gas ratio
    ``sim.ini.dust.distExp`` : initial particle size distribution ``n(a) da ∝ a^{distExp} da``
    ``sim.ini.dust.allowDriftingParticles`` : if ``True`` the particle size distribution
    will be filled up to ``aIniMax``, if ``False`` the maximum particle size will be chosen
    such that there are no drifting particles initially. This prevents a particle wave traveling
    though the simulation, that is already drifting initially."""
    exp = sim.ini.dust.distExp
    # Calculating pressure gradient
    P = sim.gas.P
    fP = interp1d(sim.grid.r, P, kind="linear", fill_value="extrapolate")
    Pi = fP(sim.grid.ri)
    gamma = (Pi[1:] - Pi[:-1]) / (sim.grid.ri[1:] - sim.grid.ri[:-1])
    gamma = np.abs(gamma)
    # Exponent of pressure gradient
    gamma *= sim.grid.r / P
    gamma = 1. / gamma
    # Maximum drift limited particle size with safety margin
    ad = 1.e-4 * 2./np.pi * sim.ini.dust.d2gRatio * sim.gas.Sigma \
        / sim.dust.fill[:, 0] * sim.dust.rhos[:, 0] * (sim.grid.OmegaK * sim.grid.r)**2. \
        / sim.gas.cs**2. / gamma
    # Set maximum particle size
    if(sim.ini.dust.allowDriftingParticles):
        aIni = sim.ini.dust.aIniMax
    else:
        aIni = np.minimum(sim.ini.dust.aIniMax, ad)[:, None]
    # Fill distribution
    ret = np.where(sim.dust.a <= aIni, sim.dust.a**(exp+4), 0.)
    s = np.sum(ret, axis=1)[..., None]
    s = np.where(s > 0., s, 1.)
    # Normalize to mass
    ret = ret / s * sim.gas.Sigma[..., None] * sim.ini.dust.d2gRatio
    return ret


def rho_midplane(sim):
    """Function calculates the midplane mass density.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    rho : Field
        Midplane mass density"""
    return sim.dust.Sigma / (np.sqrt(2 * c.pi) * sim.dust.H)


def S_coag(sim, Sigma=None):
    """Function calculates the coagulation source terms.

    Parameters
    ----------
    sim : Frame
        Parent simulation Frame
    Sigma : Field, optional, default : None
        Surface density to be used if not None

    Returns
    -------
    Scoag : Field
        Coagulation source terms"""
    Sigma = Sigma if Sigma is not None else sim.dust.Sigma
    return sim.dust.S.coag


def S_hyd(sim, Sigma=None):
    """Function calculates the hydrodynamic source terms.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame
    Sigma : Field, optional, default : None
        Surface density to be used if not None

    Returns
    -------
    Shyd : Field
        Hydrodynamic source terms"""
    if Sigma is None:
        Sigma = sim.dust.Sigma
        Fi = sim.dust.Fi.tot
    else:
        Fi = sim.dust.Fi.tot.updater.beat(sim, Sigma=Sigma)
    if Fi is None:
        Fi = sim.dust.Fi.tot
    return dust_f.s_hyd(Fi, sim.grid.ri)


def S_tot(sim, Sigma=None):
    """Function calculates the total source terms.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame
    Sigma : Field, optional, default : None
        Surface density to be used if not None

    Returns
    -------
    Stot : Field
        Total source terms of surface density"""
    if Sigma is None:
        Sigma = sim.dust.Sigma
        Scoag = sim.dust.S.coag
        Shyd = sim.dust.S.hyd
    else:
        Scoag = sim.dust.S.coag.updater.beat(sim, Sigma=Sigma)
        Shyd = sim.dust.S.hyd.updater.beat(sim, Sigma=Sigma)
    S = np.zeros_like(sim.dust.S.tot)
    if Scoag is not None:
        S += Scoag
    if Shyd is not None:
        S += Shyd
    S += sim.dust.S.ext
    return S


def Sigma_deriv(sim, t, Sigma):
    """Function calculates the derivative of the dust surface density.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame
    t : IntVar
        Current time
    Sigma : Field
        Current Surface density

    Returns
    -------
    Sigma_dot: Field
        Derivative of Surface density"""
    return sim.dust.S.tot.updater.beat(sim, Sigma=Sigma)


def SigmaFloor(sim):
    """Function calculates the floor value for the dust distribution. Floor value means that there is less than
    one particle in an annulus.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    Sigma_floor : Field
        Floor value of surface density"""
    area = c.pi * (sim.grid.ri[1:]**2. - sim.grid.ri[:-1]**2.)
    return 1 / area[..., None] * sim.grid.m


def St_Epstein_StokesI(sim):
    """Function calculates the Stokes number using the Epstein and the Stokes I drag regimes.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    St : Field
        Stokes number"""
    rho = sim.dust.rhos * sim.dust.fill
    return dust_f.st_epstein_stokes1(sim.dust.a, sim.gas.mfp, rho, sim.gas.Sigma)


def vdriftmax(sim):
    """Function calculates the maximum drift velocity of the dust including back reaction of
    the gas onto the dust, if the back reaction coefficients are set.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    vdriftmax : Field
        Maximum drift velocity"""
    A = sim.dust.backreaction.A
    B = sim.dust.backreaction.B
    return 0.5 * B * sim.gas.v.visc - A * \
        sim.gas.eta * sim.grid.r * sim.grid.OmegaK


def vrad(sim):
    """Function calculated the radial velocity of the dust.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    vrad : Field
        Radial dust velocity"""
    return dust_f.vrad(sim.dust.St, sim.dust.v.driftmax, sim.gas.v.rad)


def vrel_azimuthal_drift(sim):
    """Function calculates the relative particle velocities due to azimuthal drift.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    vrel : Field
        Relative velocities"""
    return dust_f.vrel_azimuthal_drift(sim.dust.v.driftmax, sim.dust.St)


def vrel_brownian_motion(sim):
    """Function calculates the relative particle velocities due to Brownian motion.
    The maximum value is set to the sound speed.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    vrel : Field
        Relative velocities"""
    return dust_f.vrel_brownian_motion(sim.gas.cs, sim.grid.m, sim.gas.T)


def vrel_radial_drift(sim):
    """Function calculates the relative particle velocities due to radial drift.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    vrel : Field
        Relative velocities"""
    return dust_f.vrel_radial_drift(sim.dust.v.rad)


def vrel_tot(sim):
    """Function calculates the total relative vparticle velocities by taking the root mean square
    of all individual sources.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    vrel : Field
        Relative velocities"""
    ret = sim.dust.v.rel.azi**2.
    ret += sim.dust.v.rel.brown**2
    ret += sim.dust.v.rel.rad**2
    ret += sim.dust.v.rel.turb**2
    ret += sim.dust.v.rel.vert**2
    return np.sqrt(ret)


def vrel_turbulent_motion(sim):
    """Function calculates the relative particle velocities due to turbulent motion.
    It uses the prescription of Cuzzi & Ormel (2007).

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    vrel : Field
        Relative velocities"""
    return dust_f.vrel_cuzzi_ormel_2007(
        sim.dust.delta.turb,
        sim.gas.cs,
        sim.gas.mu,
        sim.grid.OmegaK,
        sim.gas.Sigma,
        sim.dust.St)


def vrel_vertical_settling(sim):
    """Function calculates the relative particle velocities due to vertical settling.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    vrel : Field
        Relative velocities"""
    return dust_f.vrel_vertical_settling(sim.dust.H, sim.grid.OmegaK, sim.dust.St)
