import numpy as np
from scipy.interpolate import interp1d

import dustpy.constants as c
from dustpy.std import gas_f


def boundary(sim):
    """Function set the boundary conditions of the gas.
    Not implemented, yet.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame"""
    pass


def enforce_floor_value(sim):
    """Function enforces floor value to gas surface density.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame"""
    sim.gas.Sigma = np.maximum(sim.gas.Sigma, sim.gas.SigmaFloor)


def finalize(sim):
    """Function finalizes gas integration step.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame"""
    boundary(sim)
    enforce_floor_value(sim)


def cs_adiabatic(sim):
    """Function calculates the adiabatic sound speed of the gas. For the isothermal sound speed set
    ``Simulation.gas.gamma = 1``.

    Paramters
    ---------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    cs : Field
        Sound speed"""
    return np.sqrt(sim.gas.gamma * c.k_B * sim.gas.T / sim.gas.mu)


def eta_midplane(sim):
    """Function calculates the midplane pressure gradient parameter.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    eta : Field
        eta pressure gradient parameter"""
    return gas_f.eta_midplane(sim.gas.Hp, sim.gas.P, sim.grid.r, sim.grid.ri)


def Fi(sim, Sigma=None):
    """Function calculates the mass flux through the cell interfaces.
    The fluxes are calculated with donor cell assuming
    vrad_i(0) = vrad_i(1) and vrad_i(Nr) = vrad_i(Nr-1).

    Parameters
    ----------
    sim : Frame
        Parent simulation frame
    Sigma : Field, optional, default : None
        Surface density to be used if not None

    Returns
    -------
    Fi : Field
        Mass flux through grid cell interfaces"""
    if Sigma is None:
        vrad = sim.gas.v.rad
        Sigma = sim.gas.Sigma
    else:
        vrad = sim.gas.v.rad.updater.beat(sim, Sigma=Sigma)
    if vrad is None:
        vrad = sim.gas.v.rad
    return gas_f.fi(Sigma, vrad, sim.grid.r, sim.grid.ri)


def Hp(sim):
    """Function calculates the pressure scale height.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    Hp : Field
        Pressure scale height"""
    return sim.gas.cs/sim.grid.OmegaK


def lyndenbellpringle1974(r, rc, p, Mdisk):
    """Function calculates the surface density according the self similar solution of Lynden-Bell & Pringle (1974).

    Parameters
    ----------
    r : float or array of floats
        radial distance from star
    rc : float
        critical cutoff radius
    p : float
        power law exponent
    Mdist : float
        disk mass

    Returns
    -------
    Sigma : float or array of floats
        Surface density profile"""
    return Mdisk / (2.*np.pi*rc**2*(2+p)) * (r/rc)**p * np.exp(-(r/rc)**(2+p))


def mfp_midplane(sim):
    """Function calculates the midplane mean free path.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    mfp : Field
        Mean free path"""
    return 1. / (np.sqrt(2.) * sim.gas.n * c.sigma_H2)


def n_midplane(sim):
    """Function calculates the midplane number density of the gas.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    n : Field
        Midplane number density"""
    return sim.gas.rho / sim.gas.mu


def P_midplane(sim):
    """Function calculates the midplane gas pressure.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    P : Field
        Midplane pressure"""
    return sim.gas.rho * sim.gas.cs**2 / sim.gas.gamma


def rho_midplane(sim):
    """Function calculates the midplane mass density of the gas.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    rho : Field
        Midplane mass density"""
    return sim.gas.Sigma / (np.sqrt(2. * c.pi) * sim.gas.Hp)


def Sigma_deriv(sim, t, Sigma):
    """Function calculates the derivative of the gas surface density.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame
    t : IntVar
        Current time
    Sigma : Field
        Surface density

    Returns
    -------
    Sigma_dot : Field
        Derivative of surface density"""
    S = sim.gas.S.tot.updater.beat(sim, Sigma=Sigma)
    if S is not None:
        return S
    return np.zeros(np.int(sim.grid.Nr))


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
    S_hyd : Field
        Hydrodynamic source terms"""
    if Sigma is None:
        Fi = sim.gas.Fi
    else:
        Fi = sim.gas.Fi.updater.beat(sim, Sigma=Sigma)
    if Fi is None:
        Fi = sim.gas.Fi
    return gas_f.s_hyd(Fi, sim.grid.ri)


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
    S_tot : Field
        Total surface density source terms"""
    # If Sigma is None then the updater is called and not the derivative
    if Sigma is None:
        return sim.gas.S.ext + sim.gas.S.hyd
    Shyd = sim.gas.S.hyd.updater.beat(sim, Sigma=Sigma)
    if Shyd is None:
        Shyd = sim.gas.S.hyd
    return Shyd + sim.gas.S.ext


def T_passive(sim):
    """Function calculates the temperature profile of a passively irridiated disk with a constant irradiation
    angle of 0.05.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    T : Field
        Gas temperature"""
    return (sim.star.L * 0.05 / (4. * c.pi * sim.grid.r**2 * c.sigma_sb))**0.25


def vrad(sim, Sigma=None):
    """Function calculates the radial radial gas velocity.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame
    Sigma : Field, optional, default : None
        Surface density to be used if not None

    Returns
    -------
    vrad : Field
        Radial gas velocity"""
    if Sigma is None:
        return sim.gas.v.visc
    return sim.gas.v.visc.updater.beat(sim, Sigma=Sigma)


def vvisc(sim, Sigma=None):
    """Function calculates the viscous radial gas velocity.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame
    Sigma : Field, optional, default : None
        Surface density to be used if not None

    Returns
    -------
    vvisc : Field
        Viscous radial gas velocity"""
    if Sigma is None:
        Sigma = sim.gas.Sigma
    nu = sim.gas.alpha * sim.gas.cs * sim.gas.Hp
    return gas_f.v_visc(Sigma, nu, sim.grid.r, sim.grid.ri)
