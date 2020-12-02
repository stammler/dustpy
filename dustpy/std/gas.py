'''Module containing standard functions for the gas.'''

import numpy as np
from simframe.integration import Scheme

import dustpy.constants as c
from dustpy.std import gas_f


def boundary(sim):
    """Function set the boundary conditions of the gas.
    Not implemented, yet.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame"""
    sim.gas.boundary.inner.setboundary()
    sim.gas.boundary.outer.setboundary()


def enforce_floor_value(sim):
    """Function enforces floor value to gas surface density.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame"""
    sim.gas.Sigma[:] = np.maximum(sim.gas.Sigma, sim.gas.SigmaFloor)


def prepare(sim):
    """Function prepares gas integration step.
    It stores the current value of the surface density in a hidden field.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame"""
    # Setting external sources at boundaries to zero
    sim.gas.S.ext[0] = 0.
    sim.gas.S.ext[-1] = 0.
    # Storing current surface density
    sim.gas._SigmaOld[:] = sim.gas.Sigma[:]


def finalize(sim):
    """Function finalizes gas integration step.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame"""
    boundary(sim)
    enforce_floor_value(sim)
    sim.gas.v.update()
    sim.gas.Fi.update()
    sim.gas.S.hyd.update()
    set_implicit_boundaries(sim)


def set_implicit_boundaries(sim):
    """Function calculates the fluxes at the boundaries after the implicit integration step.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame"""
    # Total source terms
    sim.gas.S.tot[0] = (sim.gas.Sigma[0]-sim.gas._SigmaOld[0])/sim.t.stepsize
    sim.gas.S.tot[-1] = (sim.gas.Sigma[-1] -
                         sim.gas._SigmaOld[-1])/sim.t.stepsize
    # Hydrodynamic source terms
    sim.gas.S.hyd[0] = sim.gas.S.tot[0]
    sim.gas.S.hyd[-1] = sim.gas.S.tot[-1]
    # Fluxes
    sim.gas.Fi[0] = (0.5*sim.gas.S.hyd[0]*(sim.grid.ri[1]**2 -
                                           sim.grid.ri[0]**2) + sim.grid.ri[1]*sim.gas.Fi[1])/sim.grid.ri[0]
    sim.gas.Fi[-1] = (sim.gas.Fi[-2]*sim.grid.ri[-2] - 0.5*sim.gas.S.hyd[-1]
                      * (sim.grid.ri[-1]**2-sim.grid.ri[-2]**2))/sim.grid.ri[-1]


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


def Fi(sim):
    """Function calculates the mass flux through the cell interfaces.
    The fluxes are calculated from the implicit integration outcome.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    Fi : Field
        Mass flux through grid cell interfaces"""
    return gas_f.fi(sim.gas.Sigma, sim.gas.v.rad, sim.grid.r, sim.grid.ri)


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
    return sim.gas.cs/(np.sqrt(sim.gas.gamma)*sim.grid.OmegaK)


def jacobian(sim, x, *args, **kwargs):
    """Functions calculates the Jacobian for the gas.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame
    x : IntVar
        Integration variable
    args : additional positional arguments
    kwargs : additional keyworda arguments

    Returns
    -------
    jac : Field
        Jacobi matrix for gas evolution

    Notes
    -----
    The boundaries need information about the time step, which is only available at
    the integration stage. The boundary values are therefore meaningless and should not
    be used to calculate the source terms via matrix-vector multiplication.
    See the documentation for details."""

    # Parameters
    nu = sim.gas.nu * sim.dust.backreaction.A
    v = sim.dust.backreaction.B * 2. * sim.gas.eta * sim.grid.r * sim.grid.OmegaK

    # Helper variables for convenience
    dt = sim.t.stepsize
    r = sim.grid.r
    ri = sim.grid.ri
    area = sim.grid.A

    # Construct Jacobian
    A, B, C = gas_f.jac_abc(area, nu, r, ri, v)
    jac = np.diag(A[1:], -1) + np.diag(B) + np.diag(C[:-1], 1)

    # Right hand side
    sim.gas._rhs[:] = sim.gas.Sigma

    # Boundaries

    # Inner boundary
    if sim.gas.boundary.inner is not None:
        # Given value
        if sim.gas.boundary.inner.condition == "val":
            jac[0, 0] = 0.
            sim.gas._rhs[0] = sim.gas.boundary.inner.value
        # Constant value
        elif sim.gas.boundary.inner.condition == "const_val":
            jac[0, 0] = 0.
            jac[0, 1] = 1./dt
            sim.gas._rhs[0] = 0.
        # Given gradient
        elif sim.gas.boundary.inner.condition == "grad":
            K1 = - r[1]/r[0]
            jac[0, 0] = 0.
            jac[0, 1] = -K1/dt
            sim.gas._rhs[0] = - ri[1]/r[0] * \
                (r[1]-r[0])*sim.gas.boundary.inner.value
        # Constant gradient
        elif sim.gas.boundary.inner.condition == "const_grad":
            Di = ri[1]/ri[2] * (r[1]-r[0]) / (r[2]-r[0])
            K1 = - r[1]/r[0] * (1. + Di)
            K2 = r[2]/r[0] * Di
            jac[0, 0] = 0.
            jac[0, 1] = -K1/dt
            jac[0, 2] = -K2/dt
            sim.gas._rhs[0] = 0.
        # Given power law
        elif sim.gas.boundary.inner.condition == "pow":
            p = sim.gas.boundary.inner.value
            jac[0, 0] = 0.
            sim.gas._rhs[0] = sim.gas.Sigma[1] * (r[0]/r[1])**p
        # Constant power law
        elif sim.gas.boundary.inner.condition == "const_pow":
            p = np.log(sim.gas.Sigma[2] /
                       sim.gas.Sigma[1]) / np.log(r[2]/r[1])
            K1 = - (r[0]/r[1])**p
            jac[0, 0] = 0.
            jac[0, 1] = -K1/dt
            sim.gas._rhs[0] = 0.

    # Outer boundary
    if sim.gas.boundary.outer is not None:
        # Given value
        if sim.gas.boundary.outer.condition == "val":
            jac[-1, -1] = 0
            sim.gas._rhs[-1] = sim.gas.boundary.outer.value
        # Constant value
        elif sim.gas.boundary.outer.condition == "const_val":
            jac[-1, -1] = 0
            jac[-1, -2] = 1./dt
            sim.gas._rhs[-1] = 0.
        # Given gradient
        elif sim.gas.boundary.outer.condition == "grad":
            KNrm2 = - r[-2]/r[-1]
            jac[-1, -1] = 0.
            jac[-1, -2] = -KNrm2/dt
            sim.gas._rhs[-1] = ri[-2]/r[-1] * \
                (r[-1]-r[-2])*sim.gas.boundary.outer.value
        # Constant gradient
        elif sim.gas.boundary.outer.condition == "const_grad":
            Do = ri[-2]/ri[-3] * (r[-1]-r[-2]) / (r[-2]-r[-3])
            KNrm2 = - r[-2]/r[-1] * (1. + Do)
            KNrm3 = r[-3]/r[-1] * Do
            jac[-1, -1] = 0.
            jac[-1, -2] = -KNrm2/dt
            jac[-1, -3] = -KNrm3/dt
            sim.gas._rhs[-1] = 0.
        # Given power law
        elif sim.gas.boundary.outer.condition == "pow":
            p = sim.gas.boundary.outer.value
            jac[-1, -1] = 0.
            sim.gas._rhs[-1] = sim.gas.Sigma[-2] * (r[-1]/r[-2])**p
        # Constant power law
        elif sim.gas.boundary.outer.condition == "const_pow":
            p = np.log(sim.gas.Sigma[-2] /
                       sim.gas.Sigma[-3]) / np.log(r[-2]/r[-3])
            KNrm2 = - (r[-1]/r[-2])**p
            jac[-1, -1] = 0.
            jac[-1, -2] = -KNrm2/dt
            sim.gas._rhs[-1] = 0.

    return jac


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
    Mdisk : float
        disk mass

    Returns
    -------
    Sigma : float or array of floats
        Surface density profile"""
    return (2+p)*Mdisk / (2.*np.pi*rc**2) * (r/rc)**p * np.exp(-(r/rc)**(2+p))


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


def nu(sim):
    """Function calculates the kinematic viscocity of the gas.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    nu : Field
        Kinematic viscosity"""
    return sim.gas.alpha * sim.gas.cs * sim.gas.Hp


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


def S_hyd(sim):
    """Function calculates the hydrodynamic source terms.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    S_hyd : Field
        Hydrodynamic source terms"""
    return gas_f.s_hyd(sim.gas.Fi, sim.grid.ri)


def S_tot(sim):
    """Function calculates the total source terms.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    S_tot : Field
        Total surface density source terms"""
    return sim.gas.S.hyd + sim.gas.S.ext


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


def vrad(sim):
    """Function calculates the radial radial gas velocity.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    vrad : Field
        Radial gas velocity"""
    vback = 2. * sim.gas.eta * sim.grid.r * sim.grid.OmegaK
    return sim.dust.backreaction.A * sim.gas.v.visc + sim.dust.backreaction.B * vback


def vvisc(sim):
    """Function calculates the viscous radial gas velocity.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    vvisc : Field
        Viscous radial gas velocity"""
    return gas_f.v_visc(sim.gas.Sigma, sim.gas.nu, sim.grid.r, sim.grid.ri)


def _f_impl_1_euler_direct(x0, Y0, dx, jac=None, rhs=None, *args, **kwargs):
    """Implicit 1st-order Euler integration scheme with direct matrix inversion

    Parameters
    ----------
    x0 : Intvar
        Integration variable at beginning of scheme
    Y0 : Field
        Variable to be integrated at the beginning of scheme
    dx : IntVar
        Stepsize of integration variable
    jac : Field, optional, defaul : None
        Current Jacobian. Will be calculated, if not set
    args : additional positional arguments
    kwargs : additional keyworda arguments

    Returns
    -------
    dY : Field
        Delta of variable to be integrated

    Butcher tableau
    ---------------
     1 | 1
    ---|---
       | 1 
    """
    if jac is None:
        jac = Y0.jacobian(x0 + dx, dt=dx)
    if rhs is None:
        rhs = np.array(Y0)

    # Add external source terms to right-hand side
    rhs[1:-1] += dx*Y0._owner.gas.S.ext[1:-1]

    N = jac.shape[0]
    eye = np.eye(N)
    A = eye - dx * jac
    Y1 = np.linalg.inv(A) @ rhs

    return Y1 - Y0


impl_1_euler_direct = Scheme(
    _f_impl_1_euler_direct, description="Implicit 1st-order direct Euler method")
