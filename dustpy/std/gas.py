"""Module containing standard functions for the gas."""

import numpy as np
import scipy.sparse as sp
from simframe.integration import Scheme
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
    sim.gas.Sigma[:] = gas_f.enforce_floor(
        sim.gas.Sigma,
        sim.gas.SigmaFloor
    )


def prepare(sim):
    """Function prepares gas integration step.
    It stores the current value of the surface density in a hidden field.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame"""
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
    ret = gas_f.implicit_boundaries(
        sim.t.prevstepsize,
        sim.gas.Fi,
        sim.grid.ri,
        sim.gas.Sigma,
        sim.gas._SigmaOld
    )

    # Source terms
    sim.gas.S.tot[0] = ret[0]
    sim.gas.S.hyd[0] = ret[0]
    sim.gas.S.tot[-1] = ret[1]
    sim.gas.S.hyd[-1] = ret[1]

    # Fluxes through boundaries
    sim.gas.Fi[0] = ret[2]
    sim.gas.Fi[-1] = ret[3]


def dt(sim):
    """Function calculates the time step from the gas sources.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    dt : float
        Gas time step"""
    return gas_f.timestep(
        sim.gas.S.tot,
        sim.gas.Sigma,
        sim.gas.SigmaFloor
    )


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
    return gas_f.cs_adiabatic(
        sim.gas.gamma,
        sim.gas.mu,
        sim.gas.T
    )


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
    return gas_f.hp(
        sim.gas.cs,
        sim.gas.gamma,
        sim.grid.OmegaK
    )


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
    v += sim.gas.torque.v

    # Helper variables for convenience
    r = sim.grid.r
    ri = sim.grid.ri
    area = sim.grid.A
    Nr = int(sim.grid.Nr)

    # Construct Jacobian
    A, B, C = gas_f.jac_abc(area, nu, r, ri, v)
    row_hyd = np.hstack(
        (np.arange(Nr-1)+1, np.arange(Nr), np.arange(Nr-1)))
    col_hyd = np.hstack(
        (np.arange(Nr-1), np.arange(Nr), np.arange(Nr-1)+1))
    dat_hyd = np.hstack((A.ravel()[1:], B.ravel(), C.ravel()[:-1]))

    # Right hand side
    sim.gas._rhs[:] = sim.gas.Sigma

    # Boundaries. This is only reserving space in the sparce matrix
    row_in = [0, 0, 0]
    col_in = [0, 1, 2]
    dat_in = [0., 0., 0.]
    row_out = [Nr-1, Nr-1, Nr-1]
    col_out = [Nr-3, Nr-2, Nr-1]
    dat_out = [0., 0., 0.]

    # Stitching together the generators
    row = np.hstack((row_hyd, row_in, row_out))
    col = np.hstack((col_hyd, col_in, col_out))
    dat = np.hstack((dat_hyd, dat_in, dat_out))
    gen = (dat, (row, col))
    # Building sparse matrix of coagulation Jacobian
    J = sp.csc_matrix(
        gen,
        shape=(Nr, Nr)
    )

    return J


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
    return gas_f.mfp_midplane(
        sim.gas.n
    )


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
    return gas_f.n_midplane(
        sim.gas.mu,
        sim.gas.rho
    )


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
    return gas_f.viscosity(
        sim.gas.alpha,
        sim.gas.cs,
        sim.gas.Hp
    )


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
    return gas_f.p_midplane(
        sim.gas.cs,
        sim.gas.gamma,
        sim.gas.rho
    )


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
    return gas_f.rho_midplane(
        sim.gas.Hp,
        sim.gas.Sigma
    )


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
    return gas_f.s_tot(
        sim.gas.S.ext,
        sim.gas.S.hyd
    )


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
    return gas_f.t_pass(
        sim.star.L,
        sim.grid.r
    )


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
    return gas_f.v_rad(
        sim.dust.backreaction.A,
        sim.dust.backreaction.B,
        sim.gas.eta,
        sim.grid.OmegaK,
        sim.grid.r,
        sim.gas.v.visc
    )


def vtorque(sim):
    """Function calculates the velocity profile due by a torque profile

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    vtorque : array
        Velocity profile"""
    return 2. * sim.gas.torque.Lambda / (sim.grid.OmegaK * sim.grid.r)


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


def _f_impl_1_direct(x0, Y0, dx, *args, **kwargs):
    """Implicit 1st-order Euler integration scheme with direct matrix inversion

    Parameters
    ----------
    x0 : Intvar
        Integration variable at beginning of scheme
    Y0 : Field
        Variable to be integrated at the beginning of scheme
    dx : IntVar
        Stepsize of integration variable
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
    # Getting keyword arguments. Default is standard gas.
    boundary = kwargs.get("boundary", Y0._owner.gas.boundary)
    Sext = kwargs.get("Sext", Y0._owner.gas.S.ext)

    jac = Y0.jacobian(x0, dx)
    rhs = np.array(Y0)

    # Setting boundary values in jac and rhs

    # Inner boundary
    if boundary.inner is not None:
        # Given value
        if boundary.inner.condition == "val":
            rhs[0] = boundary.inner.value
        # Constant value
        elif boundary.inner.condition == "const_val":
            jac[0, 1] = 1./dx
            rhs[0] = 0.
        # Given gradient
        elif boundary.inner.condition == "grad":
            K1 = - boundary.inner._r[1]/boundary.inner._r[0]
            jac[0, 1] = -K1/dx
            rhs[0] = - boundary.inner._ri[1]/boundary.inner._r[0] * \
                (boundary.inner._r[1]-boundary.inner._r[0]) * \
                boundary.inner.value
        # Constant gradient
        elif boundary.inner.condition == "const_grad":
            Di = boundary.inner._ri[1]/boundary.inner._ri[2] * (
                boundary.inner._r[1]-boundary.inner._r[0]) / (boundary.inner._r[2]-boundary.inner._r[0])
            K1 = - boundary.inner._r[1]/boundary.inner._r[0] * (1. + Di)
            K2 = boundary.inner._r[2]/boundary.inner._r[0] * Di
            jac[0, :3] = 0.
            jac[0, 1] = -K1/dx
            jac[0, 2] = -K2/dx
            rhs[0] = 0.
        # Given power law
        elif boundary.inner.condition == "pow":
            p = boundary.inner.value
            rhs[0] = Y0[1] * (boundary.inner._r[0]/boundary.inner._r[1])**p
        # Constant power law
        elif boundary.inner.condition == "const_pow":
            p = np.log(Y0[2] / Y0[1]) / \
                np.log(boundary.inner._r[2]/boundary.inner._r[1])
            K1 = - (boundary.inner._r[0]/boundary.inner._r[1])**p
            jac[0, 1] = -K1/dx
            rhs[0] = 0.

    # Outer boundary
    if boundary.outer is not None:
        # Given value
        if boundary.outer.condition == "val":
            rhs[-1] = boundary.outer.value
        # Constant value
        elif boundary.outer.condition == "const_val":
            jac[-1, -2] = (1./dx)
            rhs[-1] = 0.
        # Given gradient
        elif boundary.outer.condition == "grad":
            KNrm2 = - boundary.outer._r[1]/boundary.outer._r[0]
            jac[-1, -2] = -(KNrm2/dx)
            rhs[-1] = boundary.outer._ri[1]/boundary.outer._r[0] * \
                (boundary.outer._r[0]-boundary.outer._r[1]) * \
                boundary.outer.value
        # Constant gradient
        elif boundary.outer.condition == "const_grad":
            Do = boundary.outer._ri[1]/boundary.outer._ri[2] * (
                boundary.outer._r[0]-boundary.outer._r[1]) / (boundary.outer._r[1]-boundary.outer._r[2])
            KNrm2 = - boundary.outer._r[1]/boundary.outer._r[0] * (1. + Do)
            KNrm3 = boundary.outer._r[2]/boundary.outer._r[0] * Do
            jac[-1, -2] = -KNrm2/dx
            jac[-1, -3] = -KNrm3/dx
            rhs[-1] = 0.
        # Given power law
        elif boundary.outer.condition == "pow":
            p = boundary.outer.value
            rhs[-1] = Y0[-2] * (boundary.outer._r[-0]/boundary.outer._r[1])**p
        # Constant power law
        elif boundary.outer.condition == "const_pow":
            p = np.log(Y0[-2] / Y0[-3]) / \
                np.log(boundary.outer._r[1]/boundary.outer._r[2])
            KNrm2 = - (boundary.outer._r[0]/boundary.outer._r[1])**p
            jac[-1, -2] = -KNrm2/dx
            rhs[-1] = 0.

    # Add external source terms to right-hand side
    rhs[:] = gas_f.modified_rhs(
        dx,
        rhs,
        Sext
    )

    jac.data[:] = gas_f.modified_jacobian(
        dx,
        jac.data,
        jac.indices,
        jac.indptr
    )

    A_LU = sp.linalg.splu(jac,
                          permc_spec="MMD_AT_PLUS_A",
                          diag_pivot_thresh=0.0,
                          options=dict(SymmetricMode=True))
    Y1 = A_LU.solve(rhs)

    return Y1 - Y0


class impl_1_direct(Scheme):
    """Modified class for implicit gas integration."""

    def __init__(self):
        super().__init__(_f_impl_1_direct, description="Implicit 1st-order direct solver")
