'''Module containing standard functions for the dust.'''


import dustpy.constants as c
from dustpy.std import dust_f

import numpy as np
import scipy.sparse as sp

from simframe.integration import Scheme


def boundary(sim):
    """Function sets the boundary condition of dust surface density.
    Not implemented, yet.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame"""
    sim.dust.boundary.inner.setboundary()
    sim.dust.boundary.outer.setboundary()


def enforce_floor_value(sim):
    """Function enforces floor value onto dust surface density.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame"""
    sim.dust.Sigma = np.where(
        sim.dust.Sigma > sim.dust.SigmaFloor,
        sim.dust.Sigma,
        0.1*sim.dust.SigmaFloor)


def prepare(sim):
    """Function prepares implicit dust integration step.
    It stores the current value of the surface density in a hidden field.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame"""
    # Setting external sources at boundaries to zero
    sim.dust.S.ext[0] = 0.
    sim.dust.S.ext[-1] = 0.
    # Storing current surface density
    sim.dust._SigmaOld[...] = sim.dust.Sigma[...]


def finalize_explicit(sim):
    """Function finalizes integration step.

    Parameters
    ----------
    sim : Frame
        Parent integration frame"""
    boundary(sim)
    enforce_floor_value(sim)


def finalize_implicit(sim):
    """Function finalizes implicit integration step.

    Parameters
    ----------
    sim : Frame
        Parent integration frame"""
    boundary(sim)
    enforce_floor_value(sim)
    sim.dust.v.rad.update()
    sim.dust.Fi.update()
    sim.dust.S.hyd.update()
    sim.dust.S.coag.update()
    set_implicit_boundaries(sim)


def set_implicit_boundaries(sim):
    """Function calculates the fluxes at the boundaries after the implicit integration step.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame"""
    # Total source terms
    sim.dust.S.tot[0] = (sim.dust.Sigma[0] -
                         sim.dust._SigmaOld[0])/(sim.t.prevstepsize+1.e-100)
    sim.dust.S.tot[-1] = (sim.dust.Sigma[-1] -
                          sim.dust._SigmaOld[-1])/(sim.t.prevstepsize+1.e-100)
    # Hydrodynamic source terms
    sim.dust.S.hyd[0] = sim.dust.S.tot[0]
    sim.dust.S.hyd[-1] = sim.dust.S.tot[-1]
    # Fluxes
    sim.dust.Fi.adv[0] = (0.5*sim.dust.S.hyd[0]*(sim.grid.ri[1]**2 -
                                                 sim.grid.ri[0]**2) + sim.grid.ri[1]*sim.dust.Fi.adv[1])/sim.grid.ri[0]
    sim.dust.Fi.adv[-1] = (sim.dust.Fi.adv[-2]*sim.grid.ri[-2] - 0.5*sim.dust.S.hyd[-1]
                           * (sim.grid.ri[-1]**2-sim.grid.ri[-2]**2))/sim.grid.ri[-1]
    sim.dust.Fi.tot[0] = sim.dust.Fi.adv[0]
    sim.dust.Fi.tot[-1] = sim.dust.Fi.adv[-1]


def dt_adaptive(sim):
    """Function returns the adaptive time step.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    dt : float
        Dust time step"""
    return sim.t.suggested


def dt(sim):
    """Function calculates the time step from the dust sources.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    dt : float
        Dust time step"""
    if np.any(sim.dust.S.tot[1:-1, ...] < 0.):
        mask = np.logical_and(
            sim.dust.Sigma > sim.dust.SigmaFloor,
            sim.dust.S.tot < 0.)
        mask[0, :] = False
        mask[-1:, :] = False
        rate = sim.dust.Sigma[mask] / sim.dust.S.tot[mask]
        try:
            return np.min(np.abs(rate))
        except:
            return None


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
        Dust diffusivity

    Notes
    -----
    The diffusivity at the first and last two radial
    grid cells will be set to zero to avoid unwanted
    behavior at the boundaries."""
    v2 = sim.dust.delta.rad * sim.gas.cs**2
    Diff = dust_f.d(v2, sim.grid.OmegaK, sim.dust.St)
    Diff[:2, ...] = 0.
    Diff[-2:, ...] = 0.
    return Diff


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

    Fi = dust_f.fi_diff(sim.dust.D,
                        Sigma,
                        sim.gas.Sigma,
                        sim.dust.St,
                        np.sqrt(sim.dust.delta.rad*sim.gas.cs**2),
                        sim.grid.r,
                        sim.grid.ri)
    Fi[:1, :] = 0.
    Fi[-1:, :] = 0.

    return Fi


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
    Fi = np.zeros((int(sim.grid.Nr+1), int(sim.grid.Nm)))
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


def jacobian(sim, x, dx=None, *args, **kwargs):
    """Function calculates the Jacobian for implicit dust integration.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame
    x : IntVar
        Integration variable
    dx : float, optional, default : None
        stepsize
    args : additional positional arguments
    kwargs : additional keyworda arguments

    Returns
    -------
    jac : Sparse matrix
        Dust Jacobian

    Notes
    -----
    Function returns the Jacobian for ``Simulation.dust.Sigma.ravel()``
    instead of ``Simulation.dust.Sigma``. The Jacobian is stored as
    sparse matrix."""

    # Parameters for function call
    A = sim.dust.coagulation.A
    cstick = sim.dust.coagulation.stick
    eps = sim.dust.coagulation.eps
    ilf = sim.dust.coagulation.lf_ind
    irm = sim.dust.coagulation.rm_ind
    istick = sim.dust.coagulation.stick_ind
    m = sim.grid.m
    phi = sim.dust.coagulation.phi
    Rf = sim.dust.kernel * sim.dust.p.frag
    Rs = sim.dust.kernel * sim.dust.p.stick
    SigD = sim.dust.Sigma
    SigDfloor = sim.dust.SigmaFloor

    # Helper variables for convenience
    if dx is None:
        dt = x.stepsize
    else:
        dt = dx
    r = sim.grid.r
    ri = sim.grid.ri
    area = sim.grid.A
    Nr = int(sim.grid.Nr)
    Nm = int(sim.grid.Nm)

    # Building coagulation Jacobian

    # Total problem size
    Ntot = int((Nr*Nm))
    # Getting data vector and coordinates in sparse matrix
    dat, row, col = dust_f.jacobian_coagulation_generator(
        A, cstick, eps, ilf, irm, istick, m, phi, Rf, Rs, SigD, SigDfloor)
    gen = (dat, (row, col))
    # Building sparse matrix of coagulation Jacobian
    J_coag = sp.csc_matrix(
        gen,
        shape=(Ntot, Ntot)
    )

    # Building the hydrodynamic Jacobian
    A, B, C = dust_f.jacobian_hydrodynamic_generator(
        area,
        sim.dust.D,
        r,
        ri,
        sim.gas.Sigma,
        sim.dust.v.rad
    )
    J_hyd = sp.diags(
        (A.ravel()[Nm:], B.ravel(), C.ravel()[:-Nm]),
        offsets=(-Nm, 0, Nm),
        shape=(Ntot, Ntot),
        format="csc"
    )

    # Right-hand side
    sim.dust._rhs[Nm:-Nm] = sim.dust.Sigma.ravel()[Nm:-Nm]

    # BOUNDARIES

    # Inner boundary

    # Initializing data and coordinate vectors for sparse matrix
    dat = np.zeros(int(3.*Nm))
    row0 = np.arange(int(Nm))
    col0 = np.arange(int(Nm))
    col1 = np.arange(int(Nm)) + Nm
    col2 = np.arange(int(Nm)) + 2.*Nm
    row = np.concatenate((row0, row0, row0))
    col = np.concatenate((col0, col1, col2))

    # Filling data vector depending on boundary condition
    if sim.dust.boundary.inner is not None:
        # Given value
        if sim.dust.boundary.inner.condition == "val":
            sim.dust._rhs[:Nm] = sim.dust.boundary.inner.value
        # Constant value
        elif sim.dust.boundary.inner.condition == "const_val":
            dat[Nm:2*Nm] = 1./dt
            sim.dust._rhs[:Nm] = 0.
        # Given gradient
        elif sim.dust.boundary.inner.condition == "grad":
            K1 = - r[1]/r[0]
            dat[Nm:2*Nm] = -K1/dt
            sim.dust._rhs[:Nm] = - ri[1]/r[0] * \
                (r[1]-r[0])*sim.dust.boundary.inner.value
        # Constant gradient
        elif sim.dust.boundary.inner.condition == "const_grad":
            Di = ri[1]/ri[2] * (r[1]-r[0]) / (r[2]-r[0])
            K1 = - r[1]/r[0] * (1. + Di)
            K2 = r[2]/r[0] * Di
            dat[:Nm] = 0.
            dat[Nm:2*Nm] = -K1/dt
            dat[2*Nm:] = -K2/dt
            sim.dust._rhs[:Nm] = 0.
        # Given power law
        elif sim.dust.boundary.inner.condition == "pow":
            p = sim.dust.boundary.inner.value
            sim.dust._rhs[:Nm] = sim.dust.Sigma[1] * (r[0]/r[1])**p
        # Constant power law
        elif sim.dust.boundary.inner.condition == "const_pow":
            p = np.log(sim.dust.Sigma[2] /
                       sim.dust.Sigma[1]) / np.log(r[2]/r[1])
            K1 = - (r[0]/r[1])**p
            dat[Nm:2*Nm] = -K1/dt
            sim.dust._rhs[:Nm] = 0.

    # Creating sparce matrix for inner boundary
    gen = (dat, (row, col))
    J_in = sp.csc_matrix(
        gen,
        shape=(Ntot, Ntot)
    )

    # Outer boundary

    # Initializing data and coordinate vectors for sparse matrix
    dat = np.zeros(int(3.*Nm))
    row0 = np.arange(int(Nm))
    col0 = np.arange(int(Nm))
    col1 = np.arange(int(Nm)) - Nm
    col2 = np.arange(int(Nm)) - 2.*Nm
    offset = (Nr-1)*Nm
    row = np.concatenate((row0, row0, row0)) + offset
    col = np.concatenate((col0, col1, col2)) + offset

    # Filling data vector depending on boundary condition
    if sim.dust.boundary.outer is not None:
        # Given value
        if sim.dust.boundary.outer.condition == "val":
            sim.dust._rhs[-Nm:] = sim.dust.boundary.outer.value
        # Constant value
        elif sim.dust.boundary.outer.condition == "const_val":
            dat[-2*Nm:-Nm] = 1./dt
            sim.dust._rhs[-Nm:] = 0.
        # Given gradient
        elif sim.dust.boundary.outer.condition == "grad":
            KNrm2 = -r[-2]/r[-1]
            dat[-2*Nm:-Nm] = -KNrm2/dt
            sim.dust._rhs[-Nm:] = ri[-2]/r[-1] * \
                (r[-1]-r[-2])*sim.dust.boundary.outer.value
        # Constant gradient
        elif sim.dust.boundary.outer.condition == "const_grad":
            Do = ri[-2]/ri[-3] * (r[-1]-r[-2]) / (r[-2]-r[-3])
            KNrm2 = - r[-2]/r[-1] * (1. + Do)
            KNrm3 = r[-3]/r[-1] * Do
            dat[-2*Nm:-Nm] = -KNrm2/dt
            dat[-3*Nm:-2*Nm] = -KNrm3/dt
            sim.dust._rhs[-Nm:] = 0.
        # Given power law
        elif sim.dust.boundary.outer.condition == "pow":
            p = sim.dust.boundary.outer.value
            sim.dust._rhs[-Nm:] = sim.dust.Sigma[-2] * (r[-1]/r[-2])**p
        # Constant power law
        elif sim.dust.boundary.outer.condition == "const_pow":
            p = np.log(sim.dust.Sigma[-2] /
                       sim.dust.Sigma[-3]) / np.log(r[-2]/r[-3])
            KNrm2 = - (r[-1]/r[-2])**p
            dat[-2*Nm:-Nm] = -KNrm2/dt
            sim.dust._rhs[-Nm:] = 0.

    # Creating sparce matrix for outer boundary
    gen = (dat, (row, col))
    J_out = sp.csc_matrix(
        gen,
        shape=(Ntot, Ntot)
    )

    # Adding and returning all matrix components
    return J_in + J_coag + J_hyd + J_out


def kernel(sim):
    """Function calculates the vertically integrated collision kernel.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    K : Field
        Collision kernel"""
    return dust_f.kernel(sim.dust.a,
                         sim.dust.H,
                         sim.dust.Sigma,
                         sim.dust.SigmaFloor,
                         sim.dust.v.rel.tot)


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
    # Set maximum particle size
    if(sim.ini.dust.allowDriftingParticles):
        aIni = sim.ini.dust.aIniMax
    else:
        # Calculating pressure gradient
        P = sim.gas.P
        Pi = dust_f.interp1d(sim.grid.ri, sim.grid.r, P)
        gamma = (Pi[1:] - Pi[:-1]) / (sim.grid.ri[1:] - sim.grid.ri[:-1])
        gamma = np.abs(gamma)
        # Exponent of pressure gradient
        gamma *= sim.grid.r / P
        gamma = 1. / gamma
        # Maximum drift limited particle size with safety margin
        ad = 1.e-4 * 2./np.pi * sim.ini.dust.d2gRatio * sim.gas.Sigma \
            / sim.dust.fill[:, 0] * sim.dust.rhos[:, 0] * (sim.grid.OmegaK * sim.grid.r)**2. \
            / sim.gas.cs**2. / gamma
        aIni = np.minimum(sim.ini.dust.aIniMax, ad)[:, None]
    # Fill distribution
    ret = np.where(sim.dust.a <= aIni, sim.dust.a**(exp+4), 0.)
    s = np.sum(ret, axis=1)[..., None]
    s = np.where(s > 0., s, 1.)
    # Normalize to mass
    ret = ret / s * sim.gas.Sigma[..., None] * sim.ini.dust.d2gRatio
    return ret


def p_stick(sim):
    """Function calculates the sticking probability.
    The sticking probability is simply 1 minus the
    fragmentation probability.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    ps : Field
        Sticking probability"""
    p = 1. - sim.dust.p.frag
    p[0] = 0.
    p[-1] = 0.
    return p


def p_frag(sim):
    """Function calculates the fragmentation probability.
    It assumes a linear transition between sticking and
    fragmentation.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    pf : Field
        Fragmentation propability."""
    return dust_f.pfrag(sim.dust.v.rel.tot, sim.dust.v.frag)


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
    if Sigma is None:
        Sigma = sim.dust.Sigma
    return dust_f.s_coag(sim.dust.coagulation.stick,
                         sim.dust.coagulation.stick_ind,
                         sim.dust.coagulation.A,
                         sim.dust.coagulation.eps,
                         sim.dust.coagulation.lf_ind,
                         sim.dust.coagulation.rm_ind,
                         sim.dust.coagulation.phi,
                         sim.dust.kernel * sim.dust.p.frag,
                         sim.dust.kernel * sim.dust.p.stick,
                         sim.grid.m,
                         Sigma,
                         sim.dust.SigmaFloor)


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
    Sext = sim.dust.S.ext
    if Sigma is None:
        Sigma = sim.dust.Sigma
        Scoag = sim.dust.S.coag
        Shyd = sim.dust.S.hyd
    else:
        Scoag = sim.dust.S.coag.updater.beat(sim, Sigma=Sigma)
        if Scoag is None:
            Scoag = sim.dust.S.coag
        Shyd = sim.dust.S.hyd.updater.beat(sim, Sigma=Sigma)
        if Shyd is None:
            Shyd = sim.dust.S.hyd
    return Scoag + Shyd + Sext


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


def coagulation_parameters(sim):
    """Function calculates the coagulation parameters needed for a simple
    sticking-erosion-fragmentation collision model. The sticking matrix
    is calculated with the method described in appendices A.1. and A.2. of
    Brauer et al. (2008).

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    (cstick, cstick_ind, A, eps, klf, krm, phi) : Tuple
        Coagulation parameters

    Notes
    -----
    The sticking matrix has technically a shape of ``(Nm, Nm, Nm)``.
    For example: ``cstick(k, j, i)`` describes the change of mass bin ``k``
    resulting from a sticking collision between particles ``j`` and ``k``.
    Since this matrix has at maximum four entries per particle collision,
    only the non-zero elemts are stored in ``cstick`` of shape ``(4, Nm, Nm)``.
    The positions of the non-zero elements along the first axis are stored
    in ``cstick_ind``. For details see Brauer et al. (2008)."""
    cstick, cstick_ind, A, eps, klf, krm, phi = dust_f.coagulation_parameters(sim.ini.dust.erosionMassRatio,
                                                                              sim.ini.dust.excavatedMass,
                                                                              sim.ini.dust.fragmentDistribution,
                                                                              sim.grid.m,
                                                                              int(sim.grid.Nr))
    return cstick, cstick_ind, A, eps, klf, krm, phi


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
    It uses the prescription of Ormel & Cuzzi (2007).

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    vrel : Field
        Relative velocities"""
    return dust_f.vrel_ormel_cuzzi_2007(
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


def _f_impl_1_direct(x0, Y0, dx, jac=None, rhs=None, *args, **kwargs):
    """Implicit 1st-order integration scheme with direct matrix inversion

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
        jac = Y0.jacobian(x0, dx)
    if rhs is None:
        rhs = np.array(Y0.ravel())

    Nm = Y0._owner.dust.Sigma.shape[1]

    # Add external source terms to right-hand side
    rhs[Nm:-Nm] += dx*Y0._owner.dust.S.ext[1:-1, ...].ravel()

    N = jac.shape[0]
    eye = sp.identity(N, format="csc")

    A = eye - dx[0] * jac

    A_LU = sp.linalg.splu(A,
                          permc_spec="MMD_AT_PLUS_A",
                          diag_pivot_thresh=0.0,
                          options=dict(SymmetricMode=True))
    Y1_ravel = A_LU.solve(rhs)

    Y1 = Y1_ravel.reshape(Y0.shape)

    return Y1 - Y0


class impl_1_direct(Scheme):
    """Modified class for implicit dust integration."""

    def __init__(self):
        super().__init__(_f_impl_1_direct, description="Implicit 1st-order direct solver")
