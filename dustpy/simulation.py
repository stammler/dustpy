from simframe import Frame
from simframe import Instruction
from simframe import Integrator
from simframe import schemes
from simframe.frame import Field
from simframe.frame import Group
from simframe.frame import IntVar
from simframe.io.writers import hdf5writer

import dustpy.constants as c

import dustpy.std.dust as std_dust
import dustpy.std.gas as std_gas
import dustpy.std.grid as std_grid
import dustpy.std.sim as std_sim
import dustpy.std.star as std_star

import numpy as np
from types import SimpleNamespace


class Simulation(Frame):

    __name__ = "DustPy"

    _ini = SimpleNamespace(**{"dust": SimpleNamespace(**{"aIniMax": 0.0001,
                                                         "allowDriftingParticles": False,
                                                         "crateringMassRatio": 10.,
                                                         "d2gRatio": 1.e-2,
                                                         "distExp": -3.5,
                                                         "excavatedMass": 1.,
                                                         "fragmentDistribution": -11/6,
                                                         "rhoMonomer": 1.67,
                                                         "vfrag": 100.
                                                         }
                                                      ),
                              "gas": SimpleNamespace(**{"alpha": 1.e-3,
                                                        "gamma": 1.4,
                                                        "Mdisk": 0.01*c.M_sun,
                                                        "mu": 2.3*c.m_p,
                                                        "SigmaExp": -1,
                                                        "SigmaRc": 30.*c.au
                                                        }
                                                     ),
                              "grid": SimpleNamespace(**{"Nmbpd": 7,
                                                         "mmin": 1.e-15,
                                                         "mmax": 1.e15,
                                                         "Nr": 100,
                                                         "rmin": 1.*c.au,
                                                         "rmax": 1000.*c.au
                                                         }
                                                      ),
                              "star": SimpleNamespace(**{"M": 1.*c.M_sun,
                                                         "R": 2.*c.R_sun,
                                                         "T": 5772.,
                                                         }
                                                      )
                              }
                           )

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # Initializing everything with None to get the basic frame

        # Dust quantities
        self.dust = Group(self, description="Dust quantities")
        self.dust.a = None
        self.dust.backreaction = Group(
            self, description="Backreaction coefficients")
        self.dust.backreaction.A = None
        self.dust.backreaction.B = None
        self.dust.backreaction.updater = ["A", "B"]
        self.dust.coagulation = Group(
            self, description="Coagulation quantities")
        self.dust.coagulation.matrices = Group(
            self, description="Collision matrices")
        self.dust.coagulation.matrices.coag = None
        self.dust.coagulation.matrices.frag = None
        self.dust.D = None
        self.dust.delta = Group(self, description="Mixing parameters")
        self.dust.delta.rad = None
        self.dust.delta.turb = None
        self.dust.delta.vert = None
        self.dust.delta.updater = ["rad", "turb", "vert"]
        self.dust.eps = None
        self.dust.Fi = Group(self, description="Fluxes")
        self.dust.Fi.adv = None
        self.dust.Fi.diff = None
        self.dust.Fi.tot = None
        self.dust.Fi.updater = ["adv", "diff", "tot"]
        self.dust.fill = None
        self.dust.H = None
        self.dust.rho = None
        self.dust.rhos = None
        self.dust.S = Group(self, description="Sources")
        self.dust.S.coag = None
        self.dust.S.ext = None
        self.dust.S.hyd = None
        self.dust.S.tot = None
        self.dust.S.updater = ["coag", "hyd", "ext", "tot"]
        self.dust.Sigma = None
        self.dust.SigmaFloor = None
        self.dust.St = None
        self.dust.v = Group(self, description="Velocities")
        self.dust.v.frag = None
        self.dust.v.rel = Group(self, description="Relative velocities")
        self.dust.v.rel.azi = None
        self.dust.v.rel.brown = None
        self.dust.v.rel.rad = None
        self.dust.v.rel.tot = None
        self.dust.v.rel.turb = None
        self.dust.v.rel.vert = None
        self.dust.v.rel.updater = [
            "azi", "brown", "rad", "turb", "vert", "tot"]
        self.dust.v.driftmax = None
        self.dust.v.rad = None
        self.dust.v.updater = ["frag", "driftmax", "rad", "rel"]
        self.dust.updater = ["delta", "rhos", "fill",
                             "a", "St", "H", "rho", "backreaction", "v", "D", "eps", "Fi", "S"]

        # Gas quantities
        self.gas = Group(self, description="Gas quantities")
        self.gas.alpha = None
        self.gas.cs = None
        self.gas.eta = None
        self.gas.Fi = None
        self.gas.gamma = None
        self.gas.Hp = None
        self.gas.mfp = None
        self.gas.mu = None
        self.gas.n = None
        self.gas.P = None
        self.gas.rho = None
        self.gas.S = Group(self, description="Source terms")
        self.gas.S.ext = None
        self.gas.S.hyd = None
        self.gas.S.tot = None
        self.gas.S.updater = ["hyd", "ext", "tot"]
        self.gas.Sigma = None
        self.gas.SigmaFloor = None
        self.gas.T = None
        self.gas.v = Group(self, description="Velocities")
        self.gas.v.rad = None
        self.gas.v.visc = None
        self.gas.v.updater = ["visc", "rad"]
        self.gas.updater = ["gamma", "mu", "T", "alpha", "cs", "Hp",
                            "rho", "n", "mfp", "P", "eta", "v", "Fi", "S"]

        # Grid quantities
        self.grid = Group(self, description="Grid quantities")
        self.grid.Nr = None
        self.grid.r = None
        self.grid.ri = None
        self.grid.Nm = None
        self.grid.m = None
        self.grid.OmegaK = None
        self.grid.updater = ["OmegaK"]

        # Stellar quantities
        self.star = Group(self, description="Stellar quantities")
        self.star.L = None
        self.star.M = None
        self.star.R = None
        self.star.T = None
        self.star.updater = ["M", "R", "T", "L"]

        # Updater of frame object
        self.updater = ["star", "grid", "gas", "dust"]

        self.t = None

    @property
    def ini(self):
        '''Parameter set for setting the initial conditions.'''
        return self._ini

    def makegrids(self):
        '''Function creates radial and mass grid.

        Notes
        -----
        The grids are set up with the parameters given in ``Simulation.ini``.
        If you want to have a custom radial grid you have to set the array of grid cell interfaces ``Simulation.grid.ri``,
        before calling ``Simulation.makegrids()``.
        You cannot have a customized mass grid, because the coagulation alorithm strictly needs a logarithmic grid. Always
        set your mass grid parameters with ``Simulation.ini``.'''
        self._makemassgrid()
        self._makeradialgrid()

    def _makemassgrid(self):
        '''Function sets the radial mass grid. If ``Simulation.grid.rint`` is not given it uses the parameters set in
        ``Simulation.ini``.'''
        Nm = np.int(np.ceil(np.log10(self.ini.grid.mmax /
                                     self.ini.grid.mmin)*self.ini.grid.Nmbpd))
        m = np.logspace(np.log10(self.ini.grid.mmin), np.log10(
            self.ini.grid.mmax), num=Nm, base=10.)
        self.grid.Nm = Field(
            self, Nm, description="# of mass bins", constant=True)
        self.grid.m = Field(
            self, m, description="Mass grid [g]", constant=True)

    def _makeradialgrid(self):
        '''Function sets the mass grid using the parameters set in ``Simulation.ini``.'''
        if self.grid.ri == None:
            ri = np.logspace(np.log10(self.ini.grid.rmin), np.log10(
                self.ini.grid.rmax), num=self.ini.grid.Nr+1, base=10.)
            Nr = self.ini.grid.Nr
        else:
            ri = self.grid.ri
            Nr = ri.shape[0] - 1
        r = 0.5*(ri[:-1] + ri[1:])
        self.grid.Nr = Field(
            self, Nr, description="# of radial grid cells", constant=True)
        self.grid.r = Field(
            self, r, description="Radial grid cell centers [cm]", constant=True)
        self.grid.ri = Field(
            self, ri, description="Radial grid cell interfaces [cm]", constant=True)

    def initialize(self):
        '''Function initializes the simulation frame.

        Function sets all fields that are None with a standard value.
        If the grids are not set, it will call ``Simulation.makegrids()`` first.'''
        if not isinstance(self.grid.Nm, Field) or not isinstance(self.grid.Nr, Field):
            self.makegrids()

        shape1 = (np.int(self.grid.Nr))
        shape1p1 = (np.int(self.grid.Nr)+1)
        shape2 = (np.int(self.grid.Nr), np.int(self.grid.Nm))
        shape2p1 = (np.int(self.grid.Nr)+1, np.int(self.grid.Nm))
        shape3 = (np.int(self.grid.Nr), np.int(
            self.grid.Nm), np.int(self.grid.Nm))

        # STELLAR QUANTITIES

        # Luminosity
        if self.star.L is None:
            self.star.L = Field(self, 0., description="Luminosity [erg/s]")
            self.star.L.updater = std_star.luminosity
            # self.star.L.update()
        # Mass
        if self.star.M is None:
            self.star.M = Field(self, self.ini.star.M, description="Mass [g]")
        # Radius
        if self.star.R is None:
            self.star.R = Field(self, self.ini.star.R,
                                description="Radius [cm]")
        # Effective temperature
        if self.star.T is None:
            self.star.T = Field(self, self.ini.star.T,
                                description="Effective temperature [K]")
        # Initialize stellar quantities
        self.star.update()

        # GRID QUANTITIES

        # Keplerian frequency
        if self.grid.OmegaK is None:
            self.grid.OmegaK = Field(self, np.zeros(
                shape1), description="Keplerian frequency [1/s]")
            self.grid.OmegaK.updater = std_grid.OmegaK
        # Initialize grid quantities
        self.grid.update()

        # GAS QUANTITIES

        # Turbulent alpha parameter
        if self.gas.alpha is None:
            alpha = self.ini.gas.alpha * np.ones(shape1)
            self.gas.alpha = Field(
                self, alpha, description="Turbulent alpha parameters")
        # Sound speed
        if self.gas.cs is None:
            self.gas.cs = Field(self, np.zeros(shape1),
                                description="Sound speed [cm/s]")
            self.gas.cs.updater = std_gas.cs_adiabatic
        # Pressure gradient parameter
        if self.gas.eta is None:
            self.gas.eta = Field(self, np.zeros(
                shape1), description="Pressure gradient parameter")
            self.gas.eta.updater = std_gas.eta_midplane
        # Gas flux at the cell interfaces
        if self.gas.Fi is None:
            self.gas.Fi = Field(self, np.zeros(shape1p1),
                                description="Gas flux interfaces [g/cm/s]")
            self.gas.Fi.updater = std_gas.Fi
        # Adiabatic index
        if self.gas.gamma is None:
            gamma = self.ini.gas.gamma * np.ones(shape1)
            self.gas.gamma = Field(self, gamma,
                                   description="Adiabatic index")
        # Pressure scale height
        if self.gas.Hp is None:
            self.gas.Hp = Field(self, np.zeros(shape1),
                                description="Pressure scale height [cm]")
            self.gas.Hp.updater = std_gas.Hp
        # Mean free path
        if self.gas.mfp is None:
            self.gas.mfp = Field(self, np.zeros(shape1),
                                 description="Midplane mean free path [cm]")
            self.gas.mfp.updater = std_gas.mfp_midplane
        # Mean molecular weight
        if self.gas.mu is None:
            mu = self.ini.gas.mu * np.ones(shape1)
            self.gas.mu = Field(self, mu,
                                description="Mean molecular weight [g]")
        # Midplane number density
        if self.gas.n is None:
            self.gas.n = Field(self, np.zeros(shape1),
                               description="Miplane number density [1/cm³]")
            self.gas.n.updater = std_gas.n_midplane
        # Midplane pressure
        if self.gas.P is None:
            self.gas.P = Field(self, np.zeros(shape1),
                               description="Midplane pressure [g/cm/s²]")
            self.gas.P.updater = std_gas.P_midplane
        # Midplane mass density
        if self.gas.rho is None:
            self.gas.rho = Field(self, np.zeros(shape1),
                                 description="Miplane mass density [g/cm³]")
            self.gas.rho.updater = std_gas.rho_midplane
        # Sources
        if self.gas.S.hyd is None:
            self.gas.S.hyd = Field(self, np.zeros(
                shape1), description="Hydrodynamic sources [g/cm²/s]")
            self.gas.S.hyd.updater = std_gas.S_hyd
        if self.gas.S.ext is None:
            self.gas.S.ext = Field(self, np.zeros(
                shape1), description="External sources [g/cm²/s]")
        if self.gas.S.tot is None:
            self.gas.S.tot = Field(self, np.zeros(
                shape1), description="Total sources [g/cm²/s]")
            self.gas.S.tot.updater = std_gas.S_tot
        # Floor value
        if self.gas.SigmaFloor is None:
            self.gas.SigmaFloor = Field(
                self, 1.e-100*np.ones(shape1), description="Floor value of surface density [g/cm²]")
        # Surface density
        if self.gas.Sigma is None:
            SigmaGas = np.array(std_gas.lyndenbellpringle1974(
                self.grid.r, self.ini.gas.SigmaRc, self.ini.gas.SigmaExp, self.ini.gas.Mdisk))
            SigmaGas = np.maximum(SigmaGas, self.gas.SigmaFloor)
            self.gas.Sigma = Field(self, SigmaGas,
                                   description="Surface density [g/cm²]")
            self.gas.Sigma.differentiator = std_gas.Sigma_deriv
        # Temperature
        if self.gas.T is None:
            self.gas.T = Field(self, np.zeros(shape1),
                               description="Temperature [K]")
            self.gas.T.updater = std_gas.T_passive
        # Velocities
        # Viscous accretion velocity
        if self.gas.v.visc is None:
            self.gas.v.visc = Field(self, np.zeros(shape1),
                                    description="Viscous accretion velocity [cm/s]")
            self.gas.v.visc.updater = std_gas.vvisc
        # Radial gas velocity
        if self.gas.v.rad is None:
            self.gas.v.rad = Field(self, np.zeros(shape1),
                                   description="Radial velocity [cm/s]")
            self.gas.v.rad.updater = std_gas.vrad
        # Initialize gas quantities
        self.gas.update()

        # DUST QUANTITIES

        # Particle size
        if self.dust.a is None:
            self.dust.a = Field(self, np.zeros(shape2),
                                description="Particle size [cm")
            self.dust.a.updater = std_dust.a
        # Backreaction
        if self.dust.backreaction.A is None:
            self.dust.backreaction.A = Field(
                self, np.ones(shape1), description="Pull factor")
        if self.dust.backreaction.B is None:
            self.dust.backreaction.B = Field(
                self, np.zeros(shape1), description="Push factor")
        # Diffusivity
        if self.dust.D is None:
            self.dust.D = Field(self, np.zeros(shape2),
                                description="Diffusivity [cm²/s]")
            self.dust.D.updater = std_dust.D
        # Deltas
        delta = self.ini.gas.alpha * np.ones(shape1)
        if self.dust.delta.rad is None:
            self.dust.delta.rad = Field(
                self, delta, description="Radial mixing parameter")
        if self.dust.delta.turb is None:
            self.dust.delta.turb = Field(
                self, delta, description="Turbulent mixing parameter")
        if self.dust.delta.vert is None:
            self.dust.delta.vert = Field(
                self, delta, description="Vertical mixing parameter")
        # Vertically integrated dust to gas ratio
        if self.dust.eps is None:
            self.dust.eps = Field(self, np.zeros(
                shape1), description="Dust-to-gas ratio")
            self.dust.eps.updater = std_dust.eps
        # Fluxes
        if self.dust.Fi.adv is None:
            self.dust.Fi.adv = Field(self, np.zeros(
                shape2p1), description="Advective flux [g/cm²/s")
            self.dust.Fi.adv.updater = std_dust.F_adv
        if self.dust.Fi.diff is None:
            self.dust.Fi.diff = Field(self, np.zeros(
                shape2p1), description="Diffusive flux [g/cm²/s")
            self.dust.Fi.diff.updater = std_dust.F_diff
        if self.dust.Fi.tot is None:
            self.dust.Fi.tot = Field(self, np.zeros(
                shape2p1), description="Total flux [g/cm²/s")
            self.dust.Fi.tot.updater = std_dust.F_tot
        # Filling factor
        if self.dust.fill is None:
            self.dust.fill = Field(self, np.ones(
                shape2), description="Filling factor")
        # Scale height
        if self.dust.H is None:
            self.dust.H = Field(self, np.zeros(shape2),
                                description="Scale heights [cm]")
            self.dust.H.updater = std_dust.H
        # Midplane mass density
        if self.dust.rho is None:
            self.dust.rho = Field(self, np.zeros(
                shape2), description="Midplane mass density per mass bin [g/cm³]")
            self.dust.rho.updater = std_dust.rho_midplane
        # Solid state density
        if self.dust.rhos is None:
            rhos = self.ini.dust.rhoMonomer * np.ones(shape2)
            self.dust.rhos = Field(
                self, rhos, description="Solid state density [g/cm³]")
        # Source terms
        if self.dust.S.coag is None:
            self.dust.S.coag = Field(self, np.zeros(
                shape2), description="Coagulation sources [g/cm²/s]")
            self.dust.S.coag.updater = std_dust.S_coag
        if self.dust.S.ext is None:
            self.dust.S.ext = Field(self, np.zeros(
                shape2), description="External sources [g/cm²/s]")
        if self.dust.S.hyd is None:
            self.dust.S.hyd = Field(self, np.zeros(
                shape2), description="Hydrodynamic sources [g/cm²/s]")
            self.dust.S.hyd.updater = std_dust.S_hyd
        if self.dust.S.tot is None:
            self.dust.S.tot = Field(self, np.zeros(
                shape2), description="Tot sources [g/cm²/s]")
            self.dust.S.tot.updater = std_dust.S_tot
        # Stokes number
        if self.dust.St is None:
            self.dust.St = Field(self, np.zeros(
                shape2), description="Stokes number")
            self.dust.St.updater = std_dust.St_Epstein_StokesI
        # Velocities
        if self.dust.v.frag is None:
            vfrag = self.ini.dust.vfrag * np.ones(shape2)
            self.dust.v.frag = Field(
                self, vfrag, description="Fragmentation velocity [cm/s]")
        if self.dust.v.rel.azi is None:
            self.dust.v.rel.azi = Field(self, np.zeros(
                shape3), description="Relative azimuthal velocity [cm/s]")
            self.dust.v.rel.azi.updater = std_dust.vrel_azimuthal_drift
        if self.dust.v.rel.brown is None:
            self.dust.v.rel.brown = Field(self, np.zeros(
                shape3), description="Relative Brownian motion velocity [cm/s]")
            self.dust.v.rel.brown.updater = std_dust.vrel_brownian_motion
        if self.dust.v.rel.rad is None:
            self.dust.v.rel.rad = Field(self, np.zeros(
                shape3), description="Relative radial velocity [cm/s]")
            self.dust.v.rel.rad.updater = std_dust.vrel_radial_drift
        if self.dust.v.rel.turb is None:
            self.dust.v.rel.turb = Field(self, np.zeros(
                shape3), description="Relative turbulent velocity [cm/s]")
            self.dust.v.rel.turb.updater = std_dust.vrel_turbulent_motion
        if self.dust.v.rel.vert is None:
            self.dust.v.rel.vert = Field(self, np.zeros(
                shape3), description="Relative vertical settling velocity [cm/s]")
            self.dust.v.rel.vert.updater = std_dust.vrel_vertical_settling
        if self.dust.v.rel.tot is None:
            self.dust.v.rel.tot = Field(self, np.zeros(
                shape3), description="Total relative velocity [cm/s]")
            self.dust.v.rel.tot.updater = std_dust.vrel_tot
        if self.dust.v.driftmax is None:
            self.dust.v.driftmax = Field(self, np.zeros(
                shape1), description="Maximum drift velocity [cm/s]")
            self.dust.v.driftmax.updater = std_dust.vdriftmax
        if self.dust.v.rad is None:
            self.dust.v.rad = Field(self, np.zeros(
                shape2), description="Radial velocity [cm/s]")
            self.dust.v.rad.updater = std_dust.vrad
        # Initialize dust quantities partly to calculate Sigma
        try:
            self.dust.update()
        except:
            pass
        # Floor value
        if self.dust.SigmaFloor is None:
            SigmaFloor = std_dust.SigmaFloor(self)
            self.dust.SigmaFloor = Field(
                self, SigmaFloor, description="Floor value of surface density [g/cm²]")
        # Surface density, if not set
        if self.dust.Sigma is None:
            Sigma = np.maximum(std_dust.MRN_distribution(
                self), self.dust.SigmaFloor)
            self.dust.Sigma = Field(
                self, Sigma, description="Surface density per mass bin [g/cm²]")
            self.dust.Sigma.differentiator = std_dust.Sigma_deriv
        # Fully initialize dust quantities
        self.dust.update()

        # INTEGRATION VARIABLE

        if self.t is None:
            self.t = IntVar(self, 0., description="Time [s")
            self.t.updater = std_sim.dt
            self.t.snapshots = np.logspace(3., 4., num=2, base=10.) * c.year
            self.t.suggest(1.*c.year)

        # INTEGRATOR

        if self.integrator is None:
            instructions = [
                Instruction(schemes.expl_2_heun_euler_adptv,
                            self.gas.Sigma,
                            controller={"dYdx": self.gas.S.tot,
                                        "eps": 0.01
                                        }
                            ),
                Instruction(schemes.expl_5_cash_karp_adptv,
                            self.dust.Sigma,
                            controller={"dYdx": self.dust.S.tot,
                                        "eps": 0.01
                                        }
                            )
            ]
            self.integrator = Integrator(self.t)
            self.integrator.instructions = instructions
            self.integrator.finalization = std_sim.finalize

        # WRITER

        if self.writer is None:
            self.writer = hdf5writer
