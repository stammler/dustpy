from simframe import Frame
from simframe import Instruction
from simframe import Integrator
from simframe import schemes
from simframe.frame import Field
from simframe.frame import Group
from simframe.frame import IntVar
from simframe.utils.color import colorize

import dustpy.constants as c

from dustpy import std

from dustpy.utils import hdf5writer
from dustpy.utils.boundary import Boundary

import numpy as np
from types import SimpleNamespace


class Simulation(Frame):
    """The main simulation class for running dust coagulation simulations.

    `dustpy.Simulation`` is a child of ``simframe.Frame``.

    For setting simple initial conditions use ``Simulation.ini``,
    For making the simulation grids use ``Simulation.makegrids()``,
    For initialization use ``Simulation.initialize()``,
    For running simulations use ``Simulation.run()``.

    Please have a look at the documentation of ``simframe`` for further details."""

    __name__ = "DustPy"

    def __init__(self, **kwargs):
        """Main simulation class."""

        super().__init__(**kwargs)

        # Namespace with parameters to set the initial conditions
        self._ini = SimpleNamespace(**{"dust": SimpleNamespace(**{"aIniMax": 0.0001,
                                                                  "allowDriftingParticles": False,
                                                                  "erosionMassRatio": 10.,
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
                                                                 "Mdisk": 0.05*c.M_sun,
                                                                 "mu": 2.3*c.m_p,
                                                                 "SigmaExp": -1.,
                                                                 "SigmaRc": 60.*c.au
                                                                 }
                                                              ),
                                       "grid": SimpleNamespace(**{"Nmbpd": 7,
                                                                  "mmin": 1.e-12,
                                                                  "mmax": 1.e5,
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

        # Initializing everything with None to get the basic frame

        # Dust quantities
        self.dust = Group(self, description="Dust quantities")
        self.dust.a = None
        self.dust.backreaction = Group(
            self, description="Backreaction coefficients")
        self.dust.backreaction.A = None
        self.dust.backreaction.B = None
        self.dust.backreaction.updater = ["A", "B"]
        self.dust.boundary = Group(self, description="Boundary conditions")
        self.dust.boundary.inner = None
        self.dust.boundary.outer = None
        self.dust.coagulation = Group(
            self, description="Coagulation quantities")
        self.dust.coagulation.stick = None
        self.dust.coagulation.stick_ind = None
        self.dust.coagulation.A = None
        self.dust.coagulation.eps = None
        self.dust.coagulation.lf_ind = None
        self.dust.coagulation.rm_ind = None
        self.dust.coagulation.phi = None
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
        self.dust.kernel = None
        self.dust.rho = None
        self.dust.rhos = None
        self.dust.p = Group(self, description="Probabilities")
        self.dust.p.frag = None
        self.dust.p.stick = None
        self.dust.p.updater = ["frag", "stick"]
        self.dust.S = Group(self, description="Sources")
        self.dust.S.coag = None
        self.dust.S.ext = None
        self.dust.S.hyd = None
        self.dust.S.tot = None
        self.dust.S.updater = ["ext", "tot"]
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
        self.dust.v.updater = ["frag", "driftmax", "rel"]
        self.dust.updater = ["delta", "rhos", "fill", "a", "St", "H",
                             "rho", "backreaction", "v", "D", "eps", "kernel", "p", "S"]

        # Gas quantities
        self.gas = Group(self, description="Gas quantities")
        self.gas.alpha = None
        self.gas.boundary = Group(self, description="Boundary conditions")
        self.gas.boundary.inner = None
        self.gas.boundary.outer = None
        self.gas.cs = None
        self.gas.eta = None
        self.gas.Fi = None
        self.gas.gamma = None
        self.gas.Hp = None
        self.gas.mfp = None
        self.gas.mu = None
        self.gas.n = None
        self.gas.nu = None
        self.gas.P = None
        self.gas.rho = None
        self.gas.S = Group(self, description="Source terms")
        self.gas.S.ext = None
        self.gas.S.hyd = None
        self.gas.S.tot = None
        self.gas.S.updater = ["ext", "tot"]
        self.gas.Sigma = None
        self.gas.SigmaFloor = None
        self.gas.T = None
        self.gas.v = Group(self, description="Velocities")
        self.gas.v.rad = None
        self.gas.v.visc = None
        self.gas.v.updater = ["visc", "rad"]
        self.gas.updater = ["gamma", "mu", "T", "alpha", "cs", "Hp", "nu",
                            "rho", "n", "mfp", "P", "eta", "S"]

        # Grid quantities
        self.grid = Group(self, description="Grid quantities")
        self.grid.Nr = None
        self.grid.r = None
        self.grid.ri = None
        self.grid.A = None
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

    def run(self):
        """This functions runs the simulation."""
        # Print welcome message
        if self.verbosity > 0:
            msg = ""
            msg += "\nDustPy v{}".format(self.__version__)
            msg += "\n"
            msg += "\nDocumentation: {}".format(
                "https://stammler.github.io/dustpy/")
            msg += "\nPyPI:          {}".format(
                "https://pypi.org/project/dustpy/")
            msg += "\nGitHub:        {}".format(
                "https://github.com/stammler/dustpy/")
            msg += "\n"
            msg += colorize("\nPlease cite Stammler & Birnstiel (2022).", "blue")
            print(msg)
        # Check for mass conserbation
        self.checkmassconservation()
        # Actually run the simulation
        super().run()

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
        logmmin = np.log10(self.ini.grid.mmin)
        logmmax = np.log10(self.ini.grid.mmax)
        decades = np.ceil(logmmax - logmmin)
        Nm = int(decades * self.ini.grid.Nmbpd) + 1
        m = np.logspace(np.log10(self.ini.grid.mmin), np.log10(
            self.ini.grid.mmax), num=Nm, base=10.)
        self.grid.Nm = Field(
            self, Nm, description="# of mass bins", constant=True)
        self.grid.m = Field(
            self, m, description="Mass grid [g]", constant=True)

    def _makeradialgrid(self):
        '''Function sets the mass grid using the parameters set in ``Simulation.ini``.'''
        if self.grid.ri is None:
            ri = np.logspace(np.log10(self.ini.grid.rmin), np.log10(
                self.ini.grid.rmax), num=self.ini.grid.Nr+1, base=10.)
            Nr = self.ini.grid.Nr
        else:
            ri = self.grid.ri
            Nr = ri.shape[0] - 1
        r = 0.5*(ri[:-1] + ri[1:])
        A = np.pi*(ri[1:]**2 - ri[:-1]**2)
        self.grid.Nr = Field(
            self, Nr, description="# of radial grid cells", constant=True)
        self.grid.r = Field(
            self, r, description="Radial grid cell centers [cm]", constant=True)
        self.grid.ri = Field(
            self, ri, description="Radial grid cell interfaces [cm]", constant=True)
        self.grid.A = Field(
            self, A, description="Radial grid annulus area [cm²]", constant=True)

    def checkmassconservation(self):
        """Function checks for mass conservation and prints the maximum relative mass error."""

        # Check if required fields are present
        if self.dust.coagulation.stick is None:
            raise RuntimeError(
                "'Simulation.dust.coagulation.stick' is not set.")
        if self.dust.coagulation.stick_ind is None:
            raise RuntimeError(
                "'Simulation.dust.coagulation.stick_ind' is not set.")
        if self.dust.coagulation.A is None:
            raise RuntimeError(
                "'Simulation.dust.coagulation.A' is not set.")
        if self.dust.coagulation.eps is None:
            raise RuntimeError(
                "'Simulation.dust.coagulation.eps' is not set.")
        if self.dust.coagulation.lf_ind is None:
            raise RuntimeError(
                "'Simulation.dust.coagulation.lf_ind' is not set.")
        if self.dust.coagulation.rm_ind is None:
            raise RuntimeError(
                "'Simulation.dust.coagulation.rm_ind' is not set.")
        if self.dust.coagulation.phi is None:
            raise RuntimeError(
                "'Simulation.dust.coagulation.phi' is not set.")
        if self.grid.m is None:
            raise RuntimeError("'sim.grid.m' is not set.")

        if self.verbosity > 0:

            # Maximum acceptable error
            erracc = 1.e-13

            # Checking for sticking error
            msg = "\n"
            msg += colorize("Checking for mass conservation...\n",
                            color="yellow")
            print(msg)
            msg = colorize("    - Sticking:", color="yellow")
            print(msg)
            errmax, i, j = std.dust_f.check_mass_conservation_sticking(
                self.dust.coagulation.stick, self.dust.coagulation.stick_ind, self.grid.m)
            tup = (j, i)
            color = "red"
            if(errmax < erracc):
                color = "green"
            error = "{:9.2e}".format(errmax)
            msg = "        max. rel. error: {:}\n".format(
                colorize(error, color=color))
            msg += "        for particle collision\n"
            msg += "            m[{:d}] = {:9.2e} g    with\n".format(
                tup[0], self.grid.m[tup[0]])
            msg += "            m[{:d}] = {:9.2e} g".format(
                tup[1], self.grid.m[tup[1]])
            msg = colorize(msg)
            print(msg)

            # Checking for full fragmentation error
            msg = colorize("    - Full fragmentation:", color="yellow")
            print(msg)
            A = self.dust.coagulation.A
            eps = self.dust.coagulation.eps
            klf = self.dust.coagulation.lf_ind
            krm = self.dust.coagulation.rm_ind
            m = self.grid.m
            phi = self.dust.coagulation.phi
            errmax, i, j = std.dust_f.check_mass_conservation_full_fragmentation(
                A, klf, m, phi)
            tup = (j, i)
            color = "red"
            if(errmax < erracc):
                color = "green"
            error = "{:9.2e}".format(errmax)
            msg = "        max. rel. error: {:}\n".format(
                colorize(error, color=color))
            msg += "        for particle collision\n"
            msg += "            m[{:d}] = {:9.2e} g    with\n".format(
                tup[0], self.grid.m[tup[0]])
            msg += "            m[{:d}] = {:9.2e} g".format(
                tup[1], self.grid.m[tup[1]])
            msg = colorize(msg)
            print(msg)

            # Checking for erosion error
            msg = colorize("    - Erosion:", color="yellow")
            print(msg)
            errmax, i, j = std.dust_f.check_mass_conservation_erosion(
                A, eps, klf, krm, m, phi)
            tup = (j, i)
            color = "red"
            if(errmax < erracc):
                color = "green"
            error = "{:9.2e}".format(errmax)
            msg = "        max. rel. error: {:}\n".format(
                colorize(error, color=color))
            msg += "        for particle collision\n"
            msg += "            m[{:d}] = {:9.2e} g    with\n".format(
                tup[0], self.grid.m[tup[0]])
            msg += "            m[{:d}] = {:9.2e} g\n".format(
                tup[1], self.grid.m[tup[1]])
            msg = colorize(msg)
            print(msg)

    def initialize(self):
        '''Function initializes the simulation frame.

        Function sets all fields that are None with a standard value.
        If the grids are not set, it will call ``Simulation.makegrids()`` first.'''
        if not isinstance(self.grid.Nm, Field) or not isinstance(self.grid.Nr, Field):
            self.makegrids()

        # INTEGRATION VARIABLE
        if self.t is None:
            self.t = IntVar(self, 0., description="Time [s]")
            self.t.cfl = 0.1
            self.t.updater = std.sim.dt
            self.t.snapshots = np.logspace(3., 5., num=21, base=10.) * c.year

        # STELLAR QUANTITIES
        self._initializestar()

        # GRID QUANTITIES
        self._initializegrid()

        # GAS QUANTITIES
        self._initializegas()

        # DUST QUANTITIES
        self._initializedust()

        # INTEGRATOR
        if self.integrator is None:
            instructions = [
                Instruction(std.dust.impl_1_direct,
                            self.dust.Sigma,
                            controller={"rhs": self.dust._rhs
                                        },
                            description="Dust: implicit 1st-order direct solver"
                            ),
                Instruction(std.gas.impl_1_direct,
                            self.gas.Sigma,
                            controller={"rhs": self.gas._rhs
                                        },
                            description="Gas: implicit 1st-order direct solver"
                            ),
            ]
            self.integrator = Integrator(
                self.t, description="Default integrator")
            self.integrator.instructions = instructions
            self.integrator.preparator = std.sim.prepare_implicit_dust
            self.integrator.finalizer = std.sim.finalize_implicit_dust

        # WRITER
        if self.writer is None:
            self.writer = hdf5writer()

        # Updating the entire Simulation object including integrator finalization
        self.integrator._finalize()
        self.update()

    def _initializedust(self):
        '''Function to initialize dust quantities'''
        shape1 = (int(self.grid.Nr))
        shape2 = (int(self.grid.Nr), int(self.grid.Nm))
        shape2ravel = (int(self.grid.Nr*self.grid.Nm))
        shape2p1 = (int(self.grid.Nr)+1, int(self.grid.Nm))
        shape3 = (int(self.grid.Nr), int(
            self.grid.Nm), int(self.grid.Nm))
        # Particle size
        if self.dust.a is None:
            self.dust.a = Field(self, np.zeros(shape2),
                                description="Particle size [cm]")
            self.dust.a.updater = std.dust.a
        # Coagulation parameters
        stick, stick_ind, A, eps, lf_ind, rm_ind, phi = std.dust.coagulation_parameters(
            self)
        if self.dust.coagulation.A is None:
            self.dust.coagulation.A = Field(
                self, A, description="Fragment normalization factors", constant=True)
        if self.dust.coagulation.eps is None:
            self.dust.coagulation.eps = Field(
                self, eps, description="Remnant mass distribution", constant=True)
        if self.dust.coagulation.lf_ind is None:
            self.dust.coagulation.lf_ind = Field(
                self, lf_ind, description="Index of largest fragment", constant=True)
        if self.dust.coagulation.rm_ind is None:
            self.dust.coagulation.rm_ind = Field(
                self, rm_ind, description="Smaller index of remnant", constant=True)
        if self.dust.coagulation.phi is None:
            self.dust.coagulation.phi = Field(
                self, phi, description="Fragment distribution", constant=True)
        if self.dust.coagulation.stick is None:
            self.dust.coagulation.stick = Field(
                self, stick, description="Sticking matrix", constant=True)
        if self.dust.coagulation.stick_ind is None:
            self.dust.coagulation.stick_ind = Field(
                self, stick_ind, description="Non-zero elements of sticking matrix", constant=True)
        # Diffusivity
        if self.dust.D is None:
            self.dust.D = Field(self, np.zeros(shape2),
                                description="Diffusivity [cm²/s]")
            self.dust.D.updater = std.dust.D
        # Deltas
        if self.dust.delta.rad is None:
            delta = self.ini.gas.alpha * np.ones(shape1)
            self.dust.delta.rad = Field(
                self, delta, description="Radial mixing parameter")
        if self.dust.delta.turb is None:
            delta = self.ini.gas.alpha * np.ones(shape1)
            self.dust.delta.turb = Field(
                self, delta, description="Turbulent mixing parameter")
        if self.dust.delta.vert is None:
            delta = self.ini.gas.alpha * np.ones(shape1)
            self.dust.delta.vert = Field(
                self, delta, description="Vertical mixing parameter")
        # Vertically integrated dust to gas ratio
        if self.dust.eps is None:
            self.dust.eps = Field(self, np.zeros(
                shape1), description="Dust-to-gas ratio")
            self.dust.eps.updater = std.dust.eps
        # Fluxes
        if self.dust.Fi.adv is None:
            self.dust.Fi.adv = Field(self, np.zeros(
                shape2p1), description="Advective flux [g/cm/s]")
            self.dust.Fi.adv.updater = std.dust.F_adv
        if self.dust.Fi.diff is None:
            self.dust.Fi.diff = Field(self, np.zeros(
                shape2p1), description="Diffusive flux [g/cm/s]")
            self.dust.Fi.diff.updater = std.dust.F_diff
        if self.dust.Fi.tot is None:
            self.dust.Fi.tot = Field(self, np.zeros(
                shape2p1), description="Total flux [g/cm/s]")
            self.dust.Fi.tot.updater = std.dust.F_tot
        # Filling factor
        if self.dust.fill is None:
            self.dust.fill = Field(self, np.ones(
                shape2), description="Filling factor")
        # Scale height
        if self.dust.H is None:
            self.dust.H = Field(self, np.zeros(shape2),
                                description="Scale heights [cm]")
            self.dust.H.updater = std.dust.H
        # Kernel
        if self.dust.kernel is None:
            self.dust.kernel = Field(self, np.zeros(
                shape3), description="Collision kernel [cm²/s]")
            self.dust.kernel.updater = std.dust.kernel
        # Midplane mass density
        if self.dust.rho is None:
            self.dust.rho = Field(self, np.zeros(
                shape2), description="Midplane mass density per mass bin [g/cm³]")
            self.dust.rho.updater = std.dust.rho_midplane
        # Solid state density
        if self.dust.rhos is None:
            rhos = self.ini.dust.rhoMonomer * np.ones(shape2)
            self.dust.rhos = Field(
                self, rhos, description="Solid state density [g/cm³]")
        # Probabilities
        if self.dust.p.frag is None:
            self.dust.p.frag = Field(self, np.zeros(
                shape3), description="Fragmentation probability")
            self.dust.p.frag.updater = std.dust.p_frag
        if self.dust.p.stick is None:
            self.dust.p.stick = Field(self, np.zeros(
                shape3), description="Sticking probability")
            self.dust.p.stick.updater = std.dust.p_stick
        # Source terms
        if self.dust.S.coag is None:
            self.dust.S.coag = Field(self, np.zeros(
                shape2), description="Coagulation sources [g/cm²/s]")
            self.dust.S.coag.updater = std.dust.S_coag
        if self.dust.S.ext is None:
            self.dust.S.ext = Field(self, np.zeros(
                shape2), description="External sources [g/cm²/s]")
        if self.dust.S.hyd is None:
            self.dust.S.hyd = Field(self, np.zeros(
                shape2), description="Hydrodynamic sources [g/cm²/s]")
            self.dust.S.hyd.updater = std.dust.S_hyd
        if self.dust.S.tot is None:
            self.dust.S.tot = Field(self, np.zeros(
                shape2), description="Tot sources [g/cm²/s]")
            self.dust.S.tot.updater = std.dust.S_tot
        # Stokes number
        if self.dust.St is None:
            self.dust.St = Field(self, np.zeros(
                shape2), description="Stokes number")
            self.dust.St.updater = std.dust.St_Epstein_StokesI
        # Velocities
        if self.dust.v.frag is None:
            vfrag = self.ini.dust.vfrag * np.ones(shape1)
            self.dust.v.frag = Field(
                self, vfrag, description="Fragmentation velocity [cm/s]")
        if self.dust.v.rel.azi is None:
            self.dust.v.rel.azi = Field(self, np.zeros(
                shape3), description="Relative azimuthal velocity [cm/s]")
            self.dust.v.rel.azi.updater = std.dust.vrel_azimuthal_drift
        if self.dust.v.rel.brown is None:
            self.dust.v.rel.brown = Field(self, np.zeros(
                shape3), description="Relative Brownian motion velocity [cm/s]")
            self.dust.v.rel.brown.updater = std.dust.vrel_brownian_motion
        if self.dust.v.rel.rad is None:
            self.dust.v.rel.rad = Field(self, np.zeros(
                shape3), description="Relative radial velocity [cm/s]")
            self.dust.v.rel.rad.updater = std.dust.vrel_radial_drift
        if self.dust.v.rel.turb is None:
            self.dust.v.rel.turb = Field(self, np.zeros(
                shape3), description="Relative turbulent velocity [cm/s]")
            self.dust.v.rel.turb.updater = std.dust.vrel_turbulent_motion
        if self.dust.v.rel.vert is None:
            self.dust.v.rel.vert = Field(self, np.zeros(
                shape3), description="Relative vertical settling velocity [cm/s]")
            self.dust.v.rel.vert.updater = std.dust.vrel_vertical_settling
        if self.dust.v.rel.tot is None:
            self.dust.v.rel.tot = Field(self, np.zeros(
                shape3), description="Total relative velocity [cm/s]")
            self.dust.v.rel.tot.updater = std.dust.vrel_tot
        if self.dust.v.driftmax is None:
            self.dust.v.driftmax = Field(self, np.zeros(
                shape1), description="Maximum drift velocity [cm/s]")
            self.dust.v.driftmax.updater = std.dust.vdriftmax
        if self.dust.v.rad is None:
            self.dust.v.rad = Field(self, np.zeros(
                shape2), description="Radial velocity [cm/s]")
            self.dust.v.rad.updater = std.dust.vrad
        # Initialize dust quantities partly to calculate Sigma
        try:
            self.dust.update()
        except:
            pass
        # Floor value
        if self.dust.SigmaFloor is None:
            SigmaFloor = 0.1*std.dust.SigmaFloor(self)
            self.dust.SigmaFloor = Field(
                self, SigmaFloor, description="Floor value of surface density [g/cm²]")
        # Surface density, if not set
        if self.dust.Sigma is None:
            Sigma = std.dust.MRN_distribution(self)
            Sigma = np.where(Sigma <= self.dust.SigmaFloor,
                             0.1*self.dust.SigmaFloor,
                             Sigma)
            self.dust.Sigma = Field(
                self, Sigma, description="Surface density per mass bin [g/cm²]")
        self.dust.Sigma.differentiator = std.dust.Sigma_deriv
        self.dust.Sigma.jacobinator = std.dust.jacobian
        # Fully initialize dust quantities
        self.dust.update()
        # Hidden fields
        # We store the old values of the surface density in a hidden field
        # to calculate the fluxes through the boundaries in case of implicit integration.
        self.dust._SigmaOld = Field(
            self, self.dust.Sigma, description="Previous value of surface density [g/cm²]")
        # The right-hand side of the matrix equation is stored in a hidden field
        self.dust._rhs = Field(self, np.zeros(
            shape2ravel), description="Right-hand side of matrix equation [g/cm²]")
        # Boundary conditions
        if self.dust.boundary.inner is None:
            self.dust.boundary.inner = Boundary(
                self.grid.r,
                self.grid.ri,
                self.dust.Sigma,
                condition="const_grad"
            )
        if self.dust.boundary.outer is None:
            self.dust.boundary.outer = Boundary(
                self.grid.r[::-1],
                self.grid.ri[::-1],
                self.dust.Sigma[::-1],
                condition="val",
                value=0.1*self.dust.SigmaFloor[-1]
            )

    def _initializegas(self):
        '''Function to initialize gas quantities'''
        shape1 = (int(self.grid.Nr))
        shape1p1 = (int(self.grid.Nr)+1)
        # Turbulent alpha parameter
        if self.gas.alpha is None:
            alpha = self.ini.gas.alpha * np.ones(shape1)
            self.gas.alpha = Field(
                self, alpha, description="Turbulent alpha parameter")
        # Sound speed
        if self.gas.cs is None:
            self.gas.cs = Field(self, np.zeros(shape1),
                                description="Sound speed [cm/s]")
            self.gas.cs.updater = std.gas.cs_adiabatic
        # Pressure gradient parameter
        if self.gas.eta is None:
            self.gas.eta = Field(self, np.zeros(
                shape1), description="Pressure gradient parameter")
            self.gas.eta.updater = std.gas.eta_midplane
        # Gas flux at the cell interfaces
        if self.gas.Fi is None:
            self.gas.Fi = Field(self, np.zeros(shape1p1),
                                description="Gas flux interfaces [g/cm/s]")
            self.gas.Fi.updater = std.gas.Fi
        # Adiabatic index
        if self.gas.gamma is None:
            gamma = self.ini.gas.gamma * np.ones(shape1)
            self.gas.gamma = Field(self, gamma,
                                   description="Adiabatic index")
        # Pressure scale height
        if self.gas.Hp is None:
            self.gas.Hp = Field(self, np.zeros(shape1),
                                description="Pressure scale height [cm]")
            self.gas.Hp.updater = std.gas.Hp
        # Mean free path
        if self.gas.mfp is None:
            self.gas.mfp = Field(self, np.zeros(shape1),
                                 description="Midplane mean free path [cm]")
            self.gas.mfp.updater = std.gas.mfp_midplane
        # Mean molecular weight
        if self.gas.mu is None:
            mu = self.ini.gas.mu * np.ones(shape1)
            self.gas.mu = Field(self, mu,
                                description="Mean molecular weight [g]")
        # Midplane number density
        if self.gas.n is None:
            self.gas.n = Field(self, np.zeros(shape1),
                               description="Miplane number density [1/cm³]")
            self.gas.n.updater = std.gas.n_midplane
        # Viscosity
        if self.gas.nu is None:
            self.gas.nu = Field(self, np.zeros(shape1),
                                description="Kinematic viscosity [cm²/s]")
            self.gas.nu.updater = std.gas.nu
        # Midplane pressure
        if self.gas.P is None:
            self.gas.P = Field(self, np.zeros(shape1),
                               description="Midplane pressure [g/cm/s²]")
            self.gas.P.updater = std.gas.P_midplane
        # Midplane mass density
        if self.gas.rho is None:
            self.gas.rho = Field(self, np.zeros(shape1),
                                 description="Miplane mass density [g/cm³]")
            self.gas.rho.updater = std.gas.rho_midplane
        # Sources
        if self.gas.S.hyd is None:
            self.gas.S.hyd = Field(self, np.zeros(
                shape1), description="Hydrodynamic sources [g/cm²/s]")
            self.gas.S.hyd.updater = std.gas.S_hyd
        if self.gas.S.ext is None:
            self.gas.S.ext = Field(self, np.zeros(
                shape1), description="External sources [g/cm²/s]")
        if self.gas.S.tot is None:
            self.gas.S.tot = Field(self, np.zeros(
                shape1), description="Total sources [g/cm²/s]")
            self.gas.S.tot.updater = std.gas.S_tot
        # Floor value
        if self.gas.SigmaFloor is None:
            self.gas.SigmaFloor = Field(
                self, 1.e-100*np.ones(shape1), description="Floor value of surface density [g/cm²]")
        # Surface density
        if self.gas.Sigma is None:
            SigmaGas = np.array(std.gas.lyndenbellpringle1974(
                self.grid.r, self.ini.gas.SigmaRc, self.ini.gas.SigmaExp, self.ini.gas.Mdisk))
            SigmaGas = np.maximum(SigmaGas, self.gas.SigmaFloor)
            self.gas.Sigma = Field(self, SigmaGas,
                                   description="Surface density [g/cm²]")
        self.gas.Sigma.jacobinator = std.gas.jacobian
        # Temperature
        if self.gas.T is None:
            self.gas.T = Field(self, np.zeros(shape1),
                               description="Temperature [K]")
            self.gas.T.updater = std.gas.T_passive
        # Velocities
        # Viscous accretion velocity
        if self.gas.v.visc is None:
            self.gas.v.visc = Field(self, np.zeros(shape1),
                                    description="Viscous accretion velocity [cm/s]")
            self.gas.v.visc.updater = std.gas.vvisc
        # Radial gas velocity
        if self.gas.v.rad is None:
            self.gas.v.rad = Field(self, np.zeros(shape1),
                                   description="Radial velocity [cm/s]")
            self.gas.v.rad.updater = std.gas.vrad
        # Hidden fields
        # We store the old values of the surface density in a hidden field
        # to calculate the fluxes through the boundaries.
        self.gas._SigmaOld = Field(self, np.zeros(
            shape1), description="Previous value of surface density [g/cm²]")
        self.gas._SigmaOld[:] = self.gas.Sigma
        # The right-hand side of the matrix equation is stored in a hidden field
        self.gas._rhs = Field(self, np.zeros(
            shape1), description="Right-hand side of matrix equation [g/cm²]")

        # The dust backreaction coefficients have to be initialized before the gas,
        # since the gas velocities need them.
        # Backreaction
        if self.dust.backreaction.A is None:
            self.dust.backreaction.A = Field(
                self, np.ones(shape1), description="Pull factor")
        if self.dust.backreaction.B is None:
            self.dust.backreaction.B = Field(
                self, np.zeros(shape1), description="Push factor")

        # Initialize gas quantities
        self.gas.update()
        # Boundary conditions
        if self.gas.boundary.inner is None:
            self.gas.boundary.inner = Boundary(
                self.grid.r, self.grid.ri, self.gas.Sigma)
            self.gas.boundary.inner.setcondition("const_grad")
        if self.gas.boundary.outer is None:
            self.gas.boundary.outer = Boundary(
                self.grid.r[::-1], self.grid.ri[::-1], self.gas.Sigma[::-1])
            self.gas.boundary.outer.setcondition(
                "val", 0.1*self.gas.SigmaFloor[-1])

    def _initializegrid(self):
        '''Function to initialize grid quantities'''
        shape1 = (int(self.grid.Nr))
        # Keplerian frequency
        if self.grid.OmegaK is None:
            self.grid.OmegaK = Field(self, np.zeros(shape1),
                                     description="Keplerian frequency [1/s]")
            self.grid.OmegaK.updater = std.grid.OmegaK
        # Initialize grid quantities
        self.grid.update()

    def _initializestar(self):
        '''Function to initialize the stellar quantities'''
        # Luminosity
        if self.star.L is None:
            self.star.L = Field(self, 0., description="Luminosity [erg/s]")
            self.star.L.updater = std.star.luminosity
        # Mass
        if self.star.M is None:
            self.star.M = Field(self, self.ini.star.M,
                                description="Mass [g]")
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

    def setdustintegrator(self, scheme="explicit", method="cash-karp"):
        """Function sets the dust integrator.

        Parameters
        ----------
        scheme : string, optional, default : "explicit"
            Possible values
                {"explicit", "implicit"}
        method : string, optional, default : "cash-karp"
            Possible values for explicit integration
                {"cash-karp"}
            Possible values for implicit integration
                {"direct", "gmres", "bicgstab}"""

        if not isinstance(self.grid.Nm, Field) or not isinstance(self.grid.Nr, Field):
            raise RuntimeError(
                "The simulation frame has to be initialized before calling setdustimplicit().")

        # Get index of dust instruction
        for i, inst in enumerate(self.integrator.instructions):
            if inst.Y is self.dust.Sigma:
                break

        if scheme == "implicit":

            shape2ravel = (int(self.grid.Nr*self.grid.Nm))

            # Hidden fields
            # We store the old values of the surface density in a hidden field
            # to calculate the fluxes through the boundaries in case of implicit integration.
            self.dust._SigmaOld = Field(
                self, self.dust.Sigma, description="Previous value of surface density [g/cm²]")
            # The right-hand side of the matrix equation is stored in a hidden field
            self.dust._rhs = Field(self, np.zeros(
                shape2ravel), description="Right-hand side of matrix equation [g/cm²]")

            # Setting the Jacobinator
            self.dust.Sigma.jacobinator = std.dust.jacobian

            # Time step routine
            self.t.updater = std.sim.dt

            # Updaters
            self.dust.v.updater = ["frag", "driftmax", "rel"]
            self.dust.updater = ["delta", "rhos", "fill", "a", "St", "H",
                                 "rho", "backreaction", "v", "D", "eps", "kernel", "p", "S"]
            self.dust.S.updater = ["ext", "tot"]

            # Preparation/Finalization
            self.integrator.preparator = std.sim.prepare_implicit_dust
            self.integrator.finalizer = std.sim.finalize_implicit_dust

            # Integrator
            if method == "direct":

                inst = Instruction(std.dust.impl_1_direct,
                                   self.dust.Sigma,
                                   controller={"rhs": self.dust._rhs
                                               },
                                   description="Dust: implicit 1st-order direct solver"
                                   )
                self.integrator.instructions[i] = inst

            elif method == "gmres":
                raise NotImplementedError(
                    "GMRES method is not implemented, yet.")
            elif method == "bicgstab":
                raise NotImplementedError("BiCGSTAB is not implemented, yet.")
            else:
                raise RuntimeError("Invalid method for implicit integration.")
        elif scheme == "explicit":

            # Remove hidden fields if they exist
            if hasattr(self.dust, "_SigmaOld"):
                del self.dust._SigmaOld
            if hasattr(self.dust, "_rhs"):
                del self.dust._rhs

            # Unset Jacobian
            self.dust.Sigma.jacobinator = None

            # Updaters
            self.dust.v.updater = ["frag", "driftmax", "rad", "rel"]
            self.dust.updater = ["delta", "rhos", "fill", "a", "St", "H",
                                 "rho", "backreaction", "v", "D", "eps", "Fi", "kernel", "p", "S"]
            self.dust.S.updater = ["coag", "hyd", "ext", "tot"]

            # Preparation/Finalization
            self.integrator.preparator = std.sim.prepare_explicit_dust
            self.integrator.finalizer = std.sim.finalize_explicit_dust

            if method == "cash-karp":

                # Adaptive time step routine
                self.t.updater = std.sim.dt_adaptive

                # Instruction
                inst = Instruction(schemes.expl_5_cash_karp_adptv,
                                   self.dust.Sigma,
                                   controller={"dYdx": self.dust.S.tot,
                                               "eps": 0.1,
                                               "S": 0.9,
                                               },
                                   description="Dust: explicit 5th-order adaptive Cash-Karp method"
                                   )
                self.integrator.instructions[i] = inst
                self.t.suggest(1.*c.year)

            else:
                raise RuntimeError("Invalid method for explicit integration.")
        else:
            raise RuntimeError("Unknown integration scheme.")

        self.integrator._finalize()
        self.update()

        if self.verbosity > 0:
            msg = "Setting dust integrator\n    scheme: {}\n    method: {}".format(
                colorize(scheme, "blue"), colorize(method, "blue"))
            print(msg)
