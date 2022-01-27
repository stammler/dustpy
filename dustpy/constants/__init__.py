'''Package containing constants that are used in the simulations. All constants
are in cgs units.

The constants are defined in ``dustpy/constants/constants.f90`` so they can be used in
the ``Fortran`` scripts.

Constants
---------
| ``au`` : Astronomical unit
| ``G`` : Gravitational constant
| ``k_B`` : Boltzmann constant
| ``m_p`` : Proton mass
| ``M_earth`` : Mass of the Earth
| ``M_jup`` : Mass of Jupiter
| ``M_sun`` : Mass of the Sun
| ``pi`` : Ratio of circle's circumference to diameter
| ``R_sun`` : Radius of the sun
| ``sigma_H2`` : Geometrical cross-section of H2 molecule
| ``sigma_sb`` : Stephan-Boltzmann constant
| ``year`` : Year in seconds
'''

import numpy as _np

import dustpy.constants._constants_f as _c

au = float(_c.constants.au)
G = float(_c.constants.g)
k_B = float(_c.constants.k_b)
m_p = float(_c.constants.m_p)
M_earth = float(_c.constants.m_earth)
M_jup = float(_c.constants.m_jup)
M_sun = float(_c.constants.m_sun)
pi = float(_c.constants.pi)
R_sun = float(_c.constants.r_sun)
sigma_H2 = float(_c.constants.sigma_h2)
sigma_sb = float(_c.constants.sigma_sb)
year = float(_c.constants.year)

__all__ = ["au",
           "G",
           "k_B",
           "m_p",
           "M_earth",
           "M_jup",
           "M_sun",
           "pi",
           "R_sun",
           "sigma_H2",
           "sigma_sb",
           "year"]
