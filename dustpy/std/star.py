'''Module containing standard functions for the central star.'''

import dustpy.constants as c


def luminosity(sim):
    """Function calculates the luminosity of the star.

    Parameters
    ----------
    sim : Simulation
        Parent simulation frame

    Returns
    -------
    L : Field
        Luminostiy of star"""
    return 4. * c.pi * sim.star.R**2 * c.sigma_sb * sim.star.T**4
