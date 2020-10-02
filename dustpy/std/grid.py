'''Module containing standard functions for the grids.'''

import dustpy.constants as c

import numpy as np


def OmegaK(sim):
    """Function calculates the Keplerian frequency.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    OmegaK : Field
        Keplerian frequency"""
    return np.sqrt(c.G * sim.star.M / sim.grid.r**3)
