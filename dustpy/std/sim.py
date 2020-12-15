'''Module containing standard functions for the main simulation object.'''

from dustpy import std

import numpy as np


def dt_adaptive(sim):
    """Function that returns the suggested adaptive timestep.
    By default DustPy uses adaptive integration schemes. The step
    size function is therefore simply returning the suggested
    step size.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame

    Returns
    dt : float
        Time step"""
    return sim.t.suggested


def dt(sim):
    """Function returns the timestep depending on the source terms.

    Paramters
    ---------
    sim : Frame
        Parent simulation frame

    Returns
    -------
    dt : float
        Time step"""

    dt_gas = std.gas.dt(sim) or 1.e100
    dt_dust = std.dust.dt(sim) or 1.e100
    dt = np.minimum(dt_gas, dt_dust)
    return sim.t.cfl * dt


def prepare_explicit_dust(sim):
    """This function is the preparation function that is called
    before every integration step.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame"""
    std.gas.prepare(sim)


def prepare_implicit_dust(sim):
    """This function is the preparation function that is called
    before every integration step.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame"""
    std.gas.prepare(sim)
    std.dust.prepare(sim)


def finalize_explicit_dust(sim):
    """This function is the finalization function that is called
    after every integration step. It is managing the boundary
    conditions and is enforcing floor values.

    Paramters
    ---------
    sim : Frame
        Parent simulation frame"""
    std.gas.finalize(sim)
    std.dust.finalize_explicit(sim)


def finalize_implicit_dust(sim):
    """This function is the finalization function that is called
    after every integration step. It is managing the boundary
    conditions and is enforcing floor values.

    Parameters
    ----------
    sim : Frame
        Parent simulation frame"""
    std.gas.finalize(sim)
    std.dust.finalize_implicit(sim)
