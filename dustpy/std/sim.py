'''Module containing standard functions for the main simulation object.'''

import dustpy.std.dust as std_dust
import dustpy.std.gas as std_gas

import numpy as np


def dt_adaptive(sim):
    """Function that returns the suggested adaptive timestep.
    By default DustPy uses adaptive integration schemes. The step
    size function is therefore simply returning the suggested
    step size."""
    return sim.t.suggested


def dt(sim):
    """Function returns the timestep depending on the source terms."""
    cfl = 0.9
    dt = 1.e100

    if np.any(sim.gas.S.tot[1:-1] < 0.):
        mask_gas = np.logical_and(
            sim.gas.Sigma > sim.gas.SigmaFloor,
            sim.gas.S.tot < 0.)
        mask_gas[0] = False
        mask_gas[-1] = False
        rate_gas = sim.gas.Sigma[mask_gas] / sim.gas.S.tot[mask_gas]
        try:
            dt_gas = np.min(np.abs(rate_gas))
            dt = np.minimum(dt, dt_gas)
        except:
            pass

    if np.any(sim.dust.S.tot[1:-1, ...] < 0.):
        mask_dust = np.logical_and(
            sim.dust.Sigma > sim.dust.SigmaFloor,
            sim.dust.S.tot < 0.)
        mask_dust[0, :] = False
        mask_dust[-1:, :] = False
        rate_dust = sim.dust.Sigma[mask_dust] / sim.dust.S.tot[mask_dust]
        try:
            dt_dust = np.min(np.abs(rate_dust))
            dt = np.minimum(dt, dt_dust)
        except:
            pass

    return cfl * dt


def prepare_explicit_dust(sim):
    """This function is the preparation function that is called
    before every integration step."""
    std_gas.prepare(sim)


def prepare_implicit_dust(sim):
    """This function is the preparation function that is called
    before every integration step."""
    std_gas.prepare(sim)
    std_dust.prepare(sim)


def finalize_explicit_dust(sim):
    """This function is the finalization function that is called
    after every integration step. It is managing the boundary
    conditions and is enforcing floor values."""
    std_gas.finalize(sim)
    std_dust.finalize_explicit(sim)


def finalize_implicit_dust(sim):
    """This function is the finalization function that is called
    after every integration step. It is managing the boundary
    conditions and is enforcing floor values."""
    std_gas.finalize(sim)
    std_dust.finalize_implicit(sim)
