'''Module containing standard functions for the main simulation object.'''

import dustpy.std.dust as std_dust
import dustpy.std.gas as std_gas


def dt(sim):
    """Function calculates the timestep.
    By default DustPy uses adaptive integration schemes. The step
    size function is therefore simply returning the suggested
    step size."""
    return sim.t.suggested


def prepare(sim):
    """This function is the preparation function that is called
    before every integration step."""
    std_gas.prepare(sim)


def finalize(sim):
    """This function is the finalization function that is called
    after every integration step. It is managing the boundary
    conditions and is enforcing floor values."""
    std_gas.finalize(sim)
    std_dust.finalize(sim)
