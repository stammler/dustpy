import dustpy.std.dust as std_dust
import dustpy.std.gas as std_gas


def dt(sim):
    '''Function calculates the timestep'''
    return sim.t.suggested


def finalize(sim):
    std_gas.finalize(sim)
    std_dust.finalize(sim)
