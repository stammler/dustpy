'''``DustPy`` is a package for simulating dust coagulation and disk evolution in protoplanetary disks.

The class for performing simulations is ``dustpy.Simulation``.

Furthermore ``DustPy`` contains a simple package for plotting data ``dustpy.plot`` and a package
that contains pre-defined standard functions ``dustpy.std`` that can be used in the simulation.

``dustpy.utils`` contains some helper classes/functions that are used within the simulation.

``DustPy`` is mainly written in ``Python``. Computation intensive calculations are written in ``Fortran``.

``DustPy`` is using the ``simframe`` package for setting up scientific simulations.'''

from dustpy import plot
from dustpy.simulation import Simulation
from dustpy import constants
from dustpy import utils
from dustpy.utils import hdf5writer

from simframe.io.dump import readdump

__name__ = "dustpy"
__version__ = "1.0.2"

Simulation.__version__ = __version__
plot.__version__ = __version__
utils.__version__ = __version__

__all__ = ["constants", "hdf5writer", "plot", "readdump", "Simulation"]

# Print warning if dustpy version is outdated
utils.print_version_warning()
