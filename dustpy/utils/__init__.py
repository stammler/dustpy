'''Package containing utility classes and functions used in the simulation.'''

from dustpy.utils.boundary import Boundary
from dustpy.utils.data import read_data
from dustpy.utils.simplenamespace import SimpleNamespace
from dustpy.utils.version import print_version_warning

__all__ = [
    "Boundary",
    "read_data",
    "print_version_warning",
    "SimpleNamespace",
]
__version__ = None
