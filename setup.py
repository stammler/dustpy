import setuptools
import sys

try:
    from numpy.distutils.core import Extension
    from numpy.distutils.core import setup
except ImportError as exc:  # We do not have our build deps installed
    msg = "Error: {} must be installed before running the build.".format(
        exc.name)
    msg = "Please install NumPy first. You can do this with `pip install numpy`"
    print(msg)
    sys.exit(1)

# from setuptools import setup
import pathlib

package_name = "dustpy"
here = pathlib.Path(__file__).absolute().parent

# Fortran modules
ext_const = Extension(name="dustpy.constants._constants_f", sources=[
                      "dustpy/constants/constants.f90"])
ext_dust = Extension(name="dustpy.std.dust_f", sources=["dustpy/constants/constants.f90",
                                                        "dustpy/utils/interpolation.f90",
                                                        "dustpy/std/dust.f90"])
ext_gas = Extension(name="dustpy.std.gas_f", sources=["dustpy/constants/constants.f90",
                                                      "dustpy/utils/interpolation.f90",
                                                      "dustpy/std/gas.f90"])
extensions = [ext_const, ext_dust, ext_gas]


def read_version():
    with (here / package_name / '__init__.py').open() as fid:
        for line in fid:
            if line.startswith('__version__'):
                delim = '"' if '"' in line else "'"
                return line.split(delim)[1]
        else:
            raise RuntimeError("Unable to find version string.")


setup(
    name=package_name,

    description="Dust evolution in protoplanetary disks",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    keywords="numerical,simulation,science,physics,astrophysics,astronomy",

    url="https://github.com/stammler/dustpy",
    project_urls={"Source Code": "https://github.com/stammler/dustpy/",
                  "Documentation": "https://stammler.github.io/dustpy/"
                  },

    author="Sebastian Stammler, Til Birnstiel",
    author_email="sebastian.stammler@gmail.com, til.birnstiel@lmu.de",
    maintainer="Sebastian Stammler",

    version=read_version(),
    license="GPLv3",

    ext_modules=extensions,

    classifiers=["Development Status :: 5 - Production/Stable",
                 "Environment :: Console",
                 "Intended Audience :: Developers",
                 "Intended Audience :: Science/Research",
                 "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
                 "Natural Language :: English",
                 "Operating System :: OS Independent",
                 "Programming Language :: Python",
                 "Programming Language :: Python :: 3 :: Only",
                 "Programming Language :: Fortran",
                 "Topic :: Education",
                 "Topic :: Scientific/Engineering",
                 "Topic :: Scientific/Engineering :: Physics",
                 ],

    packages=setuptools.find_packages(),
    install_requires=["matplotlib", "numpy", "requests", "simframe>=0.5.3"],
    include_package_data=True,
    zip_safe=False,
)
