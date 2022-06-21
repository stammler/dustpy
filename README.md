# DustPy

![GitHub](https://img.shields.io/github/license/stammler/dustpy) 
![PyPI - Downloads](https://img.shields.io/pypi/dm/dustpy) 
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](https://github.com/stammler/dustpy/blob/master/.github/CODE_OF_CONDUCT.md) 

### Dust Coagulation and Evolution in Protoplanetary Disks

`DustPy` is a Python package to simulate the evolution of dust in protoplanetary disks.

`DustPy` simulates the radial evolution of gas and dust in a protoplanetary disk, involving viscous evolution the the gas disk, advection and diffusion of the dust disk, as well as dust growth by solving the Smoluchowski equation.

Please read the [documentation](https://stammler.github.io/dustpy/) for a detailed description.


## Terms of Usage

Usage is only permitted under the condition of close collaboration and the option of co-authorship to publications that made use of `DustPy`, until the code paper is published.  
By using any version of `DustPy` you agree to these terms of usage.

## Installation

`DustPy` can be installed via PyPI:

`pip install dustpy`

## Requirements

`DustPy` needs a Python3 distribution and a Fortran compiler installed on your system.

## Contributing

To contribute to `DustPy` please read [Contribution Guidlines](https://github.com/stammler/dustpy/blob/master/.github/CONTRIBUTING.md).

## Documentation

[https://stammler.github.io/dustpy/](https://stammler.github.io/dustpy/)

* [1. Basics](https://stammler.github.io/dustpy/1_basics.html)
* [2. Simple Customization](https://stammler.github.io/dustpy/2_simple_customization.html)
* [3. Advanced Customization](https://stammler.github.io/dustpy/3_advanced_customization.html)
* [4. The Standard Model](https://stammler.github.io/dustpy/4_standard_model.html)
* [5. Dust Coagulation](https://stammler.github.io/dustpy/5_dust_coagulation.html)
* [6. Dust Evolution](https://stammler.github.io/dustpy/6_dust_evolution.html)
* [7. Gas Evolution](https://stammler.github.io/dustpy/7_gas_evolution.html)

* [Test: Analytical Coagulation Kernels](https://stammler.github.io/dustpy/test_analytical_coagulation_kernels.html)
* [Test: Gas Evolution](https://stammler.github.io/dustpy/test_gas_evolution.html)

* [Example: Ice Lines](https://stammler.github.io/dustpy/example_ice_lines.html)
* [Example: Planetary Gaps](https://stammler.github.io/dustpy/example_planetary_gaps.html)
* [Example: Planetesimal Formation](https://stammler.github.io/dustpy/example_planetesimal_formation.html)

* [List of Publications](https://stammler.github.io/dustpy/A_publications.html)

[Module Reference](https://stammler.github.io/dustpy/api.html)

`DustPy` is using the `simframe` framework.  
For details we refer to the [Simframe Documentation](https://simframe.rtfd.io/).

## Ackowledgements

`DustPy` has received funding from the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation programme under grant agreement No 714769.

`DustPy` was developed at the [University Observatory](https://www.usm.uni-muenchen.de/index_en.php) of the [Ludwig Maximilian University of Munich](https://www.en.uni-muenchen.de/index.html).