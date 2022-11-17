# DustPy

[![GitHub Workflow Status](https://img.shields.io/github/workflow/status/stammler/dustpy/pages%20build%20and%20deployment?label=docs)](https://stammler.github.io/dustpy/) 
[![GitHub](https://img.shields.io/github/license/stammler/dustpy) ](https://github.com/stammler/dustpy/blob/master/LICENSE) 
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](https://github.com/stammler/dustpy/blob/master/.github/CODE_OF_CONDUCT.md)  
[![The Astrophysical Journal](https://img.shields.io/badge/The%20Astrophysical%20Journal-10.3847%2F1538--4357%2Fac7d58-blue)](https://doi.org/10.3847/1538-4357/ac7d58) [![arXiv](https://img.shields.io/badge/arXiv-10.48550/arXiv.2207.00322-blue)](https://doi.org/10.48550/arXiv.2207.00322)  
[![PyPI - Downloads](https://img.shields.io/pypi/dm/dustpy?label=PyPI%20downloads)](https://pypistats.org/packages/dustpy)

### Dust Coagulation and Evolution in Protoplanetary Disks

`DustPy` is a Python package to simulate the evolution of dust in protoplanetary disks.

`DustPy` simulates the radial evolution of gas and dust in a protoplanetary disk, involving viscous evolution the the gas disk, advection and diffusion of the dust disk, as well as dust growth by solving the Smoluchowski equation.

Please read the [documentation](https://stammler.github.io/dustpy/) for a detailed description.
By using any version of `DustPy` you agree to these terms of usage.

## Installation

`DustPy` can be installed via the Python Package Index

`pip install dustpy`

## Requirements

`DustPy` needs a Python3 distribution and a Fortran compiler installed on your system.

## Documentation

[https://stammler.github.io/dustpy/](https://stammler.github.io/dustpy/)

* [1. Basics](https://stammler.github.io/dustpy/1_basics.html)
* [2. Simple Customization](https://stammler.github.io/dustpy/2_simple_customization.html)
* [3. Advanced Customization](https://stammler.github.io/dustpy/3_advanced_customization.html)
* [4. The Standard Model](https://stammler.github.io/dustpy/4_standard_model.html)
* [5. Dust Coagulation](https://stammler.github.io/dustpy/5_dust_coagulation.html)
* [6. Dust Evolution](https://stammler.github.io/dustpy/6_dust_evolution.html)
* [7. Gas Evolution](https://stammler.github.io/dustpy/7_gas_evolution.html) <br /> &nbsp;
* [Test: Analytical Coagulation Kernels](https://stammler.github.io/dustpy/test_analytical_coagulation_kernels.html)
* [Test: Gas Evolution](https://stammler.github.io/dustpy/test_gas_evolution.html) <br /> &nbsp;
* [Example: Ice Lines](https://stammler.github.io/dustpy/example_ice_lines.html)
* [Example: Planetary Gaps](https://stammler.github.io/dustpy/example_planetary_gaps.html)
* [Example: Planetesimal Formation](https://stammler.github.io/dustpy/example_planetesimal_formation.html) <br /> &nbsp;
* [A. Citation](https://stammler.github.io/dustpy/A_citation.html)
* [B. List of Publications](https://stammler.github.io/dustpy/B_publications.html)
* [C. Contributing/Bugs/Features](https://stammler.github.io/dustpy/C_contrib_bug_feature.html)
* [D. DustPy Discussions](https://stammler.github.io/dustpy/D_discussions.html)
* [E. Changelog](https://stammler.github.io/dustpy/E_changelog.html)

[Module Reference](https://stammler.github.io/dustpy/api.html)

## Framework

`DustPy` is using the [Simframe](http://github.com/stammler/simframe/) framework for scientific simulations ([Stammler & Birnstiel 2022](https://joss.theoj.org/papers/0ef61e034c57445e846b2ec383c920a6))

## Acknowledgements

`DustPy` has received funding from the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation programme under grant agreement No 714769.

`DustPy` was developed at the [University Observatory](https://www.usm.uni-muenchen.de/index_en.php) of the [Ludwig Maximilian University of Munich](https://www.en.uni-muenchen.de/index.html).