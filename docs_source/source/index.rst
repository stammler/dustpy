``DustPy`` Documentation
=========================================

| ``DustPy`` is a Python package to simulate the evolution of dust in protoplanetary disks.

| This documentation is for ``DustPy v1.0.6``.

| ``DustPy`` simulates the radial evolution of gas and dust in protoplanetary disks, including viscous evolution of the gas, advection and diffusion of the dust, as well as dust growth by solving the Smoluchowski equation.

| ``DustPy`` has been published in `Stammler & Birnstiel (2022) <https://iopscience.iop.org/article/10.3847/1538-4357/ac7d58>`_.

| ``DustPy`` can be installed via the Python Package Index
| ``pip install dustpy``

| ``DustPy`` requires a Python3 distribution and a Fortran compiler.

| ``DustPy`` is based on the ``Simframe`` framework for scientific simulation (`Stammler & Birnstiel 2022 <https://joss.theoj.org/papers/10.21105/joss.03882>`_).
| Please have a look at the `Simframe Documentation <https://simframe.rtfd.io/>`_ for details of its usage.

| ``dustpylib`` is a collection of auxiliary tools and extensions for ``DustPy``, containing for example interfaces to radiative transfer codes.
| For more details, please have a look at the `dustpylinbdocumentation <https://dustpylib.rtfd.io/>`_.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   1_basics
   2_simple_customization
   3_advanced_customization
   4_standard_model
   5_dust_coagulation
   6_dust_evolution
   7_gas_evolution
   dustpylib
   test_analytical_coagulation_kernels
   test_gas_evolution
   example_ice_lines
   example_planetary_gaps
   example_planetesimal_formation
   A_citation
   B_publications
   C_contrib_bug_feature
   D_discussions
   E_changelog
   api



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
