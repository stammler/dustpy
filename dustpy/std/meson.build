python_sources = [
  '__init__.py',
  'dust.py',
  'gas.py',
  'grid.py',
  'sim.py',
  'star.py'
  ]

py3.install_sources(python_sources, subdir: 'dustpy/std')

dust_fortran_sources = [
  '../constants/constants.f90',
  '../utils/math.f90',
  '../utils/interpolation.f90',
  'dust.f90'
  ]

gas_fortran_sources = [
  '../constants/constants.f90',
  '../utils/interpolation.f90',
  'gas.f90'
  ]

dust_source = custom_target('dust_fmodule.c',
  input: dust_fortran_sources,
  output: ['dust_fmodule.c', 'dust_f-f2pywrappers2.f90'],
  command: [py3, '-m', 'numpy.f2py', '--lower', '-m', 'dust_f', '@INPUT@', '--build-dir', 'dustpy/std']
  )

py3.extension_module('dust_f',
  dust_fortran_sources+dust_source,
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies: py_dep,
  install: true,
  subdir: 'dustpy/std')

gas_source = custom_target('gas_fmodule.c',
  input: gas_fortran_sources,
  output: ['gas_fmodule.c', 'gas_f-f2pywrappers2.f90'],
  command: [py3, '-m', 'numpy.f2py', '--lower', '-m', 'gas_f', '@INPUT@', '--build-dir', 'dustpy/std']
  )

py3.extension_module('gas_f',
  gas_fortran_sources+gas_source,
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies: py_dep,
  install: true,
  subdir: 'dustpy/std')
