python_sources = ['__init__.py']
py3.install_sources(python_sources, subdir: 'dustpy/constants')

constants_source = custom_target('_constants_fmodule.c',
  input: ['constants.f90'],
  output: ['_constants_fmodule.c', '_constants_f-f2pywrappers2.f90'],
  command: [py3, '-m', 'numpy.f2py', '--lower', '-m', '_constants_f', '@INPUT@', '--build-dir', 'dustpy/constants']
  )

py3.extension_module('_constants_f',
  ['constants.f90', constants_source],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies: py_dep,
  install: true,
  subdir: 'dustpy/constants')
