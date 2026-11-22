project('autocross', 'c', 'fortran',
        version: '0.1.0.dev0', default_options: ['optimization=3'])

py = import('python').find_installation(pure: false)
np_inc = run_command(py, '-c', 'import numpy; print(numpy.get_include())',
                     check: true).stdout().strip()
# Generate C with the build environment's Python, not the system python3.pc.
core_c = custom_target('autocross_cython', input: 'autocross/_core.pyx',
  output: '_core.c', command: [py, '-m', 'cython', '-3', '-o', '@OUTPUT@', '@INPUT@'])
py.extension_module('_core',
  core_c,
  files('native/randompack.f90',
        '../autocross/src/common_modules.f90',
        '../autocross/src/p3_modules.f90',
        '../autocross/src/p3_python.f90'),
  dependencies: py.dependency(),
  c_args: ['-I' + np_inc, '-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION'],
  link_language: 'fortran',
  install: true,
  subdir: 'autocross',
)
py.install_sources('autocross/__init__.py', 'autocross/_testcase.py', 'autocross/_t3.py',
                   subdir: 'autocross')
