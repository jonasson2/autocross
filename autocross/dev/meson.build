project('autocross-dev', 'fortran',
  default_options : ['warning_level=2', 'fortran_std=f2008'])

repo_root = join_paths(meson.project_source_root(), '..')

fill_arrays_sources = [
  join_paths(repo_root, 'src', 'common_modules.f90'),
  join_paths(repo_root, 'src', 'p3_modules.f90'),
  join_paths(repo_root, 'src', 'p3_subroutine.f90'),
  'fill_arrays.f90',
]

executable('fill_arrays', fill_arrays_sources, install : false)
