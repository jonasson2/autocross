include 'meanvar.f90'
include 'avevar.f90'

program testmeanvar
  use nr, only: avevar
  use meanvar_module
  double precision :: x(5) = [1, 2, 3, 4, 5], mean, var
  real :: xs(5) = [1, 2, 3, 4, 5], ave, varnr
  call avevar(xs, ave, varnr)
  call meanvar(x, mean, var)
  print *, 'From Numerical Recipes: ', ave, varnr
  print *, 'From meanvar.f90:       ', mean, var
end program testmeanvar
