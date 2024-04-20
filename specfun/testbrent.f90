module fun_module
  implicit none
  double precision :: y
contains
  function f(x,y) result(fxy)
    double precision, intent(in) :: x, y
    double precision fxy
    fxy = (x - y)**2 + 2d0
  end function f
  function f1(x) result(fx)
    double precision, intent(in) :: x
    double precision fx
    fx = f(x,y)
  end function f1
end module fun_module

program testbrent
  use fmin_module
  use fun_module
  double precision :: xmin, fxmin
  integer :: iy
  do iy = 1,4
    y = iy
    fxmin = fmin(f1, 0d0, 1d0, 3d0, 1d-8, xmin)
    print *, 'fxmin =', fxmin, ', xmin =', xmin
  end do
end program testbrent
