include 'nrtype.f90'
include 'nr.f90'
include 'nrutil.f90'
include 'erfcc.f90'
program testerfc
  use nr, only: erfcc
  implicit none
  double precision x
  real xs
  integer :: i, ix, ia, ifault
  print *, '    Built-in erfc     Numerical Recipes'
  do ix=1,5
    x = ix
    xs = ix
    write(*, '(i2  )', advance='no') ix
    write(*, '(1x,2f15.6)') erfc(x), erfcc(xs)
  end do
end program testerfc
