include 'as239.f90'
include 'nrtype.f90'
include 'nr.f90'
include 'nrutil.f90'
include 'gser.f90'
include 'gcf.f90'
include 'gammln.f90'
include 'gammp.f90'
program testgamma
  use nr, only: gammp
  implicit none
  double precision a, x, gammad, as239
  real numrecip, xs, as
  integer :: i, ix, ia, ifault
  character(len=*), parameter :: fmt = '(1x, f9.6)'
  print *, '    AS algorithm 239       Numerical Recipes'
  print *, '       1         5            1         5'
  do ix=1,10
    x = ix
    xs = ix
    write(*, '(i2  )', advance='no') ix
    do ia=1,5,4
      a = ia
      as239 = gammad(x, a, ifault)
      write(*, fmt, advance='no') as239
    end do
    write(*, '(3x)', advance='no')
    do ia=1,5,4
      as = ia
      numrecip = gammp(as, xs)
      write(*, fmt, advance='no') numrecip
    end do
    print *
  end do
end program testgamma
