program fill_arrays
  use precision
  implicit none
  integer, parameter :: n = 10
  integer :: i, iter
  real(dp), allocatable :: x(:), y(:), t(:)
  real(dp), allocatable :: u1(:), u2(:)
  real(dp), parameter :: two_pi = 2.0_dp*acos(-1.0_dp)
  real(dp) :: confidence_level, corr, ci(2), taux, tauy
  interface
    subroutine p3_subroutine(n, t, x, y, confidence_level, corr, ci, taux, tauy)
      import :: dp
      integer, intent(in) :: n
      real(dp), intent(in) :: t(n), x(n), y(n)
      real(dp), intent(in) :: confidence_level
      real(dp), intent(out) :: corr, ci(2), taux, tauy
    end subroutine p3_subroutine
  end interface

  allocate(x(n), y(n), t(n), u1(n), u2(n))

  do iter = 1, 3
    call random_number(u1)
    call random_number(u2)

    where (u1 <= 0.0_dp)
      u1 = tiny(1.0_dp)
    end where

    x = sqrt(-2.0_dp*log(u1))*cos(two_pi*u2)
    y = sqrt(-2.0_dp*log(u1))*sin(two_pi*u2)
    t = [(real(i, dp), i = 1, n)]
    confidence_level = 0.05_dp

    print '(A,I0)', 'Iteration ', iter
    call p3_subroutine(n, t, x, y, confidence_level, corr, ci, taux, tauy)

    print '(A,1X,*(F12.6,1X))', 'x =', x
    print '(A,1X,*(F12.6,1X))', 'y =', y
    print '(A,1X,*(F12.6,1X))', 't =', t
    print '(A,F10.6)', 'corr = ', corr
    print '(A,1X,*(F10.6,1X))', 'ci =', ci
    print '(A,F10.6)', 'taux = ', taux
    print '(A,F10.6)', 'tauy =', tauy
  end do

  deallocate(x, y, t, u1, u2)
end program fill_arrays
