! Read a test case from stdin and exercise the legacy Fortran entry point.
program pearson_t3_reference
  use precision
  use setting, only: l_mbb, dtrtype
  implicit none
  integer :: n, b1, b2, n_lambda, i
  real(dp) :: confidence, r, ci(2), taux, tauy
  real(dp), allocatable :: t(:), x(:), y(:)
  read(*,*) n, b1, b2, n_lambda, confidence, dtrtype
  allocate(t(n), x(n), y(n))
  do i=1,n
    read(*,*) t(i), x(i), y(i)
  end do
  call p3_subroutine(n,t,x,y,1.0_dp-confidence,b1,b2,n_lambda,r,ci,taux,tauy)
  write(*,'(A,1X,5(ES25.17,1X),I0)') 'RESULT',r,ci,taux,tauy,l_mbb
end program
