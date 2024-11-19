! Subroutine interface to PearsonT3. See comments in pearsont3.f90
!

subroutine p3_subroutine(n, t, x, y, alpha, corr, ci, taux, tauy)
  use result1
  use pearsont3_module
  use data1
  use data2
  use random
  use time
  use setting, only: n1, n2
  implicit none
  integer, intent(in) :: n
  real(dp), intent(in) :: t(n), x(n), y(n), alpha
  real(dp), intent(out) :: corr, ci(2), taux, tauy
  integer :: i1
  logical :: allocate0_here
  call ranseed()
  n1 = n
  allocate0_here = .not.allocated(t1)
  if (allocate0_here) call allocate0  ! t1, x1, y1
  t1 = t
  x1 = x
  y1 = y
  n1 = n
  call init1a           ! n2=n1
  call calc_t_inv_lambda      ! calculates percentage point tv(lambda) over a
  !                             lambda grid (Calibrated CI)
  call allocate1              ! t2, x2, y2, x3, y3, x3_resample1, y3_resample1,
  !                             x3_resample2, Dy3_resample2
  call allocate_resample_data
  call init1b                 ! t2, x2, y2, x3, y3, x3_resample1,
  !                             y3_resample1,
  !                             x3_resample2, y3_resample2
  call r_est       ! detrends (x2->x3, y2->y3), estimates r(x3, y3)
  call tauest      ! estimates persistence times taux3, tauy3 and rhox3 and
  !                  rhoy3)
  call chsett4     ! changes setting: l_mbb , block length
  call confidence  ! estimates [r_low; r_upp]
  corr = r
  ci = [r_low, r_upp]
  taux = taux3
  tauy = tauy3
  print *,'at (1)'
  if (allocate0_here) then
    call deallocate0
    call deallocate_resample_data
    call deallocate1
  end if
end subroutine p3_subroutine
