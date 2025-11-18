! Subroutine interface to PearsonT3. See comments in pearsont3.f90
!

subroutine p3_subroutine(n, t, x, y, confidence_level, corr, ci, taux, tauy)
  use result1
  use pearsont3_module
  use data1
  use data2
  use random
  use time
  use setting, only: n1
  use parameters, only: alpha
  implicit none
  integer, intent(in) :: n
  real(dp), intent(in) :: t(n), x(n), y(n), confidence_level
  real(dp), intent(out) :: corr, ci(2), taux, tauy
  logical :: allocate0_here
  write(*,*) 'p3_subroutine: start, n =', n
  call flush(6)
  alpha = confidence_level
  call ranseed()
  n1 = n
  allocate0_here = .not.allocated(t1)
  write(*,*) 'p3_subroutine: allocate0_here =', allocate0_here
  call flush(6)
  if (allocate0_here) then
    write(*,*) 'p3_subroutine: calling allocate0'
    call flush(6)
    call allocate0  ! t1, x1, y1
  end if
  t1 = t
  x1 = x
  y1 = y
  n1 = n
  write(*,*) 'p3_subroutine: calling init1a'
  call flush(6)
  call init1a           ! n2=n1
  write(*,*) 'p3_subroutine: calling calc_t_inv_lambda'
  call flush(6)
  call calc_t_inv_lambda      ! calculates percentage point tv(lambda) over a
  !                             lambda grid (Calibrated CI)
  write(*,*) 'p3_subroutine: calling allocate1'
  call flush(6)
  call allocate1              ! t2, x2, y2, x3, y3, x3_resample1, y3_resample1,
  !                             x3_resample2, Dy3_resample2
  write(*,*) 'p3_subroutine: calling allocate_resample_data'
  call flush(6)
  call allocate_resample_data
  write(*,*) 'p3_subroutine: calling init1b'
  call flush(6)
  call init1b                 ! t2, x2, y2, x3, y3, x3_resample1,
  !                             y3_resample1,
  !                             x3_resample2, y3_resample2
  write(*,*) 'p3_subroutine: calling r_est'
  call flush(6)
  call r_est       ! detrends (x2->x3, y2->y3), estimates r(x3, y3)
  write(*,*) 'p3_subroutine: calling tauest'
  call flush(6)
  call tauest      ! estimates persistence times taux3, tauy3 and rhox3 and
  !                  rhoy3)
  write(*,*) 'p3_subroutine: calling chsett4'
  call flush(6)
  call chsett4     ! changes setting: l_mbb , block length
  write(*,*) 'p3_subroutine: calling confidence'
  call flush(6)
  call confidence  ! estimates [r_low; r_upp]
  corr = r
  ci = [r_low, r_upp]
  taux = taux3
  tauy = tauy3
  write(*,*) 'p3_subroutine: before deallocation, allocate0_here =', allocate0_here
  call flush(6)
  if (allocate0_here) then
    write(*,*) 'p3_subroutine: deallocating temporary arrays'
    call flush(6)
    call deallocate0
    call deallocate_resample_data
    call deallocate1
  end if
  write(*,*) 'p3_subroutine: end'
  call flush(6)
end subroutine p3_subroutine
