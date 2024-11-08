subroutine pt3(n, t1, x1, y1, alpha, corr, ci, taux, tauy)
  use result1
  use data2
  use pearsont3_module
  use random
  use time
  use setting, only: n1, n2
  implicit none
  integer, intent(in) :: n
  real(dp), intent(in) :: t1(n), x1(n), y1(n), alpha
  real(dp), intent(out) :: corr, ci(2), taux, tauy
  integer :: i1
  corr = 0.5
  ci = [0.4, 0.6]
  taux = 3
  tauy = 4
  call ranseed()
  n1 = size(t1)
  call init1a           ! n2=n1
  print *, 'after init1a'
  call calc_t_inv_lambda      ! calculates percentage point tv(lambda) over a
  print *, 'after calc_t_inv_lambda'
  call allocate1              ! t2, x2, y2, x3, y3, x3_resample1, y3_resample1,
  call allocate_resample_data
  call init1b                 ! t2, x2, y2, x3, y3, x3_resample1, y3_resample1,
  ! return
  ! call r_est       ! detrends (x2->x3, y2->y3), estimates r(x3, y3)
  ! return
  ! call tauest      ! estimates persistence times taux3, tauy3 and rhox3 and
  ! call chsett4     ! changes setting: l_mbb , block length
  ! call confidence  ! estimates [r_low; r_upp]
  ! r = corr
  ! ci = [r_low, r_upp]
  ! taux = taux3
  ! tauy = tauy3
end subroutine pt3
