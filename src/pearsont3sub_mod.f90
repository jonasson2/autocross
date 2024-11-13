! Subroutine interface to PearsonT3
!
! Estimating Pearson's correlation coefficient with calibrated bootstrap
! confidence interval from serially dependent time series
!
! LICENCE (c) 2024 the authors; MIT license (https://opensource.org/license/mit)
! [except FMIN which is BSD-3 licensed]
!
! Authors
! =======
! Manfred Mudelsee
! Climate Risk Analysis
! Schneiderberg 26
! 30167 Hannover
! Germany
! Email: mudelsee@mudelsee.com
! URL: http://www.mudelsee.com
!
! Kristín B. Ólafsdóttir
! Icelandic Meteorological Office
! E-mail: kbo@vedur.is
!
! Kristján Jónasson
! University of Iceland
! jonasson@hi.is
!=============================================================================
! Change log for version 1
! ========================
!
! Version  Date           Comments
!
! 1.00     July 2002      o original version
! 1.10     November 2007  o adapted to Gnuplot 4.2
!                         o mean detrending method added
! 1.20     March 2010     o subroutine bootstrap: 'p = 1.0 commented out
!                         o subroutine plot, point 3.1: simplified (but: has
!                           to be adapted in future)
! 1.30     June 2011      o compiler switch to gfortran
!
! KBO-change log (May 2013), main changes:
!  o Bootstrap method changed from stationary bootstrap to pairwise MBB and
!    block length selector changed
!  o Confidence intervals changed from BCa to calibrated bootstrap Student's t
!    CI
!=============================================================================
! Changes in April 2024 by Kristján Jónasson (version 2.0)
! 1) Reformatting (max line length 80, cleaned comments, reindented)
! 2) Changed everything to double precision
! 3) Replaced Numerical Recipes functions (erfcc, brent, avevar)
! 4) Removed interface blocks, now everything is in modules
! 5) Updated author list
! 6) Now licensed with MIT license
!=============================================================================
!module pearsont3sub_mod
!contains

subroutine pearsont3sub(n, t, x, y, alpha, corr, ci, taux, tauy)
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
  !
  ! 1.    Welcome
  !       =======
  call ranseed()
  !call ranseed(42)
  !
  ! 2.    Data
  !       =====
  !call chsett2    ! changes setting: n1
  n1 = n
  allocate0_here = .not.allocated(t1)
  if (allocate0_here) call allocate0  ! t1, x1, y1
  !call init0      ! t1, x1, y1
  !call read1      ! reads data
  t1 = t
  x1 = x
  y1 = y
  !
  ! 3.    Time interval extraction and calculation
  !       =====================================
  n1 = n
  call init1a           ! n2=n1
  call calc_t_inv_lambda      ! calculates percentage point tv(lambda) over a
  ! lambda grid (Calibrated CI)
  call allocate1              ! t2, x2, y2, x3, y3, x3_resample1, y3_resample1,
  ! x3_resample2, Dy3_resample2
  call allocate_resample_data
  call init1b                 ! t2, x2, y2, x3, y3, x3_resample1,
  ! y3_resample1,
  ! x3_resample2, y3_resample2
  call r_est       ! detrends (x2->x3, y2->y3), estimates r(x3, y3)
  call tauest      ! estimates persistence times taux3, tauy3 and rhox3 and
  ! rhoy3)
  call chsett4     ! changes setting: l_mbb , block length
  call confidence  ! estimates [r_low; r_upp]
  r = corr         !
  ci = [r_low, r_upp]
  taux = taux3
  tauy = tauy3
  print *,'at (1)'
  if (allocate0_here) then
    call deallocate0
    call deallocate_resample_data
    call deallocate1
  end if
end subroutine pearsont3sub

!end module pearsont3sub_mod
