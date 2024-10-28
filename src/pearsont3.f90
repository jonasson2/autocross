 ! PearsonT, version 2.0, April 2024
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


program pear
  use pearsont3sub_module
  use result1
  use pearsont3_module
	use data1
  use random
  use time
  implicit none
  character(len=1) :: c1
  integer :: i1
	real(dp) :: ci(2)
  !
  ! 1.    Welcome
  !       =======
  call tic()
  call info0          ! welcome message
  call ranseed()
  !call ranseed(42)
  !
  ! 2.    Data
  !       =====
  call chsett1    ! changes setting: data file name
  call chsett2    ! changes setting: n1
  call allocate0  ! t1, x1, y1
  call init0      ! t1, x1, y1
  call read1      ! reads data
  call pearsont3sub(t1, x1, y1, 0.05d0, r, ci, taux3, tauy3) 
  r_low = ci(1)
  r_upp = ci(2)
  
  !
  ! 3.    Time interval extraction and calculation
  !       =====================================
	call info1('p')  ! extracted time interval
	call info2('p')  ! persistence times
	call info3('p')  ! detrending method 
	call info4('p')  ! r [ r_low; r_upp ]
	call tocprint('after info4')
	call decis1(c1)  ! decision tree: Part 1
  !
  ! 4.    Output and exit
  !       ==============
  call output
  call deallocate1
end program pear
