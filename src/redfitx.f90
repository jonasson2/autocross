! REDFIT-X
!=========
! Estimates cross-spectrum, coherency and phase spectrum from two unevenly
! spaced time series where the significance is evaluated with Monte Carlo
! simulations!
!
! The program is based on the previous program REDFIT (Schulz and Mudelsee,
! 2002), where autospectrum is estimated from unevenly spaced time series and
! tested against red noise (the AR(1) parameter is estimated directly from the
! unevenly spaced time series).
!
! Version 2.0, April 2024
!
! Main changes:
!       o Cross-spectral analysis has been implemented
!       o Significance measurements evaluated with Monte Carlo simulations
!       o Coherency Spectrum: Monte Carlo false alarm level
!       o Phase Spectrum: Monte Carlo confidence interval
!
! Changes in April 2024 by Kristján Jónasson:
!    1) Reformatted to have maximum line length 80
!    2) Removed some extra blank and blank comment lines
!    3) Changed tabs to spaces and reindented using emacs
!    4) Removed include statements
!    5) Changed everything to double precision
!    6) Moved modules from front to separate file modules.f90
!    7) Replaced Numerical Recipes (NR) functions:
!       NR       NOW USING           OBTAINED FROM
!       erfcc    erfc                Fortran built-in
!       sort     Latex dlasrt        netlib.org
!       gammp    AS algorithm 239    people.math.sc.edu/~jburkardt
!       brent    fmin                github.com/jacobwilliams/fmin  
!       ran      random_number       Fortran built-in
!       avevar   meanvar             own (new) subroutine
!    8) Moved subroutines and functions from end to modules in modules.f90
!    9) Updated KBO in and added KJ to author list
!    10) Added version number above
!    11) Now licensed with the MIT license
!        All functions used have a permissive license:
!        AS-239:  MIT
!        Lapack:  BSD-3
!        fmin:    BSD-3
!
! Authors: 
! ========
! Kristin B. Olafsdottir 
! Icelandic Meteorological Office
! E-mail: kbo@vedur.is
!
! Kristjan Jonasson
! University of Iceland
! jonasson@hi.is
!
! Michael Schulz, 
! MARUM and Faculty of Geosciences, Univ. Bremen
! E-mail:  mschulz@marum.de
!
! Manfred Mudelsee
! Climate Risk Analysis
! E-mail: mudelsee@climate-risk-analysis.com
!
! Reference: Olafsdottir, K.B., Schulz, M. and Mudelsee, M.
!    REDFIT-X: Cross-spectral analysis of unevenly spaced
!    paleoclimate time series. submitted.
!
! Reference: Schulz, M. and Mudelsee, M. (2002) REDFIT: Estimating
!            red-noise spectra directly from unevenly spaced paleoclimatic
!            time series. Computers and Geosciences, 28, 421-426.
!
! LICENCE (c) 2024 the authors; MIT license (https://opensource.org/license/mit)
!   [except AS239, DLASRT and FMIN; see above]


! Input:   parameter namelist; passed as configuration file
! ====
!          
!  &cfg  
!    fnin(1) =   x.dat', Input file name  for the 1st time series data
!    fnin(2) =   y.dat', Input file name  for the 2nd time series data

!    fnout =     ?????,  The results are written to files with this name (plain
!                        text files with various file extensions)

!    x_sign =    F,      [T/F] Change sign of first time series if T, default F
!    y_sign =    F,      [T/F] Change sign of second time series, default = F
!    nsim =      1000,   Number of simulations  (1000-2000 is recommended)
!    mctest =    T,      [T/F]  Estimate the significance of auto and
!                        coherency spectrum with Monte Carlo simulations
!    mctest_phi= T,      [T/F] Estimate Monte Carlo confidence interval for
!                        the phase spectrum
!    ofac =      4.0,    Oversampling value for Lomb-Scargle Fourier
!                        transform (typical values: 2.0-4.0)
!    hifac =     1.0,    Max. freq. = HIFAC * <f_Nyq> (Default = 1.0)
!    n50 =       8,      Number of WOSA segments (50 % overlap)      
!    alpha =     0.05,   Significance level
!    iwin =      1       Window type
!                        0: Rectangular
!                        1: Welch
!                        2: Hanning
!                        3: Triangular
!                        4: Blackman-Harris
!
!
! The two time series data need to be in two separated data files. 
! Format of the time-series file:
!
!            # comment lines
!            # .
!            # .
!            t(1)   x(1)
!            t(2)   x(2)
!             .       .
!             .       .
!            t(N)   x(N)
!
! where t(1) < t(2) < ... < t(N) are GEOLOGICAL AGES! The maximum number of
! data points N is only limited by the available amount of memory.

! Output:
! ======
! * Estimated parameters and spectra (including significane levels) are
!   written to FNOUT (self-explanatory!).
!
! * Error and warning messages are written to REDFIT-X.LOG.
!
! Notes:
! ------
! * A linear trend is subtracted from each WOSA segment.
!
! * tau is estimated separately for each WOSA segment and subsequently
!   averaged.
!
! * Default max. frequency = avg. Nyquist freq. (hifac = 1.0).
!
! * Input times must be provided as ages since a "reversed" geological time
!   vector is assumed in subroutine TAUEST.
!==========================================================================

program redfitx
  !
  use precision
  use const       ! some constants
  use timeser     ! defines tx, x, ty, y, npx, npy
  use param       ! defines namelist parameterxs (as redfit-x.cfg)
  use trigwindat  ! defines txcos, txsin, tycos, tysin, wxtau, wytau, wwx, wwy
  use nyquist  ! defines nsegx, nsegy, dfxy, avgtxy, fnyq, wz, nfreq, lfreq
  use phase    ! defines result vectors and more
  use random
  use meanvar_module
  use sort_module
  use redfit_x_module
  !
  implicit noene

  character(len = 80) :: cfgfile, fnout
  real(dp), allocatable, dimension(:,:)::  data_x, data_y, data_xy, data_cxy, &
    & data_phxy
  integer:: i
  
  cfgfile = 'redfit-x.cfg'
  open (10, file = cfgfile, form = 'formatted', status = 'old', iostat = iocheck)
  if (iocheck .ne. 0 ) stop 'cannot open config file'
  read(10, nml = cfg)
  close (10)
  namelist /cfg/ fnin, fnout, nsim, ofac, hifac, n50, iwin, alpha
  ! mctest, mctest_phi x_sign and y_sign have been removed

  call setdim(fnin)
  nx = npx
  ny = npy
  nout = nfreq
!
! setup workspace for input data
! ------------------------------
  allocate(x(npx), tx(npx))
  allocate(y(npy), ty(npy))
!
! retrieve time series data
! -------------------------
  call readdat(fnin)

  call rx_subroutine(nx, ny, nout, tx, ty, x, y, cfg_nsim, &
    & cfg_ofac, cfg_hifac, cfg_n50, cfg_alpha, cfg_i, rhox, rhoy, taux, tauy,&
    & df, 6dB, false_alarm, scale, data_x, data_y, data_xy, data_cxy,&
    & out_data_phxy
  !
  open (21, file = trim(fnout)//'.gxx', form = "formatted")
  open (22, file = trim(fnout)//'.gyy', form = "formatted")
  open (23, file = trim(fnout)//'.gxy', form = "formatted")
  open (24, file = trim(fnout)//'.cxy', form = "formatted")
  open (25, file = trim(fnout)//'.phxy', form = "formatted")

  write(21,'("#  ""Freq""",9x,"""Gxx""",7x,"""Gxx_corr""",4x,"""Gred_th""",6x,"""<Gred>""",5x,&
    & """CorrFac""",5x,"""90%-Chi2""",4x,"""95%-Chi2""",&
    & 4x,"""99%-Chi2""",5x,"""90%-MC""",6x,"""95%-MC""",6x,"""99%-MC""" )')
  do i = 1, nout
    write (21,'(1x,12(f12.6,2x))') data_x(i,:)
  end do

  write(22,'("#  ""Freq""",9x,"""Gyy""",7x,"""Gyy_corr""",4x,"""Gred_th""",6x,"""<Gred>""",5x,&
    & """CorrFac""",5x,"""90%-Chi2""",4x,"""95%-Chi2""",&
    & 4x,"""99%-Chi2""",5x,"""90%-MC""",6x,"""95%-MC""",6x,"""99%-MC""" )')
  do i = 1, nout
    write (22,'(1x,12(f12.6,2x))') data_y(i,:)
  end do

  write(23,'("#  ""Freq""",9x,"""Gxy""" )')    
  do i = 1, nout
    write (23,'(1x,2(f12.6,2x))') data_xy(i,:)
  end do

  write(24,'("#  ""Freq""",9x,"""Cxy""",9x,"""Csig""",7x,&
    &              """MC-Csig""",6x,"""90%-MC""",6x,"""95%-MC""",4x,"""99%-MC""" )')     
  do i = 1, nout
    write (24,'(1x,7(f12.6,2x))') data_cxy(i,:)
  end do

  write(25,'("#  ""Freq""",9x,"""Phxy""",6x,"""CI-low""",7x, """CI-up""",6x,"""CI-mc-low""",4x,&
    &    """CI-mc-up""")')  
  do i = 1, nout
    write (25,'(1x,6(f12.6,2x))') data_phxy(i,:)
  end do

  close(21)
  close(22)
  close(23)
  close(24)
  close(25)
  
end program redfitx
