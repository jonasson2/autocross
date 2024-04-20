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
! Version 2, April 2024
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
!    rhopre(1) = -999.0, Prescribed value for rho for the first time series,
!                        not used if rho < 0 (default = -999.0)
!    rhopre(2) = -999.0, Prescribed value for rho for the second time series,
!                        not used if rho < 0 (default = -999.0)
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

program redfit
  !
  use precision
  use const
  use timeser
  use param
  use trigwindat
  use error
  use nyquist
  use phase
  use random
  use meanvar_module
  use sort_module
  use redfit_x_module
  !
  implicit none
  !
  real(dp), dimension(:), allocatable :: &
    freq, &      ! frequency vector
    gxx, &       ! autospectrum of input data x
    gyy, &       ! autospectrum of input data y
    gxy, &       ! cross-spectrum 
    cxy, &       ! coherency spectrum
    phxy,&       ! phase spectrum                 
    gxxc, &      ! corrected autospectrum of input data x
    gyyc, &      ! corrected autospectrum of input data y
    grxxsum, &   ! sum of AR(1) spectra x
    gryysum, &   ! sum of AR(1) spectra y
    grxysum, &   ! sum of cross-spectrum bivariate AR(1)
    grxxavg, &   ! average AR(1) spectrum x
    gryyavg, &   ! average AR(1) spectrum y   
    grxyavg, &   ! average cross-spectrum bivariate AR(1)
    gredthx, &   ! theoretical AR(1) spectrum x
    gredthy, &   ! theoretical AR(1) spectrum y
    corrx, &     ! correction factor x
    corry        ! correction factor y
  real(dp), dimension(:,:), allocatable :: ci90, & ! 90% false-alarm level from MC
    ci95, &      ! 95% false-alarm level from MC
    ci99         ! 99% false-alarm level from MC
  real(dp), dimension(:,:),allocatable :: ephi !confidence interval for phase angle
  !                                        ! AR1 - spectrum
  real(dp), dimension(:,:), allocatable :: grxx  ! AR(1) spectra x
  real(dp), dimension(:,:), allocatable :: gryy  ! AR(1) spectra y
  real(dp), dimension(:,:), allocatable :: grxy  ! Cross-spectrum bivariate AR(1)
  real(dp), dimension(:,:), allocatable :: crxy  ! Coherency spectrum bivar. AR(1)
  real(dp), dimension(:,:), allocatable :: phrxy ! Phase spectrum bivariate AR(1)
  real(dp), dimension(:), allocatable :: redx ! AR(1) time series - based on time
  !                                       ! series x
  real(dp), dimension(:), allocatable :: redy    ! AR(1) time series - based on time 
  !                                          ! series y  
  real(dp) :: taux,tauy, rnsim, fac90, fac95, fac99, dof, neff, &
    avgdtx, avgdty, facx, facy, facxy, rhox, rhoxsq, rhoy, rhoysq, varx, vary, &
    varrx, varry, alphacritx, alphacrity, faccritx, faccrity,varxy, &
    varrxy, cobias, facphi, z, csig, se_dummy
  integer:: kstart, kstop, krate, kmax, ntime
  integer :: i, iocheck, ialloc
  integer :: idx90, idx95, idx99, iostat
  logical :: ini, biascorr
  character (len = 80) :: cfgfile,fnout
  character (len = 80), dimension(n_fnin) :: fnin
  character (len = 80) :: errorfile
  !
  namelist /cfg/ fnin, fnout, nsim, ofac, hifac, n50, iwin, mctest, alpha, &
    rhopre,mctest_phi, x_sign, y_sign
  !
  call system_clock(kstart, krate, kmax)
  !
  ! First try to direct error messages to stderr
  ! --------------------------------------------
  errorfile = '/dev/stderr'
  open(errio, file=errorfile, status='old', iostat=iostat, action='write')
  if (iostat /= 0) then 
    ! open log file instead
    ! ---------------------
    errorfile = 'redfit-x.log'
    open(errio, file = "redfit-x.log")
  endif
  print *,'Error messages written to ', errorfile
  !
  ! welcome message
  !----------------------------------------  
  print '(a)',' '
  print '(a)',' '
  print '(a)', '==============================================================='
  print '(a)', '                                                               '
  print '(a,a,a,a)', vers,'(Version 2.0',',',' April 2024)                  '
  print '(a)', '                                                               '
  print '(a)', '==============================================================='
  print '(a)', '                                                               '
  !
  ! retrieve command line arguments
  ! -------------------------------

  do i= 1,ntry
    print*
    print '(a)', &
      'Please type in the name of the configuration file [path + filename]'
    !read (5,'(a)') cfgfile
    cfgfile = 'redfit-x.cfg'
    open (10, file = cfgfile, form = 'formatted', status = 'old', &
      iostat = iocheck)
    if (iocheck .ne. 0 ) then
      print *,'Error - Cant''t open config file - try again!'
    else 
      exit
    end if
  end do
  if (i > ntry) then
    print *,'OK - REDFIT-X terminates.'
    stop
  end if
  read(10, nml = cfg)
  close (10)

  !
  ! workspace dimensions
  ! --------------------
  !
  print *,'xxxx'
  call setdim(fnin,maxdimx,maxdimy,nout)
  print *,'yyyy'
  !
  ! setup workspace for input data
  ! ------------------------------
  allocate(x(maxdimx), tx(maxdimx), stat = ialloc)
  if (ialloc .ne. 0) call allocerr("a")
  allocate(y(maxdimy), ty(maxdimy), stat = ialloc)
  if (ialloc .ne. 0) call allocerr("a")
  !
  ! retrieve time series data
  ! -------------------------
  call readdat(fnin)
  !
  ! check input time axis      
  ! ----------------------
  call check1()
  ! 
  ! in case of duplicate entries reinitialize array dimensions
  ! and retrieve averaged time series
  ! ----------------------------------------------------------
  if (errflagx .eqv. .true.) then
    write(*,'(/1x,a/)') &
      "Resetting dimensions after correcting for duplicate sampling times...x"
    deallocate(x, tx, stat = ialloc)
    if (ialloc .ne. 0) call allocerr("d")
    fnin(1) = "TimeSeriesx.avg"
    call setdim(fnin,maxdimx,maxdimy,nout)
    allocate(x(maxdimx), tx(maxdimx), stat = ialloc)
    if (ialloc .ne. 0) call allocerr("a")
    call readdat(fnin)
  end if
  !
  ! check input time axis
  ! ----------------------- 
  ! 
  call check2()
  !
  ! in case of duplicate entries reinitialize array dimensions
  ! and retrieve averaged time series
  ! ----------------------------------------------------------
  if (errflagy .eqv. .true.) then
    write(*,'(/1x,a/)') &
      "Resetting dimensions after correcting for duplicate sampling times...y"
    deallocate(y, ty, stat = ialloc)
    if (ialloc .ne. 0) call allocerr("d")
    fnin(1) = "TimeSeriesy.avg"
    call setdim(fnin,maxdimx,maxdimy,nout)
    allocate(y(maxdimy), ty(maxdimy), stat = ialloc)
    if (ialloc .ne. 0) call allocerr("a")
    call readdat(fnin)
  end if
  !
  !
  !change sign of the data if that is requested
  !--------------------------------------------------
  if (x_sign.eqv..true.)then
    x(:) = -x(:)
  end if

  if (y_sign.eqv..true.)then
    y(:) = -y(:)
  end if

  ! allocate remaining workspace
  ! ----------------------------
  allocate(redx(maxdimx), redy(maxdimy), stat = ialloc)
  if (ialloc .ne. 0) call allocerr("a")
  !
  allocate (freq(nout), gxx(nout), gyy(nout), gxy(nout), &
    cxy(nout), phxy(nout), stat = ialloc)
  if (ialloc .ne. 0) call allocerr("a")
  !
  allocate (grxxsum(nout),gryysum(nout),grxysum(nout), &
    grxxavg(nout),gryyavg(nout),grxyavg(nout), stat = ialloc)
  if (ialloc .ne. 0) call allocerr("a")
  !
  allocate (gredthx(nout), gredthy(nout), corrx(nout), corry(nout), &
    gxxc(nout), gyyc(nout), stat = ialloc)
  if (ialloc .ne. 0) call allocerr("a")
  !  
  if (mctest .eqv. .true.) then
    allocate (ci90(nspect,nout), ci95(nspect,nout), &  
      ci99(nspect,nout), stat = ialloc)
    if (ialloc .ne. 0) call allocerr("a")
  end if
  !  
  allocate (ephi(2,nout), stat = ialloc)
  if (ialloc .ne. 0) call allocerr("a")
  !
  ! average dt of entire time series
  ! --------------------------------
  avgdtx = sum(tx(2:npx)-tx(1:npx-1)) / real(npx-1,dp)
  avgdty = sum(ty(2:npy)-ty(1:npy-1)) / real(npy-1,dp)
  !
  ! determine autospectrum,crossspectrum, coherency and phase spectrum of
  ! input data !step 2
  ! ------------------------------------
  ini = .true.
  !
  call spectr(ini, tx(1:npx), x(1:npx), ty(1:npy), y(1:npy), ofac, n50, &
    iwin, freq, gxx, gyy, gxy, cxy, phxy)
  !
  ! estimate data variance from autospectrum
  ! ---------------------------------------------------
  varx = freq(2) * sum(gxx(1:nout))    ! NB: freq(2) = dfxy
  vary = freq(2) * sum(gyy(1:nout))    ! NB: freq(2) = dfxy  
  varxy = freq(2) * sum(gxy(1:nout))  
  !
  ! estimate tau unless tau is prescribed; die gracefully in case of an error
  ! -------------------------------------------------------------------------
  call gettau(rhopre(1), tx(1:npx), x(1:npx), npx, taux)
  if (ierr .eq. 1) then
    write (errio,*) ' Error in GETTAU'
    close(errio)
    write(*,*) &
      "An error has occurred. Check REDFIT-X.LOG for further information."
    stop
  end if
  call gettau(rhopre(2), ty(1:npy), y(1:npy), npy, tauy)
  if (ierr .eq. 1) then
    write (errio,*) ' Error in GETTAU'
    close(errio)
    write(*,*) &
      "An error has occurred. Check REDFIT-X.LOG for further information."
    stop
  end if

  print*
  print'(A5,F10.2)','taux',taux
  print'(A5,F10.2)','tauy',tauy
  print*
  !
  ! Generate NSim AR(1) Spectra
  ! ---------------------------------
  !
  if (mctest .eqv. .true.) then
    allocate(grxx(nsim,nout),gryy(nsim,nout),grxy(nsim,nout), &
      crxy(nsim,nout),phrxy(nsim,nout), stat = ialloc)
    if (ialloc .ne. 0) call allocerr("a")
  else
    allocate(grxx(1,nout),gryy(1,nout),grxy(1,nout),          & 
      crxy(1,nout),phrxy(1,nout), stat = ialloc)
    if (ialloc .ne. 0) call allocerr("a")
  end if
  !  
  call ranseed
  !  
  grxxsum(:) = 0
  gryysum(:) = 0
  grxysum(:) = 0
  grxxavg(:) = 0
  gryyavg(:) = 0
  grxyavg(:) = 0
  !
  call getdof(iwin, n50, dof, neff)
  !
  rnsim = real(nsim,dp)
  do i = 1, nsim  
    if ((mod(i,50) .eq. 0) .or. (i .eq. 1)) write(*,*) 'ISim =', i
    !
    !   setup AR(1) time series and estimate its spectrum
    !   -------------------------------------------------
    !
    call makear1(tx, npx, taux, redx)
    call makear1(ty, npy, tauy, redy)   
    !    
    ini = .false.
    if (mctest .eqv. .true.) then
      call spectr(ini, tx(1:npx), redx(1:npx), ty(1:npy), redy(1:npy),    &
        ofac, n50, iwin, freq, grxx(i,:), gryy(i,:),         &
        grxy(i,:), crxy(i,:), phrxy(i,:))
    else  
      call spectr(ini, tx(1:npx), redx(1:npx), ty(1:npy), redy(1:npy),    &
        ofac, n50, iwin, freq, grxx(1,:), gryy(1,:),         &
        grxy(1,:), crxy(1,:), phrxy(1,:))
    end if
    !
    !    scale and sum red-noise spectra
    !    -------------------------------
    !
    if (mctest .eqv. .true.) then
      varrx = freq(2) * sum(grxx(i,1:nout))   ! NB: freq(2) = df
      facx = varx / varrx
      grxx(i,1:nout) = facx * grxx(i,1:nout)  !scale grxx such that the area
      !under grxx is identical to the
      !area under gxx
      grxxsum(1:nout) = grxxsum(1:nout) + grxx(i,1:nout)
      !        
      varry = freq(2) * sum(gryy(i,1:nout))    ! NB: freq(2) = df
      facy = vary / varry
      gryy(i,1:nout) = facy * gryy(i,1:nout)  !scale gryy such that the area
      !under gryy is identical to the
      !area under gyy
      gryysum(1:nout) = gryysum(1:nout) + gryy(i,1:nout)
      !
      varrxy = freq(2) * sum(grxy(i,1:nout))    ! NB: freq(2) = df
      facxy = varxy / varrxy
      grxy(i,1:nout) = facxy * grxy(i,1:nout)  !scale grxy such that the area
      !under grxy is identical to the
      !area under gxy
      grxysum(1:nout) = grxysum(1:nout) + grxy(i,1:nout)
    else
      varrx = freq(2) * sum(grxx(1,1:nout))    ! NB: freq(2) = df
      facx = varx / varrx
      grxxsum(1:nout) = grxxsum(1:nout) + facx * grxx(1,1:nout)
      !         
      varry = freq(2) * sum(gryy(1,1:nout))    ! NB: freq(2) = df
      facy = vary / varry
      gryysum(1:nout) = gryysum(1:nout) + facy * gryy(1,1:nout)    
      !            
      varrxy = freq(2) * sum(grxy(i,1:nout))    ! NB: freq(2) = df
      facxy = varxy / varrxy
      grxysum(1:nout) = grxysum(1:nout) + facxy *grxy(1,1:nout)
    end if
  end do  !end of nsim loop
  !
  ! determine average red-noise spectrum; scale average again to
  ! make sure that roundoff errors do not affect the scaling
  ! ------------------------------------------------------------
  grxxavg(1:nout) = grxxsum(1:nout) / rnsim
  varrx = freq(2) * sum(grxxavg(1:nout))
  facx = varx / varrx
  grxxavg(1:nout) = facx * grxxavg(1:nout)  !Average AR(1) spectrum x
  !  
  gryyavg(1:nout) = gryysum(1:nout) / rnsim
  varry = freq(2) * sum(gryyavg(1:nout))
  facy = vary / varry
  gryyavg(1:nout) = facy * gryyavg(1:nout)  !Average AR(1) spectrum y
  !
  grxyavg(1:nout) = grxysum(1:nout) / rnsim  !Average Cross-spectrum bivariate
  !AR(1)
  varrxy = freq(2) * sum(grxyavg(1:nout))
  facxy = varxy / varrxy
  grxyavg(1:nout) = facxy * grxyavg(1:nout)  
  !
  ! determine lag-1 autocorrelation coefficient
  ! -------------------------------------------
  rhox = exp (-avgdtx / taux)              ! avg. autocorrelation coefficient x
  rhoxsq = rhox * rhox
  ! 
  rhoy = exp (-avgdty / tauy)              ! avg. autocorrelation coefficient y
  rhoysq = rhoy * rhoy
  !
  ! set theoretical spectrum (e.g., Mann and Lees, 1996, Eq. 4)
  ! make area equal to that of the input time series
  ! -----------------------------------------------------------
  fnyq = freq(nout)                   ! Nyquist freq.
  !
  !theoretical spectrum based on rho estimated from the time series x 
  gredthx(1:nout) = &
    (1-rhoxsq) / (1-2*rhox*cos(pi*freq(1:nout)/fnyq)+rhoxsq)
  varrx = freq(2) * sum(gredthx(1:nout))
  facx = varx / varrx
  gredthx(1:nout) = facx * gredthx(1:nout)
  !
  !theoretical spectrum based on rho estimated from the time series y  
  gredthy(1:nout) = &
    (1-rhoysq) / (1-2*rhoy*cos(pi*freq(1:nout)/fnyq)+rhoysq)
  varry = freq(2) * sum(gredthy(1:nout))
  facy = vary / varry        !step 5: Select G0 (for eq 2)
  gredthy(1:nout) = facy * gredthy(1:nout)      
  !
  ! determine correction factor
  ! --------------------------------
  !correction factor for the bias adjustment of the Lomb-Scargle spectrum
  corrx(1:nout) = grxxavg(1:nout) / gredthx(1:nout)  
  corry(1:nout) = gryyavg(1:nout) / gredthy(1:nout)

  ! correct for bias in autospectrum
  ! -------------------------------------
  gxxc(1:nout) = gxx(1:nout) / corrx(1:nout)
  gyyc(1:nout) = gyy(1:nout) / corry(1:nout)
  !
  !  red-noise false-alarm levels from percentiles of MC simulation
  ! --------------------------------------------------------------
  if (mctest .eqv. .true.) then
    do i = 1, nout
      call sort(grxx(1:nsim, i))
    end do
    !     
    do i = 1, nout
      call sort(gryy(1:nsim, i))
    end do
    !
    do i = 1, nout
      call sort(grxy(1:nsim, i))
    end do
    !
    do i = 1, nout
      call sort(crxy(1:nsim, i))
    end do
    !
    !    set percentil indices
    !    ---------------------
    !
    idx90 = int(0.90_dp * rnsim)
    idx95 = int(0.95_dp * rnsim)
    idx99 = int(0.99_dp * rnsim)
    !
    !   find frequency-dependent percentil and apply bias correction for
    !   autospectrum
    !   ----------------------------------------------------------------
    !   for gxx
    do i = 1, nout
      ci90(1,i) = grxx(idx90, i) / corrx(i)
      ci95(1,i) = grxx(idx95, i) / corrx(i)
      ci99(1,i) = grxx(idx99, i) / corrx(i)
    end do
    !
    !!for gyy
    do i = 1, nout
      ci90(2,i) = gryy(idx90, i) / corry(i)
      ci95(2,i) = gryy(idx95, i) / corry(i)
      ci99(2,i) = gryy(idx99, i) / corry(i)
    end do
    !
    !for gxy   (not used)
    do i = 1, nout
      ci90(3,i) = grxy(idx90, i)
      ci95(3,i) = grxy(idx95, i) 
      ci99(3,i) = grxy(idx99, i) 
    end do
    !     
    !for cxy 
    do i = 1, nout
      ci90(4,i) = crxy(idx90, i)
      ci95(4,i) = crxy(idx95, i) 
      ci99(4,i) = crxy(idx99, i) 
    end do
  end if
  !
  ! scaling factors for red noise from chi^2 distribution
  ! -----------------------------------------------------
  call getdof(iwin, n50, dof, neff)
  print '(a)'
  print'(A4,4x,F10.2)','dof', dof
  print'(A5,3x,F10.2)','neff', neff
  print'(A6,2x,F10.2)', 'alpha', alpha
  !  
  fac90 = getchi2(dof, 0.10_dp) / dof
  fac95 = getchi2(dof, 0.05_dp) / dof
  fac99 = getchi2(dof, 0.01_dp) / dof
  if (ierr .eq. 1) stop
  !
  !coherency - False alarm level (theoretical)    
  !-------------------------------------------
  !
  if(n50 == 1) then
    csig = 0
  else 
    csig = 1 - alpha**(1/(neff-1))
  end if
  !  
  ! coherency - mean Monte Carlo false alarm level
  !-----------------------------------------------
  if (mctest .eqv. .true.) then
    if (alpha==0.1)then
      csig_mc = sum(ci90(4,:))/nout
    else if (alpha==0.05_dp)then
      csig_mc = sum(ci95(4,:))/nout
    else if (alpha==0.01_dp)then
      csig_mc = sum(ci99(4,:))/nout
    else 
      csig_mc = 0
    end if
  else 
    csig_mc =-999
  end if

  print '(a)'
  print'(A43,F5.2,A1,F10.2)', &
    ' Mean MC false alarm for coherency (alpha =',alpha,')', csig_mc
  print'(A38,11x,F10.2)', 'Theoretical false alarm for coherency', csig
  !
  ! bias correction coherency spectrum  
  ! ----------------------------------
  !
  biascorr = .true.
  if(biascorr .eqv. .true.) then        
    do i = 1,nout
      if(cxy(i) < 1) then
        cobias = ((1-cxy(i)) *(1-cxy(i)) )/ neff
        cxy(i) = cxy(i) - cobias
      end if
      if (cxy(i)< 0) then       
        cxy(i) = 0
      end if
    end do
  end if
  !
  !Phase Spectrum - Confidence intervals (theoretical)
  !---------------------------------------------------
  !
  z = getz(alpha/2)  
  facphi = z *1/sqrt(2*neff) * 180/pi    
  !  
  do i= 1,nout
    !avoid overflow error if cxy = 1 or cxy = 0
    if(cxy(i) .lt. csig_mc) then
      ephi(1,i) = -999
      ephi(2,i) = -999
    else if(cxy(i).ge. csig_mc) then   
      if(cxy(i)>0 .and. cxy(i)<1) then 
        ephi(1,i) = sqrt(1-cxy(i))/ sqrt(cxy(i)) * facphi
        ephi(2,i) = ephi(1,i)
        ephi(1,i) = phxy(i)+ephi(1,i)
        ephi(2,i) = phxy(i)-ephi(2,i)
      else if (cxy(i)==1 )then
        ephi(1,i) = 0
        ephi(2,i) = 0
      else if (cxy(i) == 0) then 
        ephi(1,i) = phxy(i) + 180
        ephi(2,i) = 180 - phxy(i)
      end if
    end if
  end do
  !
  !Phase Spectrum - Monte Carlo confidence interval 
  !---------------------------------------------------------------
  ! 
  if((mctest_phi .eqv. .true.) .and. (mctest.eqv. .true.))then
    allocate (gbxx(nout), gbyy(nout), gbxy(nout), stat = ialloc)
    if (ialloc .ne. 0) call allocerr("a")
    allocate (cbxy(nsim, nout),phbxy(nsim, nout), stat = ialloc)
    if (ialloc .ne. 0) call allocerr("a")
    allocate (ephi_b(2, nout, 2), stat = ialloc)    ! 1:lowver boundary,
    ! 2:upper boundary ,
    ! 1:time-scale x
    ! 2:time-scale y
    if (ialloc .ne. 0) call allocerr("a")
    allocate (se_phbxy(2,nout,2), stat = ialloc)    !  1:lowver boundary,
    !  2:upper boundary , 1:
    !  time-scale x, 2:
    !  time-scale y
    if (ialloc .ne. 0) call allocerr("a")
    allocate (ephi_mc(2, nout), stat = ialloc)    ! Mean over the two time
    ! scales 1:lowver boundary,
    ! 2:upper boundary
    if (ialloc .ne. 0) call allocerr("a")
    allocate (se_mc_phxy(2, nout), stat = ialloc)    ! Mean over the two time
    ! scales 1:lowver
    ! boundary, 2:upper
    ! boundary
    if (ialloc .ne. 0) call allocerr("a")  
  end if
  !
  idxph_low = int((alpha/2)*nsim)
  idxph_up = int ((1-alpha/2)*nsim)
  ! 
  if((mctest_phi .eqv. .true.) .and. (mctest.eqv. .true.))then
    print '(a)'
    print*, 'Monte Carlo confidence interval for phase . . .'
    do sg = 1,2    !two loops. 1: time steps from time series x used. 2: time
      !steps from time series y used
      !
      If(sg==1)then
        np_xy = npx
      else if (sg==2)then
        np_xy = npy
      end if
      !
      allocate(t_xy(np_xy),stat = ialloc)
      if (ialloc .ne. 0) call allocerr("a")

      if(sg==1)then
        t_xy(:) = tx(:)
        tau_xy = taux
      else if(sg ==2)then
        t_xy(:) = ty(:)
        tau_xy = tauy
      end if
      !
      nsegx = int(2 * np_xy / (n50 + 1))             
      nsegy= int(2 * np_xy / (n50 + 1))             
      !
      allocate (redxb(np_xy), redyb(np_xy), stat = ialloc)
      if (ialloc .ne. 0) call allocerr("a")
      !
      do i = 1,nout            
        if(cxy(i) .ge.csig_mc) then                   !Confidence interval is
          !only formed at
          !frequencies where the
          !coherency is
          !significant
          g = sqrt(((-2 *sqrt(1-cxy(i)))/cxy(i))+(2/cxy(i))-(1))
          ! g value - used to generate two time series with known coherency.
          ! cxy(i) has been bias corrected.
          !
          deallocate (txsin, txcos,tysin, tycos, wxtau, wytau, wwx, wwy,  &
            stat = ialloc)
          if (ialloc .ne. 0) call allocerr("d")
          !    
          do i_boot = 1,nsim
            call make_coherar1(redxb,redyb)  ! 2 red noise time series with
            ! known coherency are generated
            if(i_boot==1)then
              ini = .true.
            else
              ini = .false.
            end if
            call spectr(ini, t_xy(1:np_xy), redxb(1:np_xy), t_xy(1:np_xy), &
              redyb(1:np_xy), ofac, n50, iwin, freq, gbxx, gbyy, &
              gbxy, cbxy(i_boot,:), phbxy(i_boot,:))
          end do
          print '(a6,i4,a12,i3)', 'freq' ,i, 'time-sc', sg
          !    
          call meanvar(phbxy(1:nsim,i), ave_phbxy, var_phbxy)  
          se_dummy = 0
          se_dummy = sqrt(var_phbxy) 
          z = getz(alpha/2)
          !       
          se_phbxy(1,i,sg) = phxy(i) + z * se_dummy    ! Standard error based
          ! CI. Lower boundary
          ! standard error
          ! *z(alpha).
          se_phbxy(2,i,sg) = phxy(i) -  z * se_dummy    ! Standard error
          ! bassed CI. Upper
          ! boundary

          call sort(phbxy(1:nsim,i))        ! for percentile based CI
          !      
          phbxy(:,i) = phbxy(:,i) - ave_phbxy       !The mean phase value from
          !the Monte Carlo ensemble
          !subtracted from the phase
          !values to center the
          !values around zero
          !    
          ephi_b(1,i,sg) = phxy(i) + phbxy(idxph_low,i)  !The estimated phase
          !value from the
          !observed series
          !added to the upper
          !and lower boundary
          !to form the
          !confidence interval
          ephi_b(2,i,sg) = phxy(i) + phbxy(idxph_up, i)

        else if (cxy(i).lt. csig_mc) then    
          ephi_b(1,i,sg) = -999
          ephi_b(2,i,sg) = -999
          se_phbxy(1,i,sg) = -999
          se_phbxy(2,i,sg) = -999
        end if
      end do
      !  
      deallocate (t_xy, redyb, redxb,  stat = ialloc)
      if (ialloc .ne. 0) call allocerr("d")
      !   
    end do
    !
    ! Form the final Monte Carlo confidence interval by taking the mean value
    ! for both time scalces
    !  
    do i = 1,nout
      ephi_mc(1,i)=  (ephi_b(1,i,1)+ephi_b(1,i,2))/2
      ephi_mc(2,i)=  (ephi_b(2,i,1) + ephi_b(2,i,2))/2
      se_mc_phxy(1,i) = (se_phbxy(1,i,1) + se_phbxy(1,i,2))/2 !not printed in
      !result file
      se_mc_phxy(2,i) = (se_phbxy(2,i,1) + se_phbxy(2,i,2))/2 !not printed in
      !result file
    end do
  End if
  !  
  nsegx = int(2 * npx / (n50 + 1)) ! nsegx and nsegy set to previous values
  nsegy= int(2 * npy / (n50 + 1))                 
  !
  ! critical false alarm level after Thomson (1990)
  ! -----------------------------------------------
  ! For autospectrum xx
  alphacritx = 1 / real(nsegx,dp)
  faccritx = getchi2(dof, alphacritx) / dof
  if (ierr .eq. 1) stop
  ! 
  ! For autospectrum yy
  alphacrity = 1 / real(nsegy,dp)
  faccrity = getchi2(dof, alphacrity) / dof
  if (ierr .eq. 1) stop
  !
  ! save results of AR(1) fit
  ! -------------------------
  call system_clock(kstop, krate, kmax)
  if (kstop .ge. kstart) then
    ntime = (kstop-kstart) / krate
  else ! kmax overflow
    ntime = ((kmax-kstart)+kstop) / krate
  end if

  ! Write result file - Autospectrum X
  ! ----------------------------------------
  open (20, file = trim(fnout)//'.gxx', form = "formatted", iostat = iocheck)
  if (iocheck .ne. 0 ) then
    write (errio, *) ' Error - Can''t create ', trim(fnout)
    close(errio)
    write(*,*) "An error has occurred. Check REDFIT-X.LOG for further information."
    stop
  end if
  write (20,*) '# ', vers, '    Autospectrum X'
  write (20,*) '#'
  write (20,*) '# Input:'
  write (20,*) '# ------'
  write (20,*) '# File = ', trim(fnin(1)) 
  write (20,*) '# OFAC = ', ofac
  write (20,*) '# HIFAC = ', hifac
  write (20,*) '# n50 = ', n50
  write (20,*) '# Iwin = ', iwin
  write (20,*) '# Nsim = ', Nsim
  write (20,*) '#'
  write (20,*) '# Initial values:'
  write (20,*) '# ---------------'
  write (20,*) '# Data variance (from data spectrum) = ', varx
  write (20,*) '# Avg. dtx = ', avgdtx
  write (20,*) '#'
  write (20,*) '# Results:'
  write (20,*) '# --------'
  if (rhopre(1) .lt. 0) then
    write (20,*) '# Avg. autocorr. coeff., rho = ', rhox
  else
    write (20,*) '# PRESCRIBED avg. autocorr. coeff., rho = ', rhopre(1)
  end if
  write (20,*) '# Avg. tau = ', taux
  write (20,*) '# Degrees of freedom = ', dof
  write (20,*) '# 6-dB Bandwidth = ', winbw(iwin, freq(2)-freq(1), ofac)
  write (20,*) '# Critical false-alarm level (Thomson, 1990) = ', &
    (1-alphacritx) * 100
  write (20,*) '#    ==> corresponding scaling factor for red noise = ', faccritx
  write (20,*) '#'
  write (20,*) '#'
  write (20,*) '# Elapsed time [s] = ', ntime
  write (20,*) '#'
  write (20,*) '# Data Columns:'
  write (20,*) '# -------------'
  write (20,*) '#  1: Freq = frequency'
  write (20,*) '#  2: Gxx = spectrum of input data'
  write (20,*) '#  3: Gxx_corr = bias-corrected spectrum of input data'
  write (20,*) '#  4: Gred_th = theoretical AR(1) spectrum'
  write (20,*) &
    '#  5: <Gred> = average spectrum of Nsim AR(1) time series (uncorrected)'
  write (20,*) '#  6: CorrFac = Gxx / Gxx_corr'
  write (20,*) '#  7: 90%-Chi2 = 90-% false-alarm level (Chi^2)'
  write (20,*) '#  8: 95%-Chi2 = 95-% false-alarm level (Chi^2)'
  write (20,*) '#  9: 99%-Chi2 = 99-% false-alarm level (Chi^2)'
  if (mctest .eqv. .true.) then
    write (20,*) '# 10: 90%-MC = 90-% false-alarm level (MC)'
    write (20,*) '# 11: 95%-MC = 95-% false-alarm level (MC)'
    write (20,*) '# 12: 99%-MC = 99-% false-alarm level (MC)'
    write (20,*) '#'
    write(20,'("#  ""Freq""",9x,"""Gxx""",7x,"""Gxx_corr""",4x,"""Gred_th""",&
      &               6x, """<Gred>""",5x,&
      &               """CorrFac""",5x,"""90%-Chi2""",4x,"""95%-Chi2""",&
      &               4x,"""99%-Chi2""",5x,"""90%-MC""",6x,"""95%-MC""",6x,&
      &               """99%-MC""" )')
    do i = 1, nout
      write (20,'(1x,13(e12.6,2x))') freq(i), gxx(i),gxxc(i), gredthx(i), &
        grxxavg(i), corrx(i), gredthx(i)*fac90, &
        gredthx(i)*fac95, gredthx(i)*fac99, ci90(1,i), &
        ci95(1,i), ci99(1,i)
    end do
  else
    write (20,*) '#'
    write(20,'("#  ""Freq""",9x,"""Gxx""",7x,"""Gxx_corr""",4x,"""Gred_th""",&
      &               6x,"""<Gred>""",5x,&
      &               """CorrFac""",4x,"""90%-Chi2""",4x,"""95%-Chi2""",&
      &               4x,"""99%-Chi2""")')    
    do i = 1, nout
      write (20,'(1x,9(e12.6,2x))') freq(i), gxx(i),gxxc(i), gredthx(i), &
        grxxavg(i), corrx(i), gredthx(i)*fac90, &
        gredthx(i)*fac95, gredthx(i)*fac99
    end do
  end if
  close(20)
  !  
  ! Write result file - Autospectrum Y
  ! ----------------------------------------
  open (20, file = trim(fnout)//'.gyy', form = "formatted", iostat = iocheck)
  if (iocheck .ne. 0 ) then
    write (errio, *) ' Error - Can''t create ', trim(fnout)
    close(errio)
    write(*,*) &
      "An error has occurred. Check REDFIT-X.LOG for further information."
    stop
  end if
  write (20,*) '# ', vers, '    Autospectrum Y'
  write (20,*) '#'
  write (20,*) '# Input:'
  write (20,*) '# ------'
  write (20,*) '# File = ', trim(fnin(2)) 
  write (20,*) '# OFAC = ', ofac
  write (20,*) '# HIFAC = ', hifac
  write (20,*) '# n50 = ', n50
  write (20,*) '# Iwin = ', iwin
  write (20,*) '# Nsim = ', Nsim
  write (20,*) '#'
  write (20,*) '# Initial values:'
  write (20,*) '# ---------------'
  write (20,*) '# Data variance (from data spectrum) = ', vary
  write (20,*) '# Avg. dtx = ', avgdty
  write (20,*) '#'
  write (20,*) '# Results:'
  write (20,*) '# --------'
  if (rhopre(2) .lt. 0) then
    write (20,*) '# Avg. autocorr. coeff., rho = ', rhoy
  else
    write (20,*) '# PRESCRIBED avg. autocorr. coeff., rho = ', rhopre(2)
  end if
  write (20,*) '# Avg. tau = ', tauy
  write (20,*) '# Degrees of freedom = ', dof
  write (20,*) '# 6-dB Bandwidth = ', winbw(iwin, freq(2)-freq(1), ofac)
  write (20,*) '# Critical false-alarm level (Thomson, 1990) = ', &
    (1-alphacrity) * 100
  write (20,*) '#    ==> corresponding scaling factor for red noise = ', &
    faccrity
  write (20,*) '#'
  write (20,*) '#'
  write (20,*) '# Elapsed time [s] = ', ntime
  write (20,*) '#'
  write (20,*) '# Data Columns:'
  write (20,*) '# -------------'
  write (20,*) '#  1: Freq = frequency'
  write (20,*) '#  2: Gyy = spectrum of input data'
  write (20,*) '#  3: Gyy_corr = bias-corrected spectrum of input data'
  write (20,*) '#  4: Gred_th = theoretical AR(1) spectrum'
  write (20,*) &
    '#  5: <Gred> = average spectrum of Nsim AR(1) time series (uncorrected)'
  write (20,*) '#  6: CorrFac = Gyy / Gyy_corr'
  write (20,*) '#  7: 90%-Chi2 = 90-% false-alarm level (Chi^2)'
  write (20,*) '#  8: 95%-Chi2 = 95-% false-alarm level (Chi^2)'
  write (20,*) '#  9: 99%-Chi2 = 99-% false-alarm level (Chi^2)'
  if (mctest .eqv. .true.) then
    write (20,*) '# 10: 90%-MC = 90-% false-alarm level (MC)'
    write (20,*) '# 11: 95%-MC = 95-% false-alarm level (MC)'
    write (20,*) '# 12: 99%-MC = 99-% false-alarm level (MC)'
    write (20,*) '#'     
    write(20,'("#  ""Freq""",9x,"""Gyy""",7x,"""Gyy_corr""",4x,"""Gred_th""",&
      &               6x,"""<Gred>""",5x,&
      &               """CorrFac""",5x,"""90%-Chi2""",4x,"""95%-Chi2""",&
      &               4x,"""99%-Chi2""",5x,"""90%-MC""",6x,"""95%-MC""",6x,&
      &               """99%-MC""" )')    
    do i = 1, nout
      write (20,'(1x,12(e12.6,2x))') freq(i), gyy(i),gyyc(i), gredthy(i), &
        gryyavg(i), corry(i), gredthy(i)*fac90, &
        gredthy(i)*fac95, gredthy(i)*fac99, ci90(2,i), &
        ci95(2,i), ci99(2,i)
    end do
  else
    write (20,*) '#'
    write(20,'("#  ""Freq""",9x,"""Gyy""",7x,"""Gyy_corr""",4x,"""Gred_th""",&
      &               6x,"""<Gred>""",5x,&
      &               """CorrFac""",4x,"""90%-Chi2""",4x,"""95%-Chi2""",&
      &               4x,"""99%-Chi2""")')    
    do i = 1, nout
      write (20,'(1x,9(e12.6,2x))') freq(i), gyy(i),gyyc(i), gredthy(i), &
        gryyavg(i), corry(i), gredthy(i)*fac90, &
        gredthy(i)*fac95, gredthy(i)*fac99
    end do
  end if
  close(20)
  !  
  ! Write result file - Cross-spectrum XY
  ! -------------------------------------------
  open (20, file = trim(fnout)//'.gxy', form = "formatted", iostat = iocheck)
  if (iocheck .ne. 0 ) then
    write (errio, *) ' Error - Can''t create ', trim(fnout)
    close(errio)
    write(*,*) &
      "An error has occurred. Check REDFIT-X.LOG for further information."
    stop
  end if
  write (20,*) '# ', vers, '    Cross-spectrum XY'
  write (20,*) '#'
  write (20,*) '# Input:'
  write (20,*) '# ------'
  write (20,*) '# File = ', trim(fnin(1))
  write (20,*) '# File = ', trim(fnin(2)) 
  write (20,*) '# OFAC = ', ofac
  write (20,*) '# HIFAC = ', hifac
  write (20,*) '# n50 = ', n50
  write (20,*) '# Iwin = ', iwin
  write (20,*) '# Nsim = ', Nsim
  write (20,'(a22,f6.3)') '# Level of signif. = ', alpha
  write (20,*) '#'
  write (20,*) '# Initial values:'
  write (20,*) '# ---------------'
  write (20,*) '# Avg. dtxy = ', avgdty
  write (20,*) '#'
  write (20,*) '# Results:'
  write (20,*) '# --------'
  if (rhopre(1) .lt. 0) then
    write (20,*) '# Avg. autocorr. coeff., rhox = ', rhox
  else
    write (20,*) '# PRESCRIBED avg. autocorr. coeff., rho = ', rhopre(1)
  end if
  if (rhopre(2) .lt. 0) then
    write (20,*) '# Avg. autocorr. coeff., rhoy = ', rhoy
  else
    write (20,*) '# PRESCRIBED avg. autocorr. coeff., rho = ', rhopre(2)
  end if
  write (20,*) '# Avg. taux = ', taux
  write (20,*) '# Avg. tauy = ', tauy
  write (20,*) '# Degrees of freedom = ', dof
  write (20,*) '# 6-dB Bandwidth = ', winbw(iwin, freq(2)-freq(1), ofac)
  write (20,*) '#'
  write (20,*) '# Elapsed time [s] = ', ntime
  write (20,*) '#'
  write (20,*) '# Data Columns:'
  write (20,*) '# -------------'
  write (20,*) '#  1: Freq = frequency'
  write (20,*) '#  2: Gxy = cross-spectrum of input data'
  write (20,*) '#'
  write(20,'("#  ""Freq""",9x,"""Gxy""" )')    
  do i = 1, nout
    write (20,'(1x,3(e12.6,2x))') freq(i), gxy(i)
  end do
  close(20)
  !
  ! Write result file - Coherency spectrum XY
  ! -------------------------------------------------
  open (20, file = trim(fnout)//'.cxy', form = "formatted", iostat = iocheck)
  if (iocheck .ne. 0 ) then
    write (errio, *) ' Error - Can''t create ', trim(fnout)
    close(errio)
    write(*,*) &
      "An error has occurred. Check REDFIT-X.LOG for further information."
    stop
  end if
  write (20,*) '# ', vers, '    Coherency-spectrum XY'
  write (20,*) '#'
  write (20,*) '# Input:'
  write (20,*) '# ------'
  write (20,*) '# File = ', trim(fnin(1))
  write (20,*) '# File = ', trim(fnin(2)) 
  write (20,*) '# OFAC = ', ofac
  write (20,*) '# HIFAC = ', hifac
  write (20,*) '# n50 = ', n50
  write (20,*) '# Iwin = ', iwin
  write (20,*) '# Nsim = ', Nsim
  write (20,'(a22,f6.3)') '# Level of signif. = ', alpha
  write (20,*) '#'
  write (20,*) '# Initial values:'
  write (20,*) '# ---------------'
  write (20,*) '# Avg. dtxy = ', avgdty
  write (20,*) '#'
  write (20,*) '# Results:'
  write (20,*) '# --------'
  if (rhopre(1) .lt. 0) then
    write (20,*) '# Avg. autocorr. coeff., rhox = ', rhox
  else
    write (20,*) '# PRESCRIBED avg. autocorr. coeff., rho = ', rhopre(1)
  end if
  if (rhopre(2) .lt. 0) then
    write (20,*) '# Avg. autocorr. coeff., rhoy = ', rhoy
  else
    write (20,*) '# PRESCRIBED avg. autocorr. coeff., rho = ', rhopre(2)
  end if
  write (20,*) '# Avg. taux = ', taux
  write (20,*) '# Avg. tauy = ', tauy
  write (20,*) '# Degrees of freedom = ', dof
  write (20,*) '# 6-dB Bandwidth = ', winbw(iwin, freq(2)-freq(1), ofac)
  write (20,*) '#'
  write (20,*) '# Elapsed time [s] = ', ntime
  write (20,*) '#'
  write (20,*) '# Data Columns:'
  write (20,*) '# -------------'
  write (20,*) '#  1: Freq = Frequency'
  write (20,*) '#  2: Cxy = Coherency-spectrum of input data'
  write (20,*) '#  3: Csig = Theoretical False alarm level'
  if (mctest .eqv. .true.) then
    write (20,fmt='(a64,f6.3,a)') &
      '#  4: MC-CSig = Mean Monte Carlo false-alarm level (for alpha =', &
      alpha,')'
    write (20,*) '#  5: 90%-MC = 90% Monte Carlo false-alarm level'
    write (20,*) '#  6: 95%-MC = 95% Monte Carlo false-alarm level'
    write (20,*) '#  7: 99%-MC = 99% Monte Carlo false-alarm level'
    write (20,*) '#'     
    write(20,'("#  ""Freq""",9x,"""Cxy""",9x,"""Csig""",7x,&
      &              """MC-Csig""",6x,"""90%-MC""",6x,"""95%-MC""",4x,&
      &              """99%-MC""" )')     
    do i = 1, nout
      write (20,'(1x,7(f12.6,2x))') freq(i), cxy(i), csig, &
        csig_mc, ci90(4,i), ci95(4,i), ci99(4,i)  
    end do
  else
    write (20,*) '#'
    write(20,'("#  ""Freq""",9x,"""Cxy""", 9x,"""Csig""" )')    
    do i = 1, nout
      write (20,'(1x,3(f12.6,2x))') freq(i), cxy(i), csig
    end do
  end if
  close(20)
  !  
  !
  ! Write result file - Phase spectrum XY
  ! --------------------------------------------
  open (20, file = trim(fnout)//'.phxy', form = "formatted", iostat = iocheck)
  if (iocheck .ne. 0 ) then
    write (errio, *) ' Error - Can''t create ', trim(fnout)
    close(errio)
    write(*,*) &
      "An error has occurred. Check REDFIT-X.LOG for further information."
    stop
  end if
  write (20,*) '# ', vers, '    Phase spectrum XY'
  write (20,*) '#'
  write (20,*) '# Input:'
  write (20,*) '# ------'
  write (20,*) '# File = ', trim(fnin(1))
  write (20,*) '# File = ', trim(fnin(2)) 
  write (20,*) '# OFAC = ', ofac
  write (20,*) '# HIFAC = ', hifac
  write (20,*) '# n50 = ', n50
  write (20,*) '# Iwin = ', iwin
  write (20,*) '# Nsim = ', nsim
  write (20,'(a22,f6.3)') '# Level of signif. = ', alpha
  write (20,*) '#'
  write (20,*) '# Initial values:'
  write (20,*) '# ---------------'
  write (20,*) '# Avg. dtxy = ', avgdty
  write (20,*) '#'
  write (20,*) '# Results:'
  write (20,*) '# --------'
  if (rhopre(1) .lt. 0) then
    write (20,*) '# Avg. autocorr. coeff., rhox = ', rhox
  else
    write (20,*) '# PRESCRIBED avg. autocorr. coeff., rho = ', rhopre(1)
  end if
  if (rhopre(2) .lt. 0) then
    write (20,*) '# Avg. autocorr. coeff., rhoy = ', rhoy
  else
    write (20,*) '# PRESCRIBED avg. autocorr. coeff., rho = ', rhopre(2)
  end if
  write (20,*) '# Avg. taux = ', taux
  write (20,*) '# Avg. tauy = ', tauy
  write (20,*) '# Degrees of freedom = ', dof
  write (20,*) '# 6-dB Bandwidth = ', winbw(iwin, freq(2)-freq(1), ofac)
  write (20,*) '#'
  write (20,*) '# Elapsed time [s] = ', ntime
  write (20,*) '#'
  write (20,*) '# Data Columns:'
  write (20,*) '# -------------'
  write (20,*) '#  1: Freq = Frequency'
  write (20,*) '#  2: Phxy = Phase-spectrum of input data'
  write (20,*) '#  3: CI-low = Theoretical Confidence Interval - lower'
  write (20,*) '#  4: CI-up =  Theoretical Confidence Interval - upper'  
  if ((mctest_phi.eqv..true.).and.(mctest.eqv..true.))then
    write (20,*) &
      '#  5: CI-mc-low = Monte Carlo Confidence Interval - lower  (percentiles)'
    write (20,*) &
      '#  6: CI-mc-up =  Monte Carlo Confidence Interval - upper  (percentiles)'
    write (20,*) '#' 
    write(20,'("#  ""Freq""",9x,"""Phxy""",6x,"""CI-low""",7x, """CI-up""",6x,&
      &    """CI-mc-low""",4x,&
      &    """CI-mc-up""")')  
    do i = 1, nout
      write (20,'(1x,6(f12.6,2x))') freq(i), phxy(i), ephi(1,i),ephi(2,i), &
        ephi_mc(1,i), ephi_mc(2,i)
    end do
  Else
    write (20,*) '#' 
    write(20,'("#  ""Freq""",9x,"""Phxy""",6x,"""CI-low""",6x, """CI-up""")')  
    do i = 1, nout
      write (20,'(1x,4(f12.6,2x))') freq(i), phxy(i), ephi(1,i),ephi(2,i)
    end do
  end if
  close(20)        
  !
  ! clean up
  ! --------
  if (ierr .eq. 0) then
    inquire(unit=errio, name=errorfile)
    if (errorfile /= '/dev/stderr') close(errio, status = "delete")
  else
    close(errio, status = "keep")
    write(*,*) &
      "An error has occurred. Check REDFIT-X.LOG for further information."
  end if
  deallocate(x, tx, redx, y, ty, redy, stat = ialloc)
  if (ialloc .ne. 0) call allocerr("d")
  deallocate (freq, gxx, gyy, gxy, cxy, phxy, stat = ialloc)
  if (ialloc .ne. 0) call allocerr("d")
  deallocate(grxx, gryy, grxy, crxy, phrxy, stat = ialloc)
  if (ialloc .ne. 0) call allocerr("d")
  deallocate(grxxsum, gryysum, grxysum,      &
    grxxavg, gryyavg, grxyavg, stat = ialloc)
  if (ialloc .ne. 0) call allocerr("d")
  deallocate (gredthx, gredthy, corrx , corry,    &
    gxxc, gyyc, stat = ialloc)
  if (ialloc .ne. 0) call allocerr("d")
  if (mctest .eqv. .true.) then
    deallocate (ci90, ci95, ci99, stat = ialloc)
    if (ialloc .ne. 0) call allocerr("d")
  end if
  if((mctest_phi.eqv..true.).and. (mctest.eqv..true.))then
    deallocate (gbxx,gbyy,gbxy,cbxy,phbxy,ephi_b,se_phbxy, ephi_mc, &
      se_mc_phxy,  stat = ialloc)
    if (ialloc .ne. 0) call allocerr("d")
  end if
end program redfit
