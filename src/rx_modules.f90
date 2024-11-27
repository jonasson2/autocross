! Here are modules taken from redfitx.f90 in version 1. These are the
! modules that are only used by Redfitx.

module const
  use precision
  implicit none
  public
  character (len = 11) :: vers = 'REDFIT-X'
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: tpi = 2 * pi
  !integer :: nout                 !nout = nfreq
  integer,parameter :: n_fnin = 2  !number of input files, 1: timeseries x, 2: timeseries y 
  integer,parameter :: nspect = 4  !1.Autospectrum x
  !2.Autospectrum y
  !3.Cross-spectrum xy
  !4.Coherency spectrum xy              
end module const
!
! time series data
! ----------------
module timeser
  use precision
  implicit none
  public
  real(dp), dimension(:), allocatable :: tx, x , ty, y
  integer :: npx, npy      !number of data points
end module timeser
!
! parameter
! ---------
module param
  use precision
  implicit none
  public
  integer, parameter :: ntry=3     ! user input: maximum number of errors
  real(dp) :: ofac
  real(dp) :: hifac = 1
  integer :: nsim, n50, iwin
  logical :: mctest = .false.
  logical::mctest_phi=.false.
  real(dp) :: alpha 
  logical:: x_sign = .false.
  logical:: y_sign = .false.
end module param
!
! trigonometric data for FT
! -------------------------
module trigwindat
  use precision
  implicit none
  public
  real(dp), dimension(:,:,:), allocatable :: txcos, txsin, tycos, tysin
  real(dp), dimension(:,:), allocatable :: wxtau, wytau, wwx, wwy
end module trigwindat
!
! nyquist frequency variables
! ---------------------------
module nyquist
  use precision
  implicit none
  public
  integer :: nsegx, nsegy ! data points per segment x and y
  real(dp) :: dfxy            ! max frequency spacing, max(dfx,dfy)
  ! real(dp) :: avgdtxy       ! max average sampling interval, max(avgdtx,avgdty)
  real(dp) :: fnyq            ! average Nyquist frequency
  real(dp) :: wz              ! omega = 2*pi*f
  integer :: nfreq        ! f(1) = f(0) ; f(nfreq) = fNyq  (number of frequencies
  integer :: lfreq        ! nfreq * 2
end module nyquist
!
!Monte Carlo phase confidence interval variables
!--------------------------------
module phase
  use precision
  implicit none
  public
  real(dp) :: csig_mc
  real(dp) :: g
  integer :: i_boot
  real(dp),dimension (:), allocatable :: redxb, redyb
  real(dp), dimension(:), allocatable :: gbxx, gbyy, gbxy
  real(dp), dimension(:,:), allocatable :: cbxy, phbxy
  integer :: idxph_low, idxph_up
  real(dp), dimension(:,:,:), allocatable :: ephi_b
  real(dp), dimension(:,:), allocatable :: ephi_mc
  real(dp) :: ave_phbxy, var_phbxy
  real(dp), dimension(:,:,:), allocatable :: se_phbxy
  real(dp), dimension(:,:), allocatable :: se_mc_phxy
  integer:: np_xy
  integer:: sg
  real(dp), dimension(:), allocatable :: t_xy
  real(dp) :: tau_xy
end module phase

module redfit_x_module
  use precision
  implicit none

contains
  !--------------------------------------------------------------------------
  subroutine setdim(fnin)  
    !--------------------------------------------------------------------------
    ! Analyze data file and set array dimensions.
    !--------------------------------------------------------------------------

    use timeser
    use param
    use nyquist
    use const, only : n_fnin, pi
    !
    implicit none
    !
    character (len = 80), dimension (n_fnin), intent(in) :: fnin 
    !
    real(dp) :: tdum, xdum, t1  
    real(dp) :: avgdtx, avgdty    !average sampling interval for x and y
    real(dp) :: avgdtxy
    real(dp) :: tpx, tpy      !average period of segment x and y
    real(dp) :: dfx, dfy      !frequency spacing for x and y
    integer :: i, iocheck, j
    character (len = 1) :: flag
    ! j loop - opens both inputfiles x and y and sets number of data points npx and npy
    ! --------------------------------------------------------------------------------
    do j = 1,n_fnin
      !
      !   open input file
      !   ---------------
      open (10, file = fnin(j), form = 'formatted', status = 'old')
      !
      !    skip header
      !    -----------
      print *,'a'
      do while (.true.)
        read (10, '(a1)') flag
        if (flag .ne. '#') then
          backspace (10)
          exit
        end if
      end do
      !
      !    count data
      !    ----------
      i = 1
      print *,'b'
      do while (.true.)
        read (10, *, iostat = iocheck) tdum, xdum
        if (i .eq. 1) t1 = tdum        ! save initial time
        if (iocheck .ne. 0) exit
        i = i + 1
      end do
      close(10)
      !
      !    number of input data
      !    --------------------
      !     
      print *,'c'
      if (j == 1) then
        npx = i - 1              !number of data points x
        avgdtx = (tdum - t1) / real(npx-1, dp)          ! avg. sampling interval for x
        print '(a,f11.3,a,f11.3,a,i8,a)', 'X time interval :   [ ',t1,'; ',tdum,' ]       - ',npx,' points'
      else if (j == 2) then
        npy = i -1              !number of data points y
        avgdty = (tdum - t1) / real(npy-1, dp)        ! avg. sampling interval for y
        print '(a,f11.3,a,f11.3,a,i8,a)', 'Y time interval :   [ ',t1,'; ',tdum,' ]       - ',npy,' points'
      end if
      !
    end do   !(end of j - loop)  
    !
    ! number of output data (results saved in the module nyquist)
    ! ----------------------------------------------------------
    nsegx = int(2 * npx / (n50 + 1))             ! points per segment x
    nsegy = int(2 * npy / (n50 + 1))             ! points per segment  y
    !  
    if (avgdtx >= avgdty) then
      avgdtxy = avgdtx
    else
      avgdtxy = avgdty              ! eq 7 in the SPECTRUM paper;avgdtxy = max (avgdtx,avgdty)
    end if
    tpx = avgdtx * nsegx                          ! average period of a segment x
    tpy = avgdty * nsegy                          ! average period of a segment y
    dfx = 1 / (ofac * tpx)                        ! freq. spacing x  (fundamental freq)
    dfy = 1 / (ofac * tpy)                        ! freq. spacing y
    if ( dfx >= dfy) then 
      dfxy = dfx                !eq 8 in SPECTRUM paper dfxy = max (dfx,dfy)
    else 
      dfxy = dfy
    end if
    wz = 2 * pi * dfxy            ! omega = 2*pi*f
    fnyq = hifac * 1 / (2 * avgdtxy)       ! average Nyquist freq.
    nfreq = nint(fnyq / dfxy + 1)                     ! f(1) = f0; f(nfreq) = fNyq
    lfreq = nfreq * 2
    !
    ! diagnostic output to stdout
    ! ---------------------------
    print '(a)'
    print '(a,f10.2)', "dtxy =",avgdtxy
    print '(a,i10)', "Nfreq =", nfreq
    print '(a,f10.2)', "fnyq =",fnyq
    print '(a)'
    !
  end subroutine setdim

  !--------------------------------------------------------------------------
  subroutine readdat(fnin)
    !--------------------------------------------------------------------------
    use timeser
    use const,only:n_fnin
    !
    implicit none
    !
    character (len = 80), dimension (n_fnin), intent(in) :: fnin 
    integer :: i, iocheck, j
    character (len = 1) :: flag
    !
    do j = 1,n_fnin    !Goes through inputfile x and y
      !
      !    open input file
      !    ---------------
      open (10, file = fnin(j), form = 'formatted', status = 'old')
      !
      !    skip header
      !    -----------
      do while (.true.)
        read (10, '(a1)') flag
        if (flag .ne. '#') then
          backspace (10)
          exit
        end if
      end do
      !
      !    retrieve data
      !    -------------
      if (j ==1 ) then
        do i = 1, npx
          read (10,*) tx(i), x(i)
        end do
        close(10)
      else if (j == 2) then
        do i = 1, npy
          read (10,*) ty(i), y(i)
        end do
        close(10)
      end if
      !
    end do
    !
  end subroutine readdat

  !--------------------------------------------------------------------------
  subroutine spectr(ini, tx, x, ty, y, ofac, n50, iwin, frq, gxx, &
    gyy, gxy, cxy, phxy)
    !--------------------------------------------------------------------------
    use trigwindat
    use const
    use nyquist
    !
    implicit none
    !
    real(dp), parameter :: si = 1
    real(dp), parameter :: tzero = 0
    !
    logical, intent(in) :: ini
    real(dp), dimension(:), intent(in) :: tx, x, ty, y
    real(dp), intent(in) :: ofac
    integer, intent(in) :: n50, iwin
    real(dp), dimension(:), intent(out) :: frq, gxx, gyy, gxy, cxy, phxy  
    !
    real(dp), dimension(:), allocatable :: txwk, xwk, ftrx, ftix, &  
      tywk, ywk, ftry, ftiy                                         
    integer :: i, j, istart, ialloc          
    real(dp) ::  scalx,scaly,scalxy      
    real(dp) :: rnx, rny  
    complex(dp),dimension(nfreq) :: cpxy  
    !   
    !
    ! setup workspace
    ! ---------------
    gxx(:) = 0  
    gyy(:) = 0
    gxy(:) = 0
    cxy(:) = 0
    phxy(:)= 0
    cpxy(:)= 0
    allocate(txwk(nsegx), xwk(nsegx), tywk(nsegy), ywk(nsegy))
    allocate(ftrx(lfreq), ftix(lfreq), ftry(lfreq), ftiy(lfreq))
    if (ini .eqv. .true.) then
      allocate(txcos(nsegx,nfreq,n50),tycos(nsegy,nfreq,n50))
      allocate(txsin(nsegx,nfreq,n50),tysin(nsegy,nfreq,n50))
      allocate(wxtau(nfreq,n50),wytau(nfreq,n50))
      allocate(wwx(nsegx,n50), wwy(nsegy,n50))
    end if
    !

    do i = 1, n50    ! start of the segment loop
      !    copy data of i'th segment into workspace
      !    ----------------------------------------
      istart = (i-1) * nsegx / 2
      do j = 1, nsegx
        txwk(j) = tx(istart + j)
        xwk(j) = x(istart + j)
      end do
      !     
      istart = (i-1) * nsegy / 2
      do j = 1, nsegy
        tywk(j) = ty(istart + j)
        ywk(j) = y(istart + j)
      end do
      !
      !    detrend data
      !    ------------
      call rmtrend (txwk(1:nsegx), xwk(1:nsegx), nsegx)
      call rmtrend (tywk(1:nsegy), ywk(1:nsegy), nsegy)
      !
      !    apply window to data
      !    --------------------
      if (ini .eqv. .true.) call winwgt(txwk(1:nsegx), iwin, wwx(1:nsegx, i))
      xwk(1:nsegx) = wwx(1:nsegx, i) * xwk(1:nsegx) 
      if (ini .eqv. .true.) call winwgt(tywk(1:nsegy), iwin, wwy(1:nsegy, i))
      ywk(1:nsegy) = wwy(1:nsegy, i) * ywk(1:nsegy)
      !
      !    setup trigonometric array for LSFT
      !    ----------------------------------
      if (ini .eqv. .true.) &
        call trig(i, txwk(1:nsegx), nsegx, wz, nfreq, txcos,txsin, wxtau)
      if (ini .eqv. .true.) &
        call trig(i, tywk(1:nsegy), nsegy, wz, nfreq, tycos,tysin, wytau)
      !
      !    LSFT
      !    ------
      call ftfix(i, xwk(1:nsegx), txwk(1:nsegx), nsegx, &
        txcos(1:nsegx,1:nfreq,1:n50), txsin(1:nsegx,1:nfreq,1:n50), &
        wxtau(1:nfreq,1:n50), wz, nfreq, si, lfreq, tzero, ftrx, ftix)
      !
      call ftfix(i, ywk(1:nsegy), tywk(1:nsegy), nsegy, &
        tycos(1:nsegy,1:nfreq,1:n50),tysin(1:nsegy,1:nfreq,1:n50), &
        wytau(1:nfreq,1:n50), wz, nfreq, si, lfreq, tzero, ftry, ftiy)
      !
      !    sum raw spectra
      !    -------------------
      do j = 1, nfreq
        gxx(j) = gxx(j) + (ftrx(j)*ftrx(j) + ftix(j)*ftix(j))
        gyy(j) = gyy(j) + (ftry(j)*ftry(j) + ftiy(j)*ftiy(j))
      end do
      !
      !    cross and phase spectra
      !    --------------------------
      do j = 1,nfreq
        cpxy(j) = cpxy(j) + &
          cmplx(ftrx(j), ftix(j), dp) * cmplx(ftry(j), -ftiy(j), dp)
        !the minus obtains the complex conjugate
        phxy(j) = phxy(j) + atan2(imag(cpxy(j)),real(cpxy(j))) * 180/pi
        ! see how atan2 function works atan2(y,x) = atan(y/x)
      end do    ! which means atan(imag(cpxy)/real(cpxy))
      !
    end do      ! end the segment loop
    !     
    ! scale autospectrum and setup frequency axis
    ! -------------------------------------------
    !
    ! determine smoothed spectral estimate, i.e calc. average over k spectra and
    ! scale spectrum such that the intergral of the smoothed spectrum equals the
    ! data variance, that is a (Gxx(f)*df) = a^2 ; change sign of the phase
    ! because of the geological time axis is inverted compared to the physical
    ! axis of time
    !
    scalx = 2 / (n50 * nsegx * dfxy * ofac)
    scaly = 2 / (n50 * nsegy * dfxy * ofac)
    rnx = nsegx           ! make the integer nsegx real number
    rny = nsegy
    scalxy = 2 / (n50 * sqrt(rnx * rny) * dfxy * ofac)
    do i = 1, nfreq
      gxx(i) = gxx(i) * scalx
      gyy(i) = gyy(i) * scaly
      gxy(i) = abs(cpxy(i)) * scalxy
      phxy(i) = phxy(i) / n50
      frq(i) = (i-1) * dfxy
    end do
    ! 
    ! coherency spectrum 
    ! ------------------
    !  mask exeptions in order to avoid 'division by zero' and 'overflow' errors. 
    !  If only 1 segment was used then set cxy(f) = 1 -> faster
    !
    if (n50 == 1) then
      cxy(:) = 1
    else   
      do i=1,nfreq
        if(gxy(i)==0 .or. gxx(i)==0 .or. gyy(i)==0) then
          cxy(i) = 0
        else 
          cxy(i)= (gxy(i) * gxy(i))/(gxx(i) * gyy(i))
        end if
      end do
    end if
    !
    deallocate(txwk, xwk, tywk, ywk)
    deallocate(ftrx, ftix, ftry, ftiy)
    !
  end subroutine spectr

  !--------------------------------------------------------------------------
  subroutine makear1(t, np, tau, red)
    !--------------------------------------------------------------------------
    use const
    use random
    !
    implicit none
    !
    integer, intent(in)    :: np
    real(dp), dimension(np), intent(in) :: t
    real(dp), intent(in)       :: tau
    real(dp), dimension(np), intent(out) :: red
    !
    real(dp) :: sigma, dt  
    integer :: i
    real(dp) ::  z1
    !
    ! set up AR(1) time series
    ! ------------------------

    if (tau.eq.0) then
      do i = 1,np
        call gasdev(z1)
        red(i) = z1
      end do
    else
      call gasdev(z1)
      red(1)=z1
      do i = 2, np
        dt = t(i) - t(i-1)
        sigma = 1 - exp(-2 * dt / tau)
        sigma = sqrt (sigma)
        call gasdev(z1)
        red(i) = exp(-dt/tau) * red(i-1) + sigma * z1   
      end do
    end if
    ! 
  end subroutine makear1

  !--------------------------------------------------------------------------
  subroutine make_coherar1(red1, red2)
    !--------------------------------------------------------------------------
    use phase
    use random
    !
    implicit none
    !
    ! Generates two coupled time series with prescribed coherency
    !
    real(dp), dimension(np_xy), intent(out) :: red1
    real(dp), dimension(np_xy), intent(out) :: red2
    !
    real(dp) :: sigmax, sigmay, dtx, dty  
    integer :: i, j
    real(dp) ::  z1, z2
    real(dp), dimension(np_xy) :: ex, ex_y
    real(dp), dimension(np_xy) :: ey, ey_x
    !
    !Set up noise x
    !--------------------
    do i = 1,np_xy
      call gasdev(z1)
      ex(i) = z1
    end do
    !
    !  Set up noise y
    !---------------------
    do i = 1,np_xy
      call gasdev(z2)
      ey(i) = z2
    end do
    !
    !  Couple the noise terms
    !----------------------------
    ex_y(:) = ex(:) + g*ey(:)
    ey_x(:) = ey(:) + g*ex(:)
    !
    ! Set up AR(1) series x
    !----------------------------
    !  
    if (tau_xy .le.0) then
      red1(:) = ex_y(:)
    else
      red1(1) = ex_y(1)
      do j = 2,np_xy
        dtx = t_xy(j) - t_xy(j-1)
        sigmax = 1 - exp(-2 * dtx / tau_xy)
        sigmax = sqrt (sigmax)
        red1(j) = exp(-dtx/tau_xy) * red1(j-1) + sigmax * ex_y(j)
      end do
    end if
    !
    ! Set up AR(1) series y
    !  
    if (tau_xy .le.0) then
      red2(:) = ey_x(:)
    else
      red2(1) = ey_x(1)
      do j = 2,np_xy
        dty = t_xy(j) - t_xy(j-1)
        sigmay = 1 - exp(-2 * dty / tau_xy)
        sigmay = sqrt (sigmay)
        red2(j) = exp(-dty/tau_xy) * red2(j-1) + sigmay * ey_x(j)
      end do
    end if
    !  
  end subroutine make_coherar1

  !-------------------------------------------------------------------
  subroutine ftfix(iseg, xx, tsamp, nn, tcos, tsin, wtau, wz, nfreq, &
    si, lfreq, tzero, ftrx, ftix)
    !-----------------------------------------------------------------
    ! Fourier transformation for unevenly spaced data
    ! (Scargle, 1989; ApJ 343, 874-887)
    !
    ! - folding of trigonom. and exp. arguments in a*pi disabled
    !-----------------------------------------------------------------
    !
    implicit  none
    !
    real(dp), dimension(:), intent(in) :: xx, tsamp
    real(dp), dimension(:), intent(out) :: ftrx, ftix
    real(dp), intent(in)    :: wz, si, tzero
    integer, intent(in) :: iseg, nn, nfreq, lfreq
    real(dp), dimension(:,:,:), intent(in) :: tcos, tsin
    real(dp), dimension(:,:), intent(in) :: wtau
    !
    !
    real(dp), parameter :: tol1 = 1e-4_dp
    real(dp), parameter :: tol2 = 1e-8_dp
    real(dp), parameter :: sqrt2= 1.41421356237309504880168872420969807856967_dp
    real(dp), parameter :: const1 = 1/sqrt2
    !
    real(dp)    :: const2, wdel, wrun, wuse, fnn, ftrd, ftid, phase, &
      sumt, sumx, sumr, sumi, scos2, ssin2, cross, wtnew
    integer :: i, i1, ii, iput, istop, nstop
    complex(dp) :: work
    !
    wuse = wz
    fnn  = float(nn)
    const2 = si * const1
    sumt = sum (tsamp(1:nn))
    sumx = sum (xx(1:nn))
    istop = nfreq
    !
    ! initialize for zero frequency
    ! -----------------------------
    ftrx(1) = sumx / sqrt(fnn)
    ftix(1) = 0
    wdel = wuse
    wrun = wuse
    II = 2
    !
    ! start frequency loop
    ! --------------------
    do while (.true.)
      wtnew = wtau(ii,iseg)
      !
      !   summations over the sample
      !   --------------------------
      cross = sum(tsamp(1:nn) * tcos(1:nn,ii,iseg) * tsin(1:nn,ii,iseg))
      scos2 = sum(tcos(1:nn,ii,iseg)**2)
      ssin2 = sum(tsin(1:nn,ii,iseg)**2)
      sumr = sum(xx(1:nn) * tcos(1:nn,ii,iseg))
      sumi = sum(xx(1:nn) * tsin(1:nn,ii,iseg))
      !
      ftrd = const1 * sumr / sqrt(scos2)
      if (ssin2 .le. tol1) then
        ftid = const2 * sumx / sqrt(fnn)
        if (abs(cross) .gt. tol2) ftid = 0
      else
        ftid = const2 * sumi / sqrt(ssin2)
      end if
      phase = wtnew - wrun * tzero
      work = cmplx(ftrd, ftid, dp) * exp(cmplx(0.0_dp, phase, dp))
      ftrx(ii) = real(work, dp)
      ftix(ii) = aimag(work)
      ii = ii + 1
      wrun = wrun + wdel
      if (ii .gt. istop) exit
    end do
    !
    ! zero-fill transform (oversample inverse) impose symmetry for real data
    ! ----------------------------------------------------------------------
    if (2 * nfreq .gt. lfreq) then
      write (*,*) 'Error: 2 * nfreq > lfreq'
      stop
    end if
    i1 = nfreq + 1
    do i = i1, lfreq
      ftrx(i) = 0
      ftix(i) = 0
    end do
    nstop = lfreq / 2
    do i = 2, nstop
      iput = lfreq - i + 2
      ftrx(iput) =  ftrx(i)
      ftix(iput) = -ftix(i)
    end do
    !
  end subroutine ftfix

  !------------------------------------------------------------
  subroutine trig(iseg, tsamp, nn, wz, nfreq, tcos, tsin, wtau)
    !----------------------------------------------------------
    use param
    !
    implicit  none
    !
    real(dp), dimension(*), intent(in) :: tsamp
    real(dp), intent(in)    :: wz
    integer, intent(in) :: iseg, nn, nfreq
    real(dp),dimension(:,:,:),intent(out):: tcos, tsin
    real(dp),dimension(:,:),intent(out):: wtau
    !
    real(dp), parameter :: tol1 = 1e-4_dp
    real(dp)    :: csum, ssum, wdel, wrun, wuse,  tim, sumtc, sumts, ttt, &
      tc, ts, watan, wtnew,arg
    integer :: i, ii, istop
    !
    wuse = wz
    istop = nfreq
    wdel = wuse
    wrun = wuse
    ii = 2
    !
    ! start frequency loop
    ! --------------------
    do while (.true.)
      !
      !   calc. tau
      !   ---------
      csum = 0
      ssum = 0
      sumtc = 0
      sumts = 0
      do i = 1, nn
        ttt  = tsamp(i)
        arg = 2 * wrun * ttt
        tc = cos(arg)
        ts = sin(arg)
        csum = csum + tc
        ssum = ssum + ts
        sumtc = sumtc + ttt * tc
        sumts = sumts + ttt * ts
      end do
      if (abs(ssum) .gt. tol1 .or. abs(csum) .gt. tol1) then
        watan = atan2 (ssum, csum)
      else
        watan = atan2 (-sumtc, sumts)
      end if
      wtau(ii,iseg) = 0.5_dp * watan
      wtnew = wtau(ii,iseg)
      !
      !   summations over the sample
      !   --------------------------
      do i = 1, nn
        tim  = tsamp(i)
        arg = wrun * tim - wtnew
        tcos(i,ii,iseg) = cos(arg)
        tsin(i,ii,iseg) = sin(arg)
      end do
      ii = ii + 1
      wrun = wrun + wdel
      if (ii .gt. istop) exit
    end do
    !
  end subroutine trig

  !------------------------------------------
  real(dp) function fold(arg)
    !----------------------------------------
    use const, only : pi
    !
    implicit  none
    !
    real(dp), intent(in) :: arg
    !
    real(dp), parameter :: argmax = 8000 * pi
    !
    fold = arg
    do while (.true.)
      if (fold .le. argmax) exit
      fold = fold - argmax
    end do
    do while (.true.)
      if (fold .gt. -argmax) exit
      fold = fold + argmax
    end do
    !
  end function fold

  !------------------------------------------------
  subroutine winwgt(t, iwin, ww)
    !----------------------------------------------
    ! calc. normalized window weights
    ! window type (iwin)  0: Rectangular
    !                     1: Welch 1
    !                     2: Hanning
    !                     3: Parzen (Triangular)
    !                     4: Blackman-Harris 3-Term
    !----------------------------------------------
    use const
    !
    implicit none
    !
    real(dp), intent(in), dimension(:) :: t
    real(dp), intent(out), dimension(:) :: ww
    integer, intent(in) :: iwin
    !
    real(dp) :: tlen, jeff, scal, fac1, fac2, fac3, fac4, sumw2, rnp
    integer :: i, nseg
    !
    ! useful factor for various windows
    ! ---------------------------------
    nseg = size(t)
    rnp = real(nseg, dp)
    fac1 = (rnp / 2 ) - 0.5_dp
    fac2 = 1 / ((rnp / 2 ) + 0.5_dp)
    fac3 = rnp - 1
    fac4 = tpi /(rnp - 1)
    tlen = t(nseg) - t(1)
    sumw2 = 0
    do i= 1, nseg
      jeff = rnp * (t(i)-t(1)) / tlen
      select case (iwin)
      case (0)                                            ! rectangle
        ww(i) = 1
      case (1)                                            ! welch I
        ww(i) = 1 - ((jeff - fac1) * fac2)**2
      case (2)                                            ! hanning
        ww(i) = 0.5_dp * (1 - cos(tpi * jeff / fac3))
      case (3)                                            ! triangular
        ww(i) = 1 - abs((jeff - fac1) * fac2)
      case (4)                                            ! blackman-harris
        ww(i) = 0.4243801_dp - 0.4973406_dp * cos(fac4 * jeff) &
          + 0.0782793_dp * cos(fac4 * 2 * jeff)
      end select
      sumw2 = sumw2 + ww(i) * ww(i)
    end do
    !
    ! determine scaling factor and scale window weights;
    ! NB: sumw2 = nseg for rectangular window
    ! --------------------------------------------------
    scal = sqrt(rnp / sumw2)
    do i = 1, nseg
      ww(i) = ww(i) * scal
    end do
    !
  end subroutine winwgt

  !--------------------------------------------------------------------------
  function winbw(iwin, df, ofac)
    !--------------------------------------------------------------------------
    ! Determine 6dB bandwidth from OFAC corrected fundamental frequency.
    ! Note that the bandwidth for the Blackman-Harris taper is higher than
    ! reported by Harris (1978, cf. Nuttall, 1981)}
    !
    ! window type (iwin)  0: Rectangular
    !                     1: Welch 1
    !                     2: Hanning
    !                     3: Parzen (Triangular)
    !                     4: Blackman-Harris 3-Term
    !--------------------------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: iwin
    real(dp), intent(in)    :: df, ofac
    real(dp)                :: winbw
    !
    real(dp), parameter, dimension(0:4) :: bw = &
      [1.21_dp, 1.59_dp, 2.00_dp, 1.78_dp, 2.26_dp]
    !
    winbw = df * ofac * bw(iwin)
    !
  end function winbw

  !------------------------------------------------------------------------
  subroutine rmtrend (x, y, n)
    !------------------------------------------------------------------------
    ! determine linear trend by means of least squares and subtract
    ! this trend from a data set
    !
    ! parameters:  x, y   : real arrays for data, on output y is replaced
    !                       by the detrended values
    !              n      : number of data pairs
    !
    ! ref.: after numerical recipes 14.2
    !
    ! written:       23/07/94
    ! modifications: 29/01/99 - allow n <> array dimension
    !------------------------------------------------------------------------
    implicit  none
    !
    real(dp), dimension(*), intent(in) :: x
    real(dp), dimension(*), intent(inout) :: y
    integer, intent(in) :: n
    !
    real(dp)    :: sx, sy, st2, a, b, ss, sxoss, z
    integer :: i
    !
    sx = 0
    sy = 0
    st2 = 0
    b = 0
    do i = 1, n
      sx = sx + x(i)
      sy = sy + y(i)
    end do
    ss = float(n)
    sxoss = sx / ss
    do i = 1, n
      z = x(i) - sxoss
      st2 = st2 + z*z
      b = b + z * y(i)
    end do
    b = b / st2
    a = (sy - sx * b) / ss
    do i = 1, n
      y(i) = y(i) - (a + b * x(i))
    end do
    !
  end subroutine rmtrend

  !----------------------------------------------------------------------
  subroutine getdof(iwin, n50, dof, neff)
    !--------------------------------------------------------------------------
    ! Effective number of degrees of freedom for the selected window
    ! and n50 overlappaing segments (Harris, 1978)
    !----------------------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: iwin, n50
    real(dp), intent(out) :: dof,neff
    !
    real(dp), parameter, dimension(0:4) :: c50 = &
      [0.500_dp, 0.344_dp, 0.167_dp, 0.250_dp, 0.096_dp]
    real(dp) :: c2, denom, rn
    !
    rn = real(n50, dp)  
    c2 = 2 * c50(iwin) * c50(iwin)

    denom = 1 + c2 - c2/rn
    neff = rn / denom
    dof = 2 * neff
    !
  end subroutine getdof

  !----------------------------------------------------------------------
  real(dp) function getchi2(dof, alpha)
    !----------------------------------------------------------------------
    use asa239
    implicit none
    !
    real(dp), parameter :: tol = 1e-3_dp
    integer, parameter :: itmax = 100
    real(dp) :: dof, alpha
    real(dp) :: ac, lm, rm, eps, chi2, za, x
    integer :: iter, ifault
    !
    ! use approximation for dof > 30 (Eq. 1.132 in Sachs (1984))
    ! ----------------------------------------------------------
    if (dof .gt. 30) then
      za = -getz(alpha)   ! NB: Eq. requires change of sign for percentile
      x = 2.0_dp / 9.0_dp / dof
      chi2 = dof * (1 - x + za * sqrt(x))**3
    else
      iter = 0
      lm = 0
      rm = 1000
      if (alpha .gt. 0.5_dp) then
        eps = (1 - alpha) * tol
      else
        eps = alpha * tol
      end if
      do
        iter= iter + 1
        if (iter .gt. itmax) then
          stop "Error in GETCHI2: Iter > ItMax"
        end if
        chi2 = 0.5_dp * (lm + rm)
        ac = 1 - gammad(0.5_dp*chi2, 0.5_dp*dof, ifault)
        if (ifault /= 0) stop "gammad failure"
        if (abs(ac - alpha) .le. eps) exit
        if (ac .gt. alpha) then
          lm = chi2
        else
          rm = chi2
        end if
      end do
    end if
    getchi2 = chi2
    !
  end function getchi2

  !----------------------------------------------------------------------
  real(dp) function getz(alpha)
    !----------------------------------------------------------------------
    ! Determine percentiles of the normal distribution using an approximation
    ! of the complementary error function by a Chebyshev polynom.
    !
    ! For a given values of alpha (a), the program returns z(a) such that
    ! P[Z <= z(a)] = a. Check values are in the front cover of Neter et al.
    !----------------------------------------------------------------------
    implicit none
    !
    real(dp), parameter :: tol = 1e-5_dp
    real(dp), parameter :: sq2 = 1.414213562_dp
    integer, parameter :: itmax = 100
    real(dp) :: alpha
    real(dp) :: atmp, acalc, zr, zl, zm, z
    integer :: iter
    !
    if (alpha .lt. 0.5_dp) then
      atmp = alpha * 2
      zr = -0.1_dp
      zl = 10
      iter = 0
      do while(.true.)
        iter= iter + 1
        if (iter .gt. itmax) then
          print *, "Error in GETZ: Iter > ItMax"
          return
        end if
        zm = (zl + zr) / 2
        z = zm
        acalc = erfc(z/sq2)
        if (acalc .gt. atmp) zr = zm
        if (acalc .le. atmp) zl = zm
        if (abs(acalc-atmp) .le. tol) exit
      end do
      z = -1 * z
    else if (alpha .ge. 0.5_dp) then
      atmp =(alpha - 0.5_dp) * 2
      zl = -0.1_dp
      zr = 10
      iter = 0
      do while(.true.)
        iter= iter + 1
        if (iter .gt. itmax) then
          write(*,*) "Error in GETZ: Iter > ItMax"
          return
        end if
        zm = (zl + zr) / 2
        z = zm
        acalc = 1 - erfc(zm/sq2)
        if (acalc .gt. atmp) zr = zm
        if (acalc .le. atmp) zl = zm
        if (abs(acalc-atmp) .le. tol) exit
      end do
    end if
    getz = z
    !
  end function getz

  !----------------------------------------------------------------------
  subroutine gettau(tx,x,npx,tau)
    !----------------------------------------------------------------------
    use const
    use param, only : n50
    use minls_module
    implicit none
    real(dp), dimension(:),intent(in) :: tx, x
    integer,intent(in) :: npx
    real(dp), intent(out) :: tau
    real(dp), dimension(:), allocatable :: twk, xwk
    integer :: nseg, i, j, istart, ialloc
    real(dp) :: rho, rhosum, avgdt
    !
    ! average dt of entire time series
    ! --------------------------------
    avgdt = sum(tx(2:npx)-tx(1:npx-1)) / real(npx-1, dp)
    !
    rhosum = 0
    nseg = int(2 * npx / (n50 + 1))         ! points per segment
    allocate(twk(nseg), xwk(nseg))
    do i = 1, n50
      !
      !       copy data of i'th segment into workspace
      !       ----------------------------------------
      istart = (i-1) * nseg / 2
      do j = 1, nseg
        twk(j) = tx(istart + j)
        xwk(j) = x(istart + j)
      end do
      !
      !        detrend data
      !        ------------
      call rmtrend (twk(1:nseg), xwk(1:nseg), nseg)
      !
      !        estimate and sum rho for each segment
      !        -------------------------------------
      call tauest_x(twk(1:nseg), xwk(1:nseg), nseg, tau, rho)
      !
      !        bias correction for rho (Kendall & Stuart, 1967; Vol. 3))
      !        ---------------------------------------------------------
      rho = (rho * (real(nseg, dp) - 1) + 1) / (real(nseg, dp) - 4)
      !
      rhosum = rhosum + rho
    end do
    !
    !     average rho
    !     -----------
    rho = rhosum / real(n50, dp)
    !
    !     average tau
    !     -----------
    tau = -avgdt / log(rho)
    !
    !     make sure that tau is non-negative
    !     ----------------------------------
    if (tau .lt. 0) then
      print *, 'Warning: GETTAU returned tau =', tau
      print *, '         Negative tau is forced to zero.'
      tau = 0
    end if
    !
    deallocate(twk, xwk)
    ! 
  end subroutine gettau

  !----------------------------------------------------------------------
  !  Manfred Mudelsee's code for tau estimation
  !----------------------------------------------------------------------
  ! TAUEST: Routine for persistence estimation for unevenly spaced time series
  !----------------------------------------------------------------------
  !       Main variables
  !
  !       t       :       time
  !       x       :       time series value
  !       np      :        number of points
  !      dt       :       average spacing
  !   scalt       :       scaling factor (time)
  !     rho       :       in the case of equidistance, rho = autocorr. coeff.
  !      ls       :       LS function
  !   brent       :       Brent's search, minimum LS value
  !    mult       :       flag (multiple solution)
  !    amin       :       estimated value of a = exp(-scalt/tau)
  !
  !----------------------------------------------------------------------
  ! subroutine tauest    !
  !   use const
  !   use meanvar_module
  !   !
  !   implicit none
  !   !
  !   integer, intent(in) :: np
  !   real(dp), dimension(np), intent(in) :: t, x
  !   real(dp), intent (out)  :: tau
  !   real(dp), dimension(np) :: tscal, xscal
  !   real(dp) :: fac, avg, var, dt, rho, scalt, amin, rhoavg
  !   real(dp) :: damin
  !   integer :: i, mult
  !   !
  !   ! Correct time direction; assume that ages are input
  !   ! --------------------------------------------------
  !   do i = 1, np
  !     tscal(i) = -t(np+1-i)
  !     xscal(i) = x(np+1-i)
  !   end do
  !   !
  !   ! Scaling of x
  !   ! ------------
  !   call meanvar(xscal(1:np), avg, var)
  !   fac = sqrt(var)
  !   xscal(1:np) = xscal(1:np) / fac
  !   !
  !   ! Scaling of t (=> start value of a = 1/e)
  !   ! ---------------------------------------
  !   dt = (tscal(np)-tscal(1)) / real(np-1, dp)
  !   call rhoest(np, xscal(1:np), rho)
  !   if (rho .le. 0) then
  !     rho = 0.05_dp
  !     write(errio,*) 'Warning: rho estimation: < 0'
  !     ierr = 2
  !   else if (rho .gt. 1) then
  !     rho = 0.95_dp
  !     write(errio,*) 'Warning: rho estimation: > 1'
  !     ierr = 2
  !   end if
  !   scalt = -log(rho)/dt
  !   tscal(1:np) = tscal(1:np) * scalt
  !   !
  !   ! Estimation
  !   ! ----------
  !   call minls(np, tscal(1:np), xscal(1:np), damin, mult)
  !   if (ierr .eq. 1) then
  !     write(errio,*) ' Error in MNILS'
  !     return
  !   end if
  !   amin = sngl(damin)
  !   if (mult .eq. 1) then
  !     write(errio,*) ' Estimation problem: LS function has > 1 minima'
  !     return
  !   end if
  !   if (amin .le. 0) then
  !     write(errio,*) ' Estimation problem: a_min =< 0'
  !     return
  !   else if (amin .ge. 1) then
  !     write(errio,*) ' Estimation problem: a_min >= 1'
  !     return
  !   end if
  !   !
  !   ! determine tau
  !   ! -------------
  !   tau = -1 /(scalt*log(amin))
  !   !
  !   ! determine rho, corresponding to tau
  !   ! -----------------------------------
  !   rhoavg = exp(-dt / tau)
  !   !
  ! end subroutine tauest

  ! Numerical Recipes (modified): Brent's search in one direction:
  ! function brent(ax,bx,cx,f,tol,xmin,xfunc,yfunc,nfunc)

  !----------------------------------------------------------------------
  ! Least-squares function
  !----------------------------------------------------------------------
  ! real(dp) function ls(a,t,x,n)
  !   implicit none
  !   real(dp), intent(in) :: a, t(:), x(:)
  !   integer, intent(in) :: n
  !   integer i
  !   ls=0.0_dp
  !   do i=2,n
  !     ls=ls+(x(i)-x(i-1)*sign(1.0_dp,a)* abs(a)**(t(i)-t(i-1)))**2
  !   end do
  !   return
  ! end function ls

  !----------------------------------------------------------------------
  ! Minimization of least-squares function ls.
  !----------------------------------------------------------------------
  ! subroutine minls(n, t, x, amin, nmu_)
  !   !
  !   use newbrent_module
  !   !
  !   implicit none
  !   !
  !   real(dp), parameter :: a_ar1 = 0.367879441_dp ! 1/e
  !   real(dp), parameter :: tol = 3e-8_dp           ! Brent's search, precision
  !   real(dp), parameter :: tol2 = 1e-6_dp          ! multiple solutions, precision
  !   integer n
  !   real(dp) t(1:n),x(1:n)
  !   real(dp) amin
  !   integer nmu_
  !   real(dp) dum1,dum2,dum3,dum4,a_ar11,a_ar12,a_ar13
  !   !
  !   nmu_=0
  !   dum1=newbrent(ls, -2.0_dp, a_ar1, +2.0_dp, tol, a_ar11, t, x, n)
  !   dum2=newbrent(ls, a_ar1, 0.5_dp*(a_ar1+1.0_dp), +2.0_dp, tol, a_ar12, t, x, n)
  !   dum3=newbrent(ls,-2.0_dp, 0.5_dp*(a_ar1-1.0_dp),  a_ar1, tol, a_ar13, t, x, n)
  !   if (ierr .eq. 1) then
  !     write(errio, *) ' Error in MINLS (call to brent)'
  !     return
  !   end if
  !   if  ((abs(a_ar12-a_ar11).gt.tol2.and.abs(a_ar12-a_ar1).gt.tol2) &
  !     .or.(abs(a_ar13-a_ar11).gt.tol2.and.abs(a_ar13-a_ar1).gt.tol2)) &
  !     nmu_=1
  !   dum4=min(dum1,dum2,dum3)
  !   if (dum4.eq.dum2) then
  !     amin=a_ar12
  !   else if (dum4.eq.dum3) then
  !     amin=a_ar13
  !   else
  !     amin=a_ar11
  !   end if
  !   return
  ! end subroutine minls

  !----------------------------------------------------------------------
  ! Autocorrelation coefficient estimation (equidistant data).
  !----------------------------------------------------------------------
  ! subroutine rhoest(n,x,rho)
  !   !
  !   implicit none
  !   !
  !   integer n
  !   real(dp) x(1:n)
  !   real(dp) rho
  !   integer i
  !   real(dp) sum1,sum2
  !   !
  !   sum1=0
  !   sum2=0
  !   do i=2,n
  !     sum1=sum1+x(i)*x(i-1)
  !     sum2=sum2+x(i)**2
  !   end do
  !   rho=sum1/sum2
  !   return
  ! end subroutine rhoest

end module redfit_x_module
