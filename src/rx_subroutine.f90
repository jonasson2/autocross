! Subroutine interface to RedfitX. See comments in redfitx.f90

subroutine rx_subroutine(npx, npy, nout, tx, ty, x, y, &
  &  cfg_nsim, cfg_rhopre, cfg_ofac, cfg_hifac, cfg_n50, cfg_alpha, cfg_i, &
  &  rhox, rhoy, taux, tauy, df, dB6, false_alarm, &
  &  scale, data_x, data_y, data_xy, data_cxy, data_phxy)
  use precision
  use const
  use timeser
  use param
  use trigwindat
  use nyquist
  use phase
  use random
  use meanvar_module
  use sort_module
  use redfit_x_module

  implicit none
  
  integer, intent(in) :: npx, npy, nsim, n50, iwin
  real(dp), intent(in) :: tx(npx), x(npx), ty(npy), y(npy)
  real(dp), intent(in) :: rhopre(2), ofac, hifac, alpha
  real(dp), intent(out) :: rhox, rhoy, taux, tauy, dof, &
    dB6, false_alarm, faccritx
  real(dp), intent(out) :: data_x(nout, 12), data_y(nout, 12), &
    data_xy(nout, 2), data_cxy(nout, 7), data_phxy(nout, 6)
 
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
  real(dp) :: rnsim, fac90, fac95, fac99, neff, &
    avgdtx, avgdty, facx, facy, facxy, rhoxsq, rhoysq, varx, vary, &
    varrx, varry, alphacritx, alphacrity, faccritx, faccrity,varxy, &
    varrxy, cobias, facphi, z, csig, se_dummy
  integer:: kstart, kstop, krate, kmax, ntime
  integer :: i, iocheck, ialloc
  integer :: idx90, idx95, idx99, iostat
  logical :: ini, biascorr
  !
  call system_clock(kstart, krate, kmax)
  !
  ! First try to direct error messages to stderr
  ! --------------------------------------------
  errorfile = '/dev/stderr'
  open(errio, file=errorfile, status='old', iostat=iostat, action='write')

  ! call check1()  Tékka hvort röðin er vaxandi eða minnkandi
  !                (hún á að vera vaxandi) (x og y sér)
  ! 
  ! in case of duplicate entries reinitialize array dimensions
  ! and retrieve averaged time series
  ! ----------------------------------------------------------

  ! allocate remaining workspace
  ! ----------------------------
  allocate(redx(npx), redy(npy), stat = ialloc)
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
  allocate (ci90(nspect,nout), ci95(nspect,nout), &  
    ci99(nspect,nout), stat = ialloc)
  if (ialloc .ne. 0) call allocerr("a")
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
  call gettau(rhopre(2), ty(1:npy), y(1:npy), npy, tauy)
  !
  ! Generate NSim AR(1) Spectra
  ! ---------------------------------
  !
  allocate(grxx(nsim,nout),gryy(nsim,nout),grxy(nsim,nout), &
    crxy(nsim,nout),phrxy(nsim,nout), stat = ialloc)
  if (ialloc .ne. 0) call allocerr("a")
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
    call spectr(ini, tx(1:npx), redx(1:npx), ty(1:npy), redy(1:npy),    &
      ofac, n50, iwin, freq, grxx(i,:), gryy(i,:),         &
      grxy(i,:), crxy(i,:), phrxy(i,:))
    !
    !    scale and sum red-noise spectra
    !    -------------------------------
    !
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
  !
  ! scaling factors for red noise from chi^2 distribution
  ! -----------------------------------------------------
  call getdof(iwin, n50, dof, neff)
  !  
  fac90 = getchi2(dof, 0.10_dp) / dof
  fac95 = getchi2(dof, 0.05_dp) / dof
  fac99 = getchi2(dof, 0.01_dp) / dof
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
  if (alpha==0.1)then
    csig_mc = sum(ci90(4,:))/nout
  else if (alpha==0.05_dp)then
    csig_mc = sum(ci95(4,:))/nout
  else if (alpha==0.01_dp)then
    csig_mc = sum(ci99(4,:))/nout
  else 
    csig_mc = 0
  end if
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
  !
  idxph_low = int((alpha/2)*nsim)
  idxph_up = int ((1-alpha/2)*nsim)
  ! 
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
  ! 
  ! For autospectrum yy
  alphacrity = 1 / real(nsegy,dp)
  faccrity = getchi2(dof, alphacrity) / dof
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
  dB6 = winbw(iwin, freq(2)-freq(1), ofac)  ! 6-dB Bandwidth
  false_alarm = (1-alphacritx) * 100 ! Critical false-alarm level
  !                                      ! (Thomson 1990)
  data_x(:, 1) = freq
  data_x(:, 2) = gxx
  data_x(:, 3) = gxxc
  data_x(:, 4) = gredthx
  data_x(:, 5) = grxxavg
  data_x(:, 6) = corrx
  data_x(:, 7) = gredthx*fac90
  data_x(:, 8) = gredthx*fac95
  data_x(:, 9) = gredthx*fac99
  data_x(:, 10) = ci90
  data_x(:, 11) = ci95
  data_x(:, 12) = ci99

  data_y(:, 1) = freq
  data_y(:, 2) = gyy
  data_y(:, 3) = gyyc
  data_y(:, 4) = gredthy
  data_y(:, 5) = gryyavg
  data_y(:, 6) = corry
  data_y(:, 7) = gredthy*fac90
  data_y(:, 8) = gredthy*fac95
  data_y(:, 9) = gredthy*fac99
  data_y(:, 10) = ci90
  data_y(:, 11) = ci95
  data_y(:, 12) = ci99

  data_xy(:, 1) = freq
  data_xy(:, 2) = gxy
  
  data_cxy(:, 1) = freq
  data_cxy(:, 2) = cxy
  data_cxy(:, 3) = csig
  data_cxy(:, 4) = csig_mc
  data_cxy(:, 5) = ci90
  data_cxy(:, 6) = ci95
  data_cxy(:, 7) = ci99
  
  data_phxy(:, 1) = freq
  data_phxy(:, 2) = phxy
  data_phxy(:, 3) = ephi(1, :)
  data_phxy(:, 4) = ephi(2, :)
  data_phxy(:, 5) = ephi_mc(1, :)
  data_phxy(:, 6) = ephi_mc(2, :)

  ! clean up
  ! --------
  deallocate(redx, redy, stat = ialloc)
  deallocate (freq, gxx, gyy, gxy, cxy, phxy)
  deallocate(grxx, gryy, grxy, crxy, phrxy)
  deallocate(grxxsum, gryysum, grxysum,      &
    grxxavg, gryyavg, grxyavg)
  deallocate (gredthx, gredthy, corrx , corry,    &
    gxxc, gyyc)
  deallocate (ci90, ci95, ci99)
  deallocate (gbxx,gbyy,gbxy,cbxy,phbxy,ephi_b,se_phbxy, ephi_mc, &
    se_mc_phxy)
end program redfitx
