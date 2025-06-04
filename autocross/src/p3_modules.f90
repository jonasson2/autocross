! This file contains modules used only by Pearsont3.
!
! These were in pearsont3.f90 in version 1. Also there was a module
! "own_interfaces" with several free-standing subroutines at the end of
! pearsont3.f90. These are now in the module pearsont3_module, except for minls,
! lstau, rhoest and tauest_x that have been moved to common_modules.
!
! Moreover, the module data1 has been split into data1 and data2.

module data1
  ! Data file (original data): n1 rows of (t1, x1, y1).
  use precision
  implicit none
  real(dp), allocatable, dimension (:) :: t1          ! time (original)
  real(dp), allocatable, dimension (:) :: x1          ! data
  real(dp), allocatable, dimension (:) :: y1          ! data
end module data1

!---------------------------------------------------------------------------
module data2
  ! Data file (original data): n1 rows of (t1, x1, y1).
  use precision
  implicit none
  real(dp), allocatable, dimension (:) :: t2          ! time (interval)
  real(dp), allocatable, dimension (:) :: x2          ! data
  real(dp), allocatable, dimension (:) :: y2          ! data
  real(dp), allocatable, dimension (:) :: x3          ! data (detrended)
  real(dp), allocatable, dimension (:) :: y3          ! data (detrended)
  real(dp), allocatable, dimension (:,:) :: x3_resample1  ! resampled x3 - 
                                                          ! formed in the 1st !
                                                          ! bootstrap loop
  real(dp), allocatable, dimension (:,:) :: y3_resample1  ! resampled y3 -
                                                          ! formed in the 1st
                                                          ! bootstrap loop
  real(dp), allocatable, dimension (:) :: x3_resample2  ! resampled x3 - formed
                                                        ! in the 2nd bootstrap
                                                        ! loop
  real(dp), allocatable, dimension (:) :: y3_resample2  ! resampled y3 - formed
                                                        ! in the 2nd bootstrap
                                                        ! loop
end module data2
!
!=============================================================================
!
module parameters
  use precision
  implicit none
  integer, parameter :: nmin=10      ! Minimum number of points
  integer, parameter :: b1=500      ! Number of bootstrap simulations 1
  integer, parameter :: b2=500      ! Number of bootstrap simulations 2	 
  real(dp) :: alpha=0.025_dp         ! Confidence level (1 - 2 * alpha)
  integer, parameter :: imax=10000   ! Big integer
  integer, parameter :: ntry=5       ! user input: maximum number of errors
  integer, parameter :: n_lambda=500 ! Number of grid pts lambda (Calibrated CI)
  !                                  ! n_lambda = 500, gives:
  !                                  ! lambda = 0.001, 0.002,...,0.499, 0.500
  character(len=13), parameter :: &  !
    outputfile='pearsont3.dat'       ! Output file name
end module parameters

!=============================================================================

module inv_tlambda
  use precision
  use parameters, only : n_lambda
  implicit none
  !  Inverse Student's t distribution over lambda range (2nd bootstrap loop)
  real(dp), dimension(n_lambda) :: & ! percentage points of Student's t
    t_inv_lambda=-999.0_dp           ! distribution with v-degrees of freedom
end module inv_tlambda
!
!=============================================================================

module resample_data
  use precision
  use parameters, only: b1,b2,n_lambda
  implicit none
  !       Resampling data.
  integer:: i_lambda
  real(dp), dimension(b1) :: &
    r_resample1=-999.0_dp       ! r_xy for resamples formed in loop1
  real(dp), dimension(b2) :: &
    r_resample2=-999.0_dp       ! r_xy for resamples formed in loop2
  real(dp):: se_r_resample1=-999.0_dp  ! standard error for r_resample1
  real(dp),dimension(b1) :: se_r_resample2=-999.0_dp ! standard error for
  ! r_resample2 (b1 times)
  real(dp), dimension(:,:) ,allocatable :: &
    r_low_resample1, r_upp_resample1
  ! b1 times lower/upper bound (r) for the resamples1 over grid of lambda
  ! Made allocateable to save compilation time and o-file size
  ! (will be allocated size (b1, n_lambda).
  real(dp),dimension(n_lambda),parameter ::         &
    lambda=(/ (0.5_dp*i_lambda/n_lambda, &
    i_lambda=1,n_lambda) /)  ! lambda over grid of values from 0-0.5. how
  !                          ! tight the values are, is decided with n_lambda
end module resample_data
!
!=============================================================================
!
module result1
  use precision
  implicit none
  real (dp) :: r=-999.0_dp     ! r_XY (Pearson's) 
  real (dp) :: r_low=-999.0_dp ! CI lower bound 
  real (dp) :: r_upp=-999.0_dp ! CI upper bound   
  real (dp) :: taux3=-999.0_dp ! persistence time (x3)
  real (dp) :: tauy3=-999.0_dp ! persistence time (y3)
  real (dp) :: rhox3=-999.0_dp ! equivalent autocorrel. coeff. for tau (x3)
  real (dp) :: rhoy3=-999.0_dp ! equivalent autocorrel. coeff. for tau (y3)
  !
end module result1
!
!=============================================================================
!
module setting
  use precision
  implicit none
  character (len=79) :: datafile
  character (len=1) :: dtrtype='m' ! detrending type put to 'm', mean    
  integer :: n1=-999               ! data size (original)
  integer :: n2=-999               ! data size (time interval)      
  integer :: l_mbb = -999.0_dp     ! block length for MBB-moving block bootstrap
end module setting

!=============================================================================
!
!=============================================================================
!
module pearsont3_module
contains
  subroutine allocate_resample_data
    use setting
    use resample_data
    implicit none
    integer error
    allocate(r_low_resample1(b1, n_lambda), r_upp_resample1(b1, n_lambda), &
      stat=error)
    if (error /= 0) stop 'Allocation failed'
  end subroutine allocate_resample_data

  subroutine deallocate_resample_data
    use resample_data
    implicit none
    deallocate(r_low_resample1, r_upp_resample1)
  end subroutine deallocate_resample_data
  
  subroutine allocate0
    use precision
    use setting, only: n1
    use data1, only: t1,x1,y1
    implicit none
    !       Allocates t1, x1, y1.
    allocate(t1(n1))
    allocate(x1(n1))
    allocate(y1(n1)) 
    ! if (error /= 0) then
    !   print *,'Subroutine allocate0:'
    !   print *,'  space requested not possible - PearsonT terminates.'
    !   stop
    ! end if
  end subroutine allocate0
  !
  !=============================================================================
  !
  subroutine allocate1
    use setting, only: n2
    use data2, only: t2,x2,y2,x3,y3,x3_resample1,y3_resample1,      &
      x3_resample2,y3_resample2
    use parameters, only: b1
    implicit none
    !       Allocates t2, x2, y2, x3, y3, x3_resample1,y3_resample1,
    !       x3_resample2, y3_resample2
    integer :: error=0
    allocate(t2(n2),x2(n2),y2(n2),x3(n2),y3(n2),                    &
      x3_resample1(b1,n2),y3_resample1(b1,n2),               &
      x3_resample2(n2),y3_resample2(n2), stat=error)
    if (error /= 0) then
      print *,'Subroutine allocate1: space requested not possible - PearsonT terminates.'
      stop
    end if
  end subroutine allocate1
  !
  !=============================================================================
  !
  subroutine bootstrap(n,x,y,l,x_resample,y_resample)
    use precision
    use random
    implicit none
    integer, intent(in) :: n          ! number of data points
    real(dp), dimension(:), intent(in)  :: x      ! original x time series
    real(dp), dimension(:), intent(in)  :: y      ! original y time series
    integer, intent(in) :: l          ! block length
    real(dp), dimension(:), intent(out) :: x_resample   ! bootstrap resampled x
    ! time series
    real(dp), dimension(:), intent(out) :: y_resample   ! bootstrap resampled y
    ! time series
    !
    ! Pairwise Moving Block bootstrap resampling Algorithm 7.2 and 3.1 in
    ! (Mudelsee, 2010) Block length = l = l_mbb calculated in subroutine
    ! chsett4. Blocks are randomly selected and concatenated until n data are
    ! resampled (resize last block by removing last observation in it). By using
    ! the same random bootstrap index for x_resample and y_resample, (x(i),y(i))
    ! pairs are resampled.
    !
    integer, dimension(n) :: indxx
    integer :: i=0
    integer :: j=0
    integer :: k=0
    integer :: n_block_start
    n_block_start=n-l+1
    if (l == 1) then         ! block length less or equal to average spacing,
      ! ordinary bootstrap is used
      do i=1,n
        indxx(i)=int(n*uniform_random())+1
      end do
    else
      k=1            ! set counter
      outer:     do i=1,n
        indxx(k)=int(n_block_start*uniform_random())+1  ! draw random block start
        k=k+1
        if (k > n) exit outer
        do j=1,l-1        ! fill up the block of length l
          indxx(k)=indxx(k-1)+1
          k=k+1
          if (k > n) exit outer
        end do          ! the loop keeps on going until k>n or k==n
      end do outer
    end if
    x_resample(:)=x(indxx(:))
    y_resample(:)=y(indxx(:))
    !
  end subroutine bootstrap

  !============================================================================

  subroutine calc_t_inv_lambda
    use precision
    use inv_tlambda
    use parameters, only: n_lambda
    use setting, only: n2
    use resample_data, only: lambda
    !       Calculates inverse of t(alpha, dof).
    !       dof = n -1 (mean estimation) (von Storch and Zwiers 1999: p. 92)
    !
    !   calculates for lambda grid; dof = 2n-5 (x, y; E_x, E_y, S_x,S_y, r_xy
    !   estimation).
    !
    integer:: i
    !
    do i=1,n_lambda
      t_inv_lambda(i)=t_inv(lambda(i),2*n2-5)
    end do
    !    
  end subroutine calc_t_inv_lambda
  !
  !=============================================================================
  !
  subroutine chsett1
    use parameters, only: ntry
    use setting, only: datafile
    implicit none
    !       Changes setting: datafile.
    integer :: i,open_error
    do i=1,ntry
      print *
      print '(a)','data (t, x, y)              [path + filename]'
      !read (5,'(a)') datafile
      datafile = "test_data.txt"
      open (unit=1, file=datafile, status='old',                   &
        form='formatted', action='read', iostat=open_error)
      if (open_error /= 0 ) then
        print *,'Error during data file opening - try again'
      else
        exit
      end if
    end do
    if (i > ntry) then
      print *,'OK - PearsonT terminates.'
      stop
    end if
  end subroutine chsett1
  !
  !=============================================================================
  !
  subroutine chsett2
    use precision
    use parameters, only: nmin
    use setting, only: datafile,n1
    implicit none
    !       Determines n, see REDFIT35.F90 subroutine setdim
    !       (Schulz and Mudelsee 2002).
    character (len = 1) :: flag
    integer :: i,iocheck,open_error
    real(dp) :: tdum, xdum, ydum  ! dummy
    open (unit=1, file=datafile, status='old',                      &
      form='formatted', action='read', iostat=open_error)
    if (open_error /= 0 ) then
      print *,'Error during data file opening - PearsonT terminates.'
      stop
    end if
    !
    ! 1.    Skip header
    !       ==========
    !
    do while (.true.)
      read (1, '(a1)') flag
      if (flag .ne. '#') then
        backspace (1)
        exit
      end if
    end do
    !
    ! 2.    Count data
    !       =========
    i = 1
    do while (.true.)
      read (1, *, iostat = iocheck) tdum, xdum, ydum
      if (iocheck .ne. 0) exit
      i = i + 1
    end do
    close (unit=1, status='keep')
    n1 = i - 1
    if (n1 < nmin) then
      print '(a,i2)','data size (original) : ',n1
      print '(a,i2)',' minimum required is : ',nmin
      print *
      print *,'PearsonT terminates.'
      stop
    end if
  end subroutine chsett2
  !
  !=============================================================================
  !
  !subroutine chsett3
  !        use parameters, only: ntry
  !        use setting, only: dtrtype
  !        implicit none
  !       Changes setting: detrending type.
  !        integer :: i
  !        do i=1,ntry
  !           print *
  !           print '(a)','detrending type:                   linear [l]'
  !           print '(a)','                                     mean [m]'
  !           read (5,'(a)') dtrtype
  !           if (dtrtype /= 'l' .and.dtrtype /= 'm') then
  !              print *,'Not possible - try again.'
  !           else
  !              exit
  !           end if
  !        end do
  !        if (i > ntry) then
  !           print *,'OK - PearsonT terminates.'
  !           stop
  !        end if
  !end subroutine chsett3
  !
  !=============================================================================
  !
  subroutine chsett4
    use precision
    use result1, only: rhox3,rhoy3   
    use setting, only: n2, l_mbb
    implicit none
    !   
    !  Changes setting: Calculates block length (formula from Carlstein et al,
    !  1986. Sherman et al., 1998 adapted it to MBB) Equations 7.31 and 7.32 in
    !  Mudelsee, 2010
    !
    real(dp) :: axy      
    real(dp), parameter :: onethird=0.3333333333333333333333333333333333333333_dp
    real(dp), parameter :: twothird=0.6666666666666666666666666666666666666666_dp
    real(dp), parameter :: half=0.5_dp
    real(dp), parameter :: one=1.0_dp
    real(dp), parameter :: two=2.0_dp
    !
    axy = sqrt(rhox3*rhoy3) !Eq 7.31 in (Mudelsee, 2010)
    !  
    l_mbb = nint(((6.0_dp**half*axy)/(one-axy**2))**twothird*n2**onethird)        
    !                          ! Eq 7.32 in (Mudelsee,2010)
    if (l_mbb < 1) then
      print*
      print'(a)','Note: block length less than average spacing -'
      print'(a)','      PearsonT uses ordinary bootstrap.'
      l_mbb = 1
    end if
    if (l_mbb == 1) then
      print*
      print'(a)','Note: block length equal to average spacing - '
      print'(a)','      PearsonT uses ordinary bootstrap.'
    end if
    if (l_mbb> n2) then
      l_mbb = n2-1
      print'(a)','Note: block length longer than number of data points n '
      print'(a)','      PearsonT uses block lenght =  n -1.'
    end if
    !
  end subroutine chsett4
  !
  !=============================================================================
  !
  subroutine confidence
    use precision
    use meanvar_module
    use data2, only: x3,x3_resample1,x3_resample2,                  &
      y3,y3_resample1,y3_resample2
    use inv_tlambda
    use parameters, only: alpha, b1, b2, n_lambda
    use resample_data
    use result1, only: r,r_low,r_upp
    use setting, only: n2,l_mbb
    use time
    implicit none
    !
    !       Calibrated Student's t confidence interval for r_xy 
    !
    integer :: j=0      ! for j=1,b1
    integer :: k=0      ! for k=1,b2
    integer :: l=0      ! for l=1,n_lambda
    !  
    real(dp) :: ave_se_1=-999.0_dp
    real(dp) :: var_se_1=-999.0_dp
    real(dp) :: ave_se_2=-999.0_dp
    real(dp) :: var_se_2=-999.0_dp
    !  
    real(dp),dimension(n_lambda)::plambda=-999.0_dp
    real(dp):: alphatest=-999.0_dp        ! test alpha (1 - 2 * alpha)
    integer :: l_loc=0 ! gives the location of lambda value, within the lambda,
    !                  ! which gives the coverage closest to 1-2*alpha  array
    character(len=1),parameter::flag_fisher='y'    ! 'y' = Fisher's transform
    !                ! 'n' = without Fisher's transform

    real(dp) :: bootstrap_time, pearsn_time
    
    !
    !   1. First bootstrap loop (Student's t confidence interval)
    !   ========================================================
    ! 
    !  1.1 Bootstrap
    !  ============
    !  Form bootstrap resamples x3_resample1 and y3_resample1 and store it
    !  Estimate Pearson's correlation coefficient, r_resample1 and store it
    !
    do j=1,b1
      call bootstrap(n2,x3,y3,l_mbb,x3_resample1(j,1:n2),y3_resample1(j,1:n2))
      call pearsn(x3_resample1(j,1:n2),y3_resample1(j,1:n2),n2,r_resample1(j))
    end do
    !
    !  1.2 Fisher's z-transformation
    !  ============================
    !  z-transform the r_resamples z=invtanh(r_resample1)
    !
    if (flag_fisher == 'y')then
      do j=1,b1
        r_resample1(j)=invtanh(r_resample1(j))
      end do

    end if
    !
    !  1.3 Bootstrap standard error
    !  ===========================
    !  Estimate the bootstrap se from all the replications
    !  
    call meanvar(r_resample1(1:b1),ave_se_1,var_se_1)
    se_r_resample1 = sqrt(var_se_1)  
    !
    ! 2. Second bootstrap loop (Forms CI for each resample1)
    ! ========================================================
    !
    !  2.1 Big bootstrap loop
    !  =====================
    !  Goes b1=2000 times through the 2nd bootstrap loop

    !  
    boot:  do j=1,b1
      !
      !  2.2 2nd Bootstrap
      !  ================
      !  Forms bootstrap resamples x3_resample2 and y3_resample2 from the
      !  resamples from bootstrap loop 1. Forms b1*b2 resamples2 (no need to
      !  store them between loops). The block length is overtaken from the
      !  original samples x3 and y3. Estimates Pearson's correlation coefficient
      !  r_resample2 and stores it b2 times (just within this loop, not within
      !  the big loop).
      !
      call tic()
      bootstrap_time = 0
      pearsn_time = 0
      do k=1,b2
        call bootstrap(n2,x3_resample1(j,1:n2),y3_resample1(j,1:n2),l_mbb, &
          x3_resample2(1:n2),y3_resample2(1:n2))
        bootstrap_time = bootstrap_time + toc()
        call pearsn(x3_resample2(1:n2),y3_resample2(1:n2),n2,r_resample2(k))
        pearsn_time = pearsn_time + toc()
      end do
      !
      !  2.3 Fisher's z-transform
      !  =======================
      !
      if (flag_fisher == 'y')then
        do k=1,b2
          r_resample2(k) = invtanh(r_resample2(k))
        end do
      end if
      !
      !  2.4 Confidence intervals for resamples 1
      !  =======================================
      !
      !  2.4.1 Bootstrap standard error
      !  =============================
      !  Estimates the bootstrap standard error from the replications,
      !  r_resample2. Results with b1 times se, one for each big loop or one for
      !  each of the 2000 CI for resample1
      !  
      call meanvar(r_resample2(1:b2),ave_se_2,var_se_2)
      se_r_resample2(j) = sqrt(var_se_2)    
      !  
      !
      !  2.4.2 Student's t CI for each resample1 (b1 times r_upp and r_low)
      !  =================================================================
      !  r_upp and r_low estimated for each bootstrap resample1 (b1 times)
      !  over a grid of confidence levels (lambda)
      !  
      if (flag_fisher == 'y')then
        do l=1,n_lambda
          r_low_resample1(j,l) = tanh(r_resample1(j) + &
            t_inv_lambda(l)*se_r_resample2(j))
          r_upp_resample1(j,l) = tanh(r_resample1(j) - &
            t_inv_lambda(l)*se_r_resample2(j))
        end do
      else if (flag_fisher == 'n')then
        do l=1,n_lambda
          r_low_resample1(j,l) = r_resample1(j) + &
            t_inv_lambda(l)*se_r_resample2(j)
          r_upp_resample1(j,l) = r_resample1(j) - &
            t_inv_lambda(l)*se_r_resample2(j)
        end do
      end if
      !    
      !
      !  2.5 End of big Bootstrap loop
      !  ============================
    end do boot
    !  
    !
    ! 3. Determination of two-sides plambda
    ! ====================================
    !   Equation 3.47 in Mudelsee(2010)
    !
    do l = 1,n_lambda
      plambda(l)= 1.0_dp *count(r .ge. r_low_resample1(1:b1,l) .and.  &
        r .le. r_upp_resample1(1:b1,l))/b1
    end do
    !  
    ! 4. Calibration - calibrated Student's t confidence interval
    ! ==========================================================
    !  
    !  4.1 Find new confidence points
    !  =============================
    !  Find the lambda value which gives coverage closest to the nominal value
    !  (1-2*alpha)
    !    
    alphatest=1.0_dp-(2.0_dp*alpha)
    !
    l_loc = maxval(minloc(abs(alphatest-plambda(:))))
    !
    !  4.2 Calculate calibrated Student's t CI
    !  ======================================
    !
    if (flag_fisher == 'y')then    
      r_low = tanh(invtanh(r) + t_inv_lambda(l_loc) * se_r_resample1)
      r_upp = tanh(invtanh(r) - t_inv_lambda(l_loc) * se_r_resample1)  
    else if (flag_fisher == 'n')then 
      r_low = r + t_inv_lambda(l_loc) * se_r_resample1
      r_upp = r - t_inv_lambda(l_loc) * se_r_resample1     
    end if
    !
  end subroutine confidence
  !
  !=============================================================================
  !
  subroutine deallocate0
    use data1, only: t1, x1, y1
    implicit none
    integer :: error=0
    deallocate(t1, x1, y1, stat=error)
    if (error /= 0) then
      print *,'Subroutine deallocate0:'
      print *,'  Unexpected deallocation error - PearsonT terminates.'
      stop
    end if
  end subroutine deallocate0

  subroutine deallocate1
    use data2, only: t2,x2,y2,x3,y3,x3_resample1,y3_resample1,  &
      x3_resample2,y3_resample2
    implicit none
    !       Deallocates t2, x2, y2, x3, y3, x3_resample1, y3_resample1,
    !       x3_resample2, y3_remsample2
    integer :: error=0
    !  
    deallocate(t2, x2, y2, x3, y3, x3_resample1, y3_resample1,      &
      x3_resample2,y3_resample2, stat=error)
    !
    if (error /= 0) then
      print *,'Subroutine deallocate1:'
      print *,'  Unexpected deallocation error - PearsonT terminates.'
      stop
    end if
  end subroutine deallocate1
  !
  !=============================================================================
  !
  subroutine decis1(c)
    use parameters, only: ntry
    implicit none
    character(len=1), intent(out) :: c
    !       Part 1 : Time interval extraction.
    integer :: i
    print *
    print *
    do i=1,ntry
      print *
      print '(a)','PearsonT3'
      print *
      print '(a)','New time interval                        [n]'
      print '(a)','Original time interval                   [o]'
      print '(a)','Output and exit                          [x]'
      !read (5,'(a1)') c
      c = 'x'
      if (c == 'n' .or. c == 'o' .or. c == 'x') exit
    end do
    if (i > ntry) then
      print *,'OK - PearsonT3 terminates.'
      stop
    end if
  end subroutine decis1
  !
  !=============================================================================
  !
  subroutine fit(x,y,ndata,a,b)
    use precision
    implicit none
    integer, intent(in) :: ndata
    real(dp), dimension(ndata), intent(in) :: x,y
    real(dp), intent(out) :: a,b
    real(dp) :: ss,sx,sxoss,sy,st2
    real(dp), dimension(ndata) :: t
    ss=real(ndata,dp)
    sx=sum(x)
    sy=sum(y)
    sxoss=sx/ss
    t(:)=x(:)-sxoss
    b=dot_product(t,y)
    st2=dot_product(t,t)
    b=b/st2
    a=(sy-sx*b)/ss
  end subroutine fit
  !
  !=============================================================================
  !
  subroutine info0
    implicit none
    !       Welcomes.
    integer :: i
    do i=1,10
      print '(a)',' '
    end do
    print '(a)','=============================================================='
    print '(a)','                                                              '
    print '(a,a,a)', 'PearsonT3 (Version 1.0',',',' May 2013)                  '
    print '(a)','                                                              '
    print '(a)','=============================================================='
    print '(a)','                                                              '
    print '(a)','  Reference:                                                  '
    print '(a)','    Olafsdottir and Mudelsee: More accurate, calibrated       '
    print '(a)','    bootstrap confidence intervals for correlating two        '
    print '(a)','    serially dependent time series (submitted to Mathematical '
    print '(a)','    Geosciences) '
    print '(a)','                                                              '
    print '(a)','  Download:                                                   '
    print '(a)','    http://www.climate-risk-analysis.com                      '
    print '(a)','                                                              '
    print '(a)','=============================================================='
    print '(a)','                                                              '
  end subroutine info0
  !
  !=============================================================================
  !
  subroutine info1(c,filec)
    use setting, only: datafile,n1,n2
    use data1, only: t1
    use data2, only: t2
    implicit none
    character(len=1), intent(in) :: c
    character(len=13), intent(in), optional :: filec
    !       Informs about: filename, datatype, number of points:
    !       observed and extracted time interval.
    integer :: open_error
    if (c == 'p') then
      print '(a,a)',                                                       &
        'Data file name:              ',datafile
      print '(a,f11.3,a,f11.3,a,i8,a)',                                  &
        'Time interval - original:    [ ',t1(1),'; ',t1(n1),' ]     - ',n1, &
      ' points'
      if (n1 /= n2)                                                      &
        print '(a,f11.3,a,f11.3,a,i8,a)',                               &
        'Time interval - extracted:   [ ',t2(1),'; ',t2(n2),' ]     - ',n2, &
        ' points'
    else if (c == 'w') then
      open (unit=1, file=filec, status='replace',                         &
        form='formatted', action='write', iostat=open_error)
      if (open_error /= 0 ) then
        print *,'Error during data file opening - terminate'
        stop
      end if
      write (unit=1, fmt='(a)')                                          &
        'PearsonT3 (www.climate-risk-analysis.com):  Output'
      write (unit=1, fmt='(a)')                                          &
        ' '
      write (unit=1, fmt='(a,a)')                                        &
        'Data file name:              ',datafile
      write (unit=1, fmt='(a,f11.3,a,f11.3,a,i8,a)')                 &
        'Time interval - original:    [ ',t1(1),'; ',t1(n1),' ]     - ',n1,&
        ' points'
      if (n1 /= n2)                                                      &
        write (unit=1, fmt='(a,f11.3,a,f11.3,a,i8,a)')              &
        'Time interval - extracted:   [ ',t2(1),'; ',t2(n2),' ]     - ',n2,&
        ' points'
      close (unit=1, status='keep')
    end if
  end subroutine info1
  !
  !============================================================================
  !
  subroutine info2(c,filec)
    use result1, only: taux3,tauy3
    implicit none
    character(len=1), intent(in) :: c
    character(len=13), intent(in), optional :: filec
    !       Informs about: taux3, tauy3.
    integer :: open_error
    if (c == 'p') then
      print *
      print '(a,23x,f11.3)',                                           &
        'Persistence time (x):',taux3
      print '(a,23x,f11.3)',                                           &
        'Persistence time (y):',tauy3
    else if (c == 'w') then
      open (unit=1, file=filec, status='old', position='append',        &
        form='formatted', action='write', iostat=open_error)
      if (open_error /= 0 ) then
        print *,'Error during data file opening - PearsonT terminates.'
        stop
      end if
      write (unit=1, fmt='(a)')                                        &
        ' '
      write (unit=1, fmt='(a,23x,f11.3)')                              &
        'Persistence time (x):',taux3
      write (unit=1, fmt='(a,23x,f11.3)')                              &
        'Persistence time (y):',tauy3
      close (unit=1, status='keep')
    end if
  end subroutine info2
  !
  !============================================================================
  !
  subroutine info3(c,filec)
    use setting, only: dtrtype
    implicit none
    character(len=1), intent(in) :: c
    character(len=13), intent(in), optional :: filec
    !       Informs about: dtrtype.
    integer :: open_error
    if (c == 'p') then
      print *
      if (dtrtype .eq. 'l') then
        print '(a,29x,a)',                                            &
          'Detrending type:     ','linear'
      else if (dtrtype .eq. 'm') then
        print '(a,29x,a)',                                            &
          'Detrending type:     ','mean'
      end if
    else if (c == 'w') then
      open (unit=1, file=filec, status='old', position='append',        &
        form='formatted', action='write', iostat=open_error)
      if (open_error /= 0 ) then
        print *,'Error during data file opening - PearsonT terminates.'
        stop
      end if
      write (unit=1, fmt='(a)')                                        &
        ' '
      if (dtrtype .eq. 'l') then
        write (unit=1, fmt='(a,29x,a)')                               &
          'Detrending type:     ','linear'
      else if (dtrtype .eq. 'm') then
        write (unit=1, fmt='(a,29x,a)')                               &
          'Detrending type:     ','mean'
      end if
      close (unit=1, status='keep')
    end if
  end subroutine info3
  !
  !============================================================================
  !
  subroutine info4(c,filec)
    use precision
    use data2, only: t2,x2,y2,x3,y3
    use parameters, only: alpha
    use result1, only: r,r_low,r_upp
    use setting, only: n2
    implicit none
    character(len=1), intent(in) :: c
    character(len=13), intent(in), optional :: filec
    !       Informs about: r, r_low, r_upp.
    integer :: open_error
    integer :: i
    if (c == 'p') then
      print *
      print '(a,a,i2,a,11x,4(f6.3,a))',                                &
        'Pearson''s r',', ',                      &
        nint(100*(1-2.0_dp*alpha)),' % confidence interval:',       &
        r, '  [ ',r_low,'; ',r_upp,' ]'
    else if (c == 'w') then
      open (unit=1, file=filec, status='old', position='append',        &
        form='formatted', action='write', iostat=open_error)
      if (open_error /= 0 ) then
        print *,'Error during data file opening - PearsonT terminates.'
        stop
      end if
      write (unit=1, fmt='(a)') ' '
      write (unit=1, fmt='(a,a,i2,a,11x,4(f6.3,a))')                     &
        'Pearson''s r',', ',                                       &
        nint(100*(1-2.0_dp*alpha)),' % confidence interval:',        &
        r, '  [ ',r_low,'; ',r_upp,' ]'
      write (unit=1, fmt='(a)') ' '
      write (unit=1,fmt='(2x,7(a))')                                   &
        ' time          ',                                         &
        ' x             ',                                         &
        ' y             ',                                           &
        ' trend(x)      ',                                         &
        ' trend(y)      ',                                         &
        ' x-detrended   ',                                         &
        ' y-detrended   '
      do i=1,n2
        write (unit=1,fmt='(7(1x,es14.6))') t2(i),x2(i),y2(i),        &
          x2(i)-x3(i),              &
          y2(i)-y3(i),              &
          x3(i),y3(i)
      end do
      close (unit=1, status='keep')
    end if
  end subroutine info4
  !
  !=============================================================================
  !
  subroutine init0
    use precision
    use data1, only: t1, x1, y1
    implicit none
    !          Initializes t1, x1, y1.
    t1=-999.0_dp
    x1=-999.0_dp
    y1=-999.0_dp
  end subroutine init0
  !
  !=============================================================================
  !
  subroutine init1a
    use setting, only: n1,n2
    implicit none
    !         Initializes: n2=n1.
    n2=n1
  end subroutine init1a
  !
  !=============================================================================
  !
  subroutine init1b
    use precision
    use data1, only: t1,x1,y1
    use data2, only: t2,x2,x3,y2,y3,x3_resample1,y3_resample1, x3_resample2,y3_resample2
    implicit none
    !       Initializes t2, x2, y2, x3, y3, x3_resample1, y3_resample1,
    !       x3_resample2, y3_resample2
    t2=t1
    x2=x1
    y2=y1
    x3=x2
    y3=y2
    x3_resample1=-999.0_dp
    y3_resample1=-999.0_dp
    x3_resample2=-999.0_dp
    y3_resample2=-999.0_dp
    !
  end subroutine init1b
  !
  !=============================================================================
  !
  ! subroutine interval(c)
  !   use precision
  !   use data1, only: t1,x1,y1
  !       	use data2, only: t2,x2,y2,x3,y3,x3_resample1,y3_resample1,  &
  !     x3_resample2,y3_resample2
  !   use parameters, only: nmin,ntry
  !   use setting, only: n1,n2
  !   implicit none
  !   character(len=1), intent(in) :: c
  !   !       Defines new: n2, t2, x2, y2.
  !   integer :: i1,i2,j,k1
  !   real(dp) :: tl,tr
  !   if (c == 'n') then
  !     do i1=1,ntry
  !       do i2=1,ntry
  !         print '(a,i6,a)','New interval must contain at least',nmin,' points'
  !         print '(3(a))',  '[left boundary',',',' right boundary]: '
  !         read (5,*) tl, tr
  !         tl=max(tl,t1(1))
  !         tr=min(tr,t1(n1))
  !         if (tl < t1(n1) .and. tr > t1(1) .and. tl < tr) exit
  !       end do
  !       if (i2 > ntry) then
  !         print *,'OK - PearsonT terminates.'
  !         stop
  !       end if
  !       if (tl == t1(1)) j=1
  !       if (tr == t1(n1)) k1=n1
  !       do i2=2,n1
  !         if (tl >  t1(i2-1) .and. tl <= t1(i2)) j=i2
  !         if (tr >= t1(i2-1) .and. tr <  t1(i2)) k1=i2-1
  !       end do
  !       n2=k1-j+1
  !       if (n2 >= nmin) exit
  !       print '(a,i6)','Time interval contains too few points: ',n2
  !     end do
  !     if (i1 > ntry) then
  !       print *,'OK - PearsonT terminates.'
  !       stop
  !     end if
  !     call deallocate1  ! t2, x2, y2, x3, y3, x3_resample1, y3_resample1,
  !     ! x3_resample2, y3_resample2
  !     call allocate1     ! t2, x2, y2, x3, y3, x3_resample1, y3_resample1,
  !     ! x3_resample2, y3_resample2
  !     do i2=j,k1
  !       t2(i2-j+1)=t1(i2)
  !       x2(i2-j+1)=x1(i2)
  !       y2(i2-j+1)=y1(i2)
  !     end do
  !     x3=x2
  !     y3=y2
  !     x3_resample1=-999.0_dp
  !     y3_resample1=-999.0_dp
  !     x3_resample2=-999.0_dp
  !     y3_resample2=-999.0_dp
  !   else if (c == 'o') then
  !     call init1a        ! n2=n1
  !     call deallocate1  ! t2, x2, y2, x3, y3, x3_resample1, y3_resample1,
  !     ! x3_resample2, y3_resample2
  !     call allocate1     ! t2, x2, y2, x3, y3, x3_resample1, y3_resample1,
  !     ! x3_resample2, y3_resample2
  !     call init1b        ! t2, x2, y2, x3, y3, x3_resample1, y3_resample1,
  !     ! x3_resample2, y3_resample2
  !   end if
  ! end subroutine interval
  !
  !=============================================================================
  !
  function invtanh(x)
    use precision
    implicit none
    real(dp) :: invtanh
    !       Inverse hyperbolic tangens.
    !
    real(dp),  intent(in) :: x
    real(dp), parameter :: tinyx=1.0e-13
    real(dp), parameter :: half=0.5_dp
    real(dp), parameter :: one=1.0_dp
    real(dp), parameter :: big=16.0_dp
    if (x.ge.one-tinyx) then
      invtanh=big
    else if (x.le.-one+tinyx)then
      invtanh=-big  
    else
      invtanh=half*log((one+x)/(one-x))
    end if
  end function invtanh
  !
  !=============================================================================
  !
  subroutine output
    use parameters, only: outputfile
    implicit none
    !       Writes output to file.
    print *
    print '(a,a)','Output written to file ',outputfile
    call info1('w',outputfile)
    call info2('w',outputfile)
    call info3('w',outputfile)
    call info4('w',outputfile)
  end subroutine output
  !
  !=============================================================================
  !
  subroutine pearsn(x,y,n,r)
    use precision
    implicit none
    !       Pearson's correlation coefficient (Numerical Recipes, modified).
    integer, intent(in) :: n
    real(dp), dimension(:), intent(in) :: x,y
    real(dp), intent(out) :: r
    real(dp), parameter :: tiny=1.0e-08_dp
    real(dp), dimension(n) :: xt,yt
    real(dp) :: ax,ay,sxx,sxy,syy
    ax=sum(x)/n
    ay=sum(y)/n
    xt(:)=x(:)-ax
    yt(:)=y(:)-ay
    sxx=dot_product(xt,xt)
    syy=dot_product(yt,yt)
    sxy=dot_product(xt,yt)
    r=sxy/(sqrt(sxx*syy)+tiny)
  end subroutine pearsn
  !
  !=============================================================================
  !
  function phi(x)
    use precision
    implicit none
    real(dp), parameter :: sqrt2 = sqrt(2.0_dp)
    real(dp) :: phi
    real(dp), intent(in) :: x
    !       Cumulative normal density distribution,
    !       calculated using the built-in erfc
    !       NR: Fractional error [phi] < 1.2e-7.
    phi=1.0_dp-0.5_dp*erfc(x/sqrt2)
  end function phi
  !
  !=============================================================================
  !
  real(dp) function phi_inv(prob)
    use precision
    implicit none
    !       Inverse cumulative normal density approximation after
    !       Odeh and Evans (1974, Appl. Stat. 23, 96-97).
    !       Error [phi_inv] < 1.5e-8.
    real(dp), intent(in) :: prob
    real(dp), parameter ::  big=1.0e+38_dp
    real(dp), parameter :: plim=1.0e-20_dp
    real(dp), parameter :: p0=-0.322232431088e0_dp
    real(dp), parameter :: p1=-1.0e0_dp
    real(dp), parameter :: p2=-0.342242088547e0_dp
    real(dp), parameter :: p3=-0.0204231210245e0_dp
    real(dp), parameter :: p4=-0.453642210148e-04_dp
    real(dp), parameter :: q0= 0.0993484626060e0_dp
    real(dp), parameter :: q1= 0.588581570495e0_dp
    real(dp), parameter :: q2= 0.531103462366e0_dp
    real(dp), parameter :: q3= 0.103537752850e0_dp
    real(dp), parameter :: q4= 0.38560700634e-02_dp
    real(dp) :: p=-999.0_dp
    real(dp) :: y=-999.0_dp
    real(dp), parameter :: half=0.5_dp
    real(dp), parameter ::  one=1.0_dp
    real(dp), parameter :: zero=0.0_dp
    phi_inv=-big
    p=prob
    !
    !       Ranges for p
    !       ===========
    !
    if (p < plim) then
      phi_inv=-big
    else if (p > one-plim) then
      phi_inv=big
    else if (p == half) then
      phi_inv=zero
    else if (p >= plim .and. p < half) then
      y=sqrt(log(1.0_dp/(p*p)))
      phi_inv=y+((((y*p4+p3)*y+p2)*y+p1)*y+p0)/((((y*q4+q3)*y+q2)*y+q1)*y+q0)
      phi_inv=-1.0_dp*phi_inv
    else if (p > half .and. p <= one-plim) then
      p=one-p
      y=sqrt(log(1.0_dp/(p*p)))
      phi_inv=y+((((y*p4+p3)*y+p2)*y+p1)*y+p0)/((((y*q4+q3)*y+q2)*y+q1)*y+q0)
    end if
  end function phi_inv
  !
  !=============================================================================
  !
  subroutine plot(type)
    use precision
    use data2, only: t2,x2,x3,y2,y3
    use setting, only: n1,n2
    implicit none
    integer, intent(in) :: type             ! Plot type
    !
    !       Plots data.
    !
    !       Plot type                   Description
    !       ======================================
    !       1               x2          vs      t2
    !                       y2          vs      t2
    !                 trend(x2)         vs      t2
    !                 trend(y2)         vs      t2
    !       2               y2          vs      x2
    !
    character(len=5), parameter :: file1='1.tmp'
    character(len=5), parameter :: file2='2.tmp'
    character(len=11) :: title_1
    character(len=26) :: title_2
    character(len=12) :: title_3
    integer :: i,open_error
    !
    ! 1.    Check type
    !       =========
    !
    if (type /=1 .and. type /= 2) then
      print *,'Plot-type value out of range - PearsonT terminates.'
      stop
    end if
    !
    ! 2.    Write plot data
    !       ==============
    !
    open (unit=1, file=file1, status='unknown',                     &
      form='formatted', action='write', iostat=open_error)
    if (open_error /= 0 ) then
      print *,'Error during plot data file opening - PearsonT terminates.'
      stop
    end if
    if (type == 1 .or. type == 2) then
      write (unit=1,fmt='(2x,7(a))')                                 &
        ' time          ',                                         &
        ' x             ',                                         &
        ' y             ',                                         &
        ' trend(x)      ',                                         &
        ' trend(y)      ',                                         &
        ' x-detrended   ',                                         &
        ' y-detrended   '
      do i=1,n2
        write (unit=1,fmt='(7(1x,es14.6))') t2(i),x2(i),y2(i), &
          x2(i)-x3(i),              &
          y2(i)-y3(i),              &
          x3(i),y3(i)
      end do
    end if
    close (unit=1, status='keep')
    !
    ! 3.    Write gnuplot scriptfile
    !       =======================
    !
    !
    ! 3.1   Key and grid
    !       ===========
    !
    open (unit=2, file=file2, status='unknown',              &
      form='formatted', action='write', iostat=open_error)
    !
    !***************************************************************************
    !       IBM-PC, Windows 9x/NT, DOS window             AbSoft 6.2
    !
    !       IBM-PC, Linux                                 Absoft 6.0
    !
    !        Sun workstation, Unix, X11                   EPCF77 2.7
    !       <=> following lines uncommented
    !

    write (unit=2, fmt='(a)') ' set key default'
    if (type == 2) then
      write (unit=2, fmt='(a)') ' set grid'
    end if
    !***************************************************************************
    !
    !       See notes to subroutine plotgn.
    !
    !***************************************************************************
    !       PC, Suse Linux 8.0                            Portland PGF90 3.2-4
    !       <=> following lines uncommented
    !
    !       if (type == 1) then
    !          write (unit=2, fmt='(a)') ' set key'
    !       else if (type == 2) then
    !          write (unit=2, fmt='(a,i8,a,f6.3,a,f6.3,a,f6.3,a)')          &
    !                ' set key right Left samplen 0.0  title  "n  = ',n2,   &
    !                '\\nr  =   ',r, '  [ ',r_low,'; ',r_upp,' ]"'
    !          write (unit=2, fmt='(a)') ' set grid'
    !       end if
    !***************************************************************************
    !
    ! 3.2   Titles
    !       ======
    !
    if (type == 1 .or. type ==2) then
      title_1='Time series'
    end if
    if (n1 == n2) then
      title_2=' - original time interval '
    else
      title_2=' - extracted time interval'
    end if
    if (type == 1) then
      title_3='            '
    else if (type == 2) then
      title_3=' - detrended'
    end if
    if (type == 1 .or. type ==2) then
      write (unit=2, fmt='(3(a))')                                 &
        ' set title "',title_1//title_2//title_3,'"'
    end if
    !
    ! 3.3   Axis labels
    !       ==========
    !
    if (type == 1) then
      write (unit=2, fmt='(a)') ' set xlabel "time  t"'
    else if (type == 2) then
      write (unit=2, fmt='(a)') ' set xlabel "time series value  x"'
    end if
    if (type == 1) then
      write (unit=2, fmt='(a)') ' set ylabel "time series value  x"'
      write (unit=2, fmt='(a)') ' set y2label "time series value  y"'
    else if (type == 2) then
      write (unit=2, fmt='(a)') ' set ylabel "time series value  y"'
    end if
    !
    ! 3.4   y-axes
    !       ======
    !
    if (type == 1) then
      write (unit=2, fmt='(a)') ' set ytics nomirror'
      write (unit=2, fmt='(a)') ' set y2tics nomirror'
    end if
    !
    ! 3.5   xrange, yrange
    !       =============
    !
    if (type == 1) then
      write (unit=2, fmt='(a)') ' set autoscale y'
      write (unit=2, fmt='(a)') ' set autoscale y2'
    else if (type == 2) then
      write (unit=2, fmt='(a)') ' set autoscale x'
      write (unit=2, fmt='(a)') ' set autoscale y'
    end if
    !
    ! 3.6   Plot command
    !       ===========
    !
    if (type == 1)                                                      &
      then
      write (unit=2, fmt='(a,a,a,a,a,a,a,a,a)')                            &
        ' plot ''',file1,''' u 1:2 axes x1y1 title ''x''         w lp 12 ',', ',   &
        '''1.tmp'' u 1:4 axes x1y1 title ''trend(x) '' w l  11 ',', ',         &
        '''1.tmp'' u 1:3 axes x1y2 title ''y''         w lp  3 ',', ',           &
        '''1.tmp'' u 1:5 axes x1y2 title ''trend(y) '' w l   4 '
    else if (type == 2)                                                  &
      then
      write (unit=2, fmt='(a,a,a)')                                     &
        ' plot ''',file1,''' u 6:7 axes x1y1 notitle w p 3'
    end if
    !
    !       See notes to subroutine plotgn.
    !***************************************************************************
    !       Sun workstation, Unix, X11                    EPCF77 2.7
    !       <=> following lines uncommented
    !      write (unit=2, fmt='(1x,a)') 'pause -1'
    !***************************************************************************
    !
    !***************************************************************************
    !       PC, Suse Linux 8.0                            Portland PGF90 3.2-4
    !       <=> following lines uncommented
    !     write (unit=2, fmt='(1x,a)') 'pause -1'
    !***************************************************************************
    !
    !       IBM workstation, Windows XPP 64 bit           Absoft 10.0
    !       <=> following lines uncommented
    write (unit=2, fmt='(1x,a)') 'pause mouse key'
    close (unit=2, status ='keep')
    !
    ! 4.    Run gnuplot
    !       ==========
    !
    call plotgn
  end subroutine plot
  !
  !=============================================================================
  !
  subroutine plotgn
    !       Runs gnuplot.
    !       Notes:  (1)     Have Gnuplot (version 3.6 or higher) executable
    !                       in your path.
    !               (2)     Set, if applicable, x11 window attributes.
    !               (3)     Shell execution is non-standard to
    !                       Fortran 90. Therefore, some example
    !                       implementations are listed here.
    !
    !***************************************************************************
    !       IBM-PC, Windows 9x/NT, DOS window            AbSoft 6.2, AbSoft 10.0
    !       <=> following lines uncommented
    !         integer*4 plot, system
    !         external system
    !         plot=system( "gnuplot.exe 2.tmp" )
    !***************************************************************************
    !       IBM-PC, Linux                                 Absoft 6.0
    !       <=> following lines uncommented
    !       integer*4 plot, system
    !       external system
    !       plot=system( ' gnuplot 2.tmp ' )
    !***************************************************************************
    !        Sun workstation, Unix, X11                   EPCF77 2.7
    !       <=> following lines uncommented
    !        call system('/u2/users/mm/gnuplot/gnuplot
    !       +             -geometry -5-335 -background white
    !       +             /u2/users/mm/xtrend/2.tmp')
    ! **************************************************************************
    !
    !       IBM/Samsung/PC, Windows XP gfortran
    !
    !       <=> following lines uncommented
    !
    character(len=50) :: env
    integer :: ierr, stat
    logical :: windows
    call get_environment_variable("OS", env, status=ierr)
    windows = ierr == 0 .and. env == 'Windows_NT'
    if (windows) then
      call execute_command_line("where wgnuplot.exe > nul 2>&1", exitstat=stat)
    else
      call execute_command_line("which gnuplot > nul 2>&1", exitstat=stat)
    endif    
    if (stat /= 0) then
      print *,'Gnuplot not found'
      return
    end if
    if (windows) then
      call execute_command_line("wgnuplot.exe 2.tmp")
    else
      call execute_command_line("gnuplot 2.tmp")
    end if
  end subroutine plotgn
  !
  !=============================================================================
  !
  subroutine read1
    use precision
    use data1, only: t1,x1,y1
    use setting, only: datafile,n1
    implicit none
    !       Reads t1, x1, y1.
    character (len = 1) :: flag
    integer :: i
    open (unit=1, file=datafile, status='old',                      &
      form='formatted', action='read')
    do while (.true.)
      read (1, '(a1)') flag
      if (flag .ne. '#') then
        backspace (1)
        exit
      end if
    end do
    do i=1,n1
      read (unit=1, fmt=*) t1(i),x1(i),y1(i)
    end do
    close (unit=1, status='keep')
  end subroutine read1
  !
  !=============================================================================
  !
  subroutine r_est
    use precision
    use data2, only: t2,x2,y2,x3,y3
    use result1, only: r
    use setting, only: dtrtype,n2
    implicit none
    !       Estimates Pearson's correlation coefficient r(x2, y2). x2 and y2 are
    !       detrended (linearly or mean), renamed x3, y3, prior to estimation.
    integer :: i
    real(dp) :: a=-999.0_dp
    real(dp) :: b=-999.0_dp
    real(dp) :: mux=-999.0_dp        
    real(dp) :: muy=-999.0_dp
    if (dtrtype .eq. 'l') then
      call fit(t2,x2,n2,a,b)
      do i=1,n2
        x3(i)=x2(i)-(a+b*t2(i))
      end do
      call fit(t2,y2,n2,a,b)
      do i=1,n2
        y3(i)=y2(i)-(a+b*t2(i))
      end do
    else if (dtrtype .eq. 'm') then
      mux=sum(x2)/n2
      do i=1,n2
        x3(i)=x2(i)-mux
      end do
      muy=sum(y2)/n2
      do i=1,n2
        y3(i)=y2(i)-muy
      end do
    end if
    call pearsn(x3,y3,n2,r)
   
  end subroutine r_est
  !
  !=============================================================================
  !
  !
  !=============================================================================
  !
  subroutine tauest
    use data2, only: t2,x3,y3
    use result1, only: taux3,tauy3,rhox3,rhoy3
    use setting, only: n2
    use minls_module
    implicit none
    !
    call tauest_x(t2, x3, n2, taux3, rhox3)
    call tauest_x(t2, y3, n2, tauy3, rhoy3)
    !
  end subroutine tauest
  !
  !=============================================================================
  !
  !
  !=============================================================================
  !
  function t_inv(alpha,dof)
    use precision
    implicit none
    real(dp) :: t_inv
    !         Inverse Student's t distribution function, approximation
    !         after Abramowitz and Stegun (1964).
    !
    !         Highest inaccuracy occurs when:         o dof is small,
    !                                                       o alpha is high.
    !         For example: dof = 4, alpha = 0.025.
    !         There results: t = 2.7756 (instead of 2.776),
    !         that is, a negligible error.
    !
    real(dp),  intent(in) :: alpha
    integer,  intent(in) :: dof
    real(dp) :: u=-999.0_dp
    integer  :: d=0
    u=phi_inv(alpha)
    d=dof
    t_inv=u                                                            + &
      (1.0_dp/    4)*(                                u**3+    u) / d    + &
      (1.0_dp/   96)*(                    5*u**5+  16*u**3+  3*u) / d**2 + &
      (1.0_dp/  384)*(          3*u**7+  19*u**5+  17*u**3- 15*u) / d**3 + &
      (1.0_dp/92160)*(79*u**9+776*u**7+1482*u**5-1920*u**3-945*u) / d**4
  end function t_inv
end module pearsont3_module
