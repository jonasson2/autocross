! This file contains modules used by version 2 of redfit-x and pearsont3
! See comments at the top of refit-x.f90 and pearsont3.f90 for details.
!--------------------------------------------------------------------------
!         M O D U L E   D E F I N I T O N S
!--------------------------------------------------------------------------
! global array dimensions and constants
! -------------------------------------

module precision  
  integer, parameter :: dp = kind(1d0)
end module precision

module const
  use precision
  implicit none
  public
  character (len = 11) :: vers = 'REDFIT-X'
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: tpi = 2 * pi
  integer :: maxdimx, maxdimy    !the length of x and y
  integer :: nout                !nout = nfreq 
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
  real(dp),dimension(2) :: rhopre = -99  !1.rhopre x , 2. rhopre y
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
! error handling
! --------------
module error
  use precision
  implicit none
  public
  integer :: ierr = 0              ! error flag; 1=error in gettau; 2=warning in gettau
  integer, parameter :: errio = 99 ! i/o unit for log file
  logical :: errflagx = .false.    ! flags existence of duplicate times in input
  logical :: errflagy = .false.    ! flags existence of duplicate times in input
end module error
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

module meanvar_module
  use precision
contains
  subroutine meanvar(x, mean, var)
    ! Use Welford's algorithm, cf. Wikipedia
    implicit none
    real(dp), dimension(:), intent(in):: x
    real(dp), intent(out):: mean, var
    real(dp) delta, M2
    integer :: i, n
    n = size(x)
    mean = 0
    M2 = 0
    do i = 1, n
      delta = x(i) - mean
      mean = mean + delta/i
      M2 = M2 + delta*(x(i) - mean)
    end do
    if (n > 1) then
      var = M2 / (n - 1)
    else
      var = 0
    endif
  end subroutine meanvar
end module meanvar_module

module gammad_module
  interface
    function gammad(x, p, ifault)
      double precision, intent(in) :: x, p
      integer, intent(out) :: ifault
      double precision :: gammad
    end function gammad
  end interface
end module gammad_module

module random
  use precision
contains
  !-------------------------------------------------------------------------------------
  subroutine ranseed(seed)
  !-------------------------------------------------------------------------------------
    ! Seeds random number generator (based on the Fortran standard)
    implicit none
    integer, intent(in), optional :: seed
    integer, allocatable :: seedvector(:)
    integer :: n, i
    call random_seed(size = n)
    allocate(seedvector(n))
    if (present(seed)) then
      seedvector(1) = seed
    else
      seedvector(1) = -206761862
    end if
    do i=2,n
      seedvector(i) = 0
    end do
    call random_seed(put = seedvector)
  end subroutine ranseed

  function uniform_random() result(r)
    implicit none
    real(dp) :: r
    call random_number(r)
  end function uniform_random
  
  !---------------------------------------------------------------------------------------
  subroutine gasdev(harvest)
    !---------------------------------------------------------------------------------------
    implicit none
    real(dp), intent(out) :: harvest
    !     Gaussian distribution (Numerical Recipes, modified: uses ran).
    real(dp) :: rsq,v1,v2
    real(dp), save :: g
    logical, save :: gaus_stored=.false.
    if (gaus_stored) then
      harvest=g
      gaus_stored=.false.
    else
      do
        call random_number(v1)
        call random_number(v2)
        v1=2*v1-1
        v2=2*v2-1
        rsq=v1**2+v2**2
        if (rsq > 0 .and. rsq < 1) exit
      end do
      rsq=sqrt(-2*log(rsq)/rsq)
      harvest=v1*rsq
      g=v2*rsq
      gaus_stored=.true.
    end if
  end subroutine gasdev
end module random

!*******************************************************************************
!     
!  Module for [[fmin]] 1D derative-free function minimizer.
!     
!###  License
!  * [BSD-3](https://github.com/jacobwilliams/fmin/blob/master/LICENSE)
     
!  This module was obtained from https://github.com/jacobwilliams/fmin
!  and modified by Kristján Jónasson, April 2024 as follows:
!  a) changed the name to newbrent
!  b) removed ifdefs and parameterized the precision with dp
!  c) indented with emacs
!  d) added parameter cx, [ax, cx] is bracket and ax < bx < cx
!  e) added parameters xfunc, yfunc and nfunc that are passed to func
!  f) changed to return the min value and get the minimizer in parameter xmin
!
!  NOTE: According to a meld comparison of the two brent functions in
!  redfit-x.f90 and pearsont3.f90 in the 2016 version of the programs they
!  implement exactly the same algorithm, the latter using more modern Fortran.

module newbrent_module

  use iso_fortran_env
  use precision

  implicit none

  private

  ! integer,parameter,public :: fmin_rk = real64 !! real kind used by this
  ! module [8 bytes]
  ! integer,parameter :: wp = fmin_rk !! local copy of `fmin_rk` with a shorter
  ! name
  integer, parameter :: wp = dp  ! KJ modification

  abstract interface
    function func(x, xfunc, yfunc, nfunc) result(f)  ! KJ modification
      !     ! interface for user function
      import :: wp
      implicit none
      real(wp),intent(in) :: x  !! indep. variable
      real(wp)            :: f  !! function value `f(x)`
      integer,intent(in)  :: nfunc
      real(wp), dimension(:), intent(in) :: xfunc, yfunc
    end function func
  end interface

  public :: newbrent

contains
  !*****************************************************************************

  !*****************************************************************************
  !     >
  !     An approximation x to the point where `f` attains a minimum on
  !     the interval `(ax,bx)` is determined.
  !     
  !     the method used is a combination of golden section search and
  !     successive parabolic interpolation. convergence is never much slower
  !     than that for a fibonacci search. if `f` has a continuous second
  !     derivative which is positive at the minimum (which is not at `ax` or
  !     `bx`), then convergence is superlinear, and usually of the order of
  !     about 1.324.
  !     
  !     the function `f` is never evaluated at two points closer together
  !     than `eps*abs(fmin) + (tol/3)`, where `eps` is approximately the square
  !     root of the relative machine precision. if `f` is a unimodal
  !     function and the computed values of `f` are always unimodal when
  !     separated by at least `eps*abs(x) + (tol/3)`, then fmin approximates
  !     the abcissa of the global minimum of `f` on the interval `ax,bx` with
  !     an error less than `3*eps*abs(fmin) + tol`. if `f` is not unimodal,
  !     then `fmin` may approximate a local, but perhaps non-global, minimum to
  !     the same accuracy.
  !     
  !###  Reference
  !     * Richard brent, "algorithms for minimization without derivatives",
  !     prentice - hall, inc. (1973).
  !     
  !###  See also
  !     * [fmin from Netlib](http://www.netlib.org/fmm/fmin.f)

  function newbrent(f,ax,bx,cx,tol,xmin,xfunc,yfunc,nfunc) result(fxmin)
    ! KJ modification

    implicit none

    procedure(func)     :: f  !! the function to minimize
    real(wp),intent(in) :: ax !! left endpoint of initial interval
    real(wp),intent(in) :: bx !! initial guess, ax < bx < cx   - KJ modification
    real(wp),intent(in) :: cx !! right endpoint of initial interval  - KJ modif
    real(wp),intent(in) :: tol !! desired length of the interval of
    !                           ! uncertainty of the final result (>=0)
    integer,intent(in)  :: nfunc  ! KJ modification
    real(wp),intent(in), dimension(nfunc) :: xfunc, yfunc  ! KJ modification
    real(wp)            :: xmin !! abcissa approximating the point where
    !                            ! f attains a minimum
    real(wp)            :: fxmin !! function value at xmin  - KJ modification

    real(wp) :: a,b,d,e,xm,p,q,r,tol1,tol2,u,v,w
    real(wp) :: fu,fv,fw,fx,x

    real(wp),parameter :: c = (3.0_wp-sqrt(5.0_wp))/2.0_wp
    !! squared inverse of golden ratio
    real(wp),parameter :: half = 0.5_wp
    real(wp),parameter :: sqrteps = sqrt(epsilon(1.0_wp))

    !     initialization

    a = ax
    ! b = bx
    ! v = a + c*(b - a)
    b = cx  ! KJ modification
    v = bx  ! KJ modification
    w = v
    x = v
    e = 0.0_wp
    fx = f(x, xfunc, yfunc, nfunc)
    fv = fx
    fw = fx

    do                        !  main loop starts here

      xm = half*(a + b)
      tol1 = sqrteps*abs(x) + tol/3.0_wp
      tol2 = 2.0_wp*tol1

      !     check stopping criterion

      if (abs(x - xm) <= (tol2 - half*(b - a))) then
        !     write(*,*) 'x             = ', x
        !     write(*,*) 'xm            = ', xm
        !     write(*,*) 'abs(x - xm)   = ', abs(x - xm)
        !     write(*,*) 'tol2          = ', tol2
        !     write(*,*) 'half*(b - a)  = ', half*(b - a)
        exit
      end if

      !     is golden-section necessary

      if (abs(e) <= tol1) then

        !     a golden-section step

        if (x >= xm) then
          e = a - x
        else
          e = b - x
        end if
        d = c*e

      else

        !     fit parabola

        r = (x - w)*(fx - fv)
        q = (x - v)*(fx - fw)
        p = (x - v)*q - (x - w)*r
        q = 2.0_wp*(q - r)
        if (q > 0.0_wp) p = -p
        q =  abs(q)
        r = e
        e = d

        !     is parabola acceptable

        if ((abs(p) >= abs(half*q*r)) .or. (p <= q*(a - x)) .or. &
          (p >= q*(b - x))) then

          !     a golden-section step

          if (x >= xm) then
            e = a - x
          else
            e = b - x
          end if
          d = c*e

        else

          !     a parabolic interpolation step

          d = p/q
          u = x + d

          !     f must not be evaluated too close to ax or bx

          if (((u - a) < tol2) .or. ((b - u) < tol2)) d = sign(tol1, xm - x)

        end if

      end if

      !     f must not be evaluated too close to x

      if (abs(d) >= tol1) then
        u = x + d
      else
        u = x + sign(tol1, d)
      end if
      fu = f(u, xfunc, yfunc, nfunc)

      !     update a, b, v, w, and x

      if (fu <= fx) then
        if (u >= x) a = x
        if (u < x) b = x
        v = w
        fv = fw
        w = x
        fw = fx
        x = u
        fx = fu
      else
        if (u < x) a = u
        if (u >= x) b = u
        if (fu <= fw .or. w == x) then
          v = w
          fv = fw
          w = u
          fw = fu
        else if (fu <= fv .or. v == x .or. v == w ) then
          v = u
          fv = fu
        end if
      end if

    end do                    !  end of main loop

    xmin = x
    fxmin = fx  ! KJ modification

  end function newbrent
  !*****************************************************************************

  !*****************************************************************************
end module newbrent_module
!*******************************************************************************

module redfit_x_module
  use precision
  implicit none

contains
  !--------------------------------------------------------------------------
  subroutine allocerr(ichar)
    !--------------------------------------------------------------------------
    ! Report errors during allocation or deallocation to errio and terminate
    ! program.
    !
    ! ichar = a : Allocation Error
    ! ichar = d : Deallocation Error
    !--------------------------------------------------------------------------
    use error
    !
    implicit none
    !
    character(len=1), intent(in)  :: ichar
    !
    if (ichar .eq. "a") then
      write (errio, '(1x,a)') 'Error - Insufficient RAM; memory allocation failed!'
      close(errio)
      write(*,*) "An error has occurred. Check REDFIT-X.LOG for further information."
      stop
    end if
    if (ichar .eq. "d") then
      write (errio, '(1x,a)') 'Error - memory deallocation failed!'
      close(errio)
      write(*,*) "An error has occurred. Check REDFIT-X.LOG for further information."
      stop
    end if
    !
  end subroutine allocerr

  !--------------------------------------------------------------------------
  subroutine setdim(fnin,maxdimx,maxdimy,nout)  
    !--------------------------------------------------------------------------
    ! Analyze data file and set array dimensions.
    !--------------------------------------------------------------------------

    use timeser
    use param
    use error
    use nyquist
    use const, only : n_fnin, pi
    !
    implicit none
    !
    character (len = 80), dimension (n_fnin), intent(in) :: fnin 
    integer, intent(out) :: maxdimx,maxdimy,nout
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
      open (10, file = fnin(j), form = 'formatted', status = 'old', &
        iostat = iocheck)
      if (iocheck .ne. 0 ) then
        write (errio, *) ' Error - Can''t open ', trim(fnin(j))
        close(errio)
        stop
      end if
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
    ! set max. array dimension
    ! ------------------------
    maxdimx = npx
    maxdimy = npy   
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
    nout = nfreq
    !
    ! diagnostic output to stdout
    ! ---------------------------
    print '(a)'
    print '(a,f10.2)', "dtxy =",avgdtxy
    print '(a,i10)', "Nout =", nout
    print '(a,f10.2)', "fnyq =",fnyq
    print '(a)'
    !
  end subroutine setdim

  !--------------------------------------------------------------------------
  subroutine readdat(fnin)
    !--------------------------------------------------------------------------
    use timeser
    use error
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
      open (10, file = fnin(j), form = 'formatted', status = 'old', &
        iostat = iocheck)
      if (iocheck .ne. 0 ) then
        write (errio, *) ' Error - Can''t open ', trim(fnin(j))
        close(errio)
        stop
      end if
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
  subroutine check1()
    !--------------------------------------------------------------------------
    use timeser
    use error
    implicit none
    !
    real(dp)    :: ave
    integer :: i, j, idx, icnt
    logical :: err_order = .false.
    logical :: err_dupl = .false.
    !
    ! check for descending time axis
    ! ------------------------------
    do i = 1, npx-1
      if (tx(i) .gt. tx(i+1)) then
        err_order = .true.
        write(errio,'(1x,a,1x,e12.6)') 'Descending time axis at t =', tx(i)
      end if
    end do
    if (err_order .eqv. .true.) stop
    !
    ! check for duplicates time and replace values of time-dependent variable by their mean
    ! -------------------------------------------------------------------------------------
    idx = 1
    do while (.true.)
      if (tx(idx+1) .eq. tx(idx)) then
        err_dupl = .true.
        write(errio,'(1x,a,1x,e12.6,1x,a)') 'Duplicate time at t =', tx(idx), '...averaging data'
        !
        !       search for multiple occurrences
        !       -------------------------------
        icnt = 0
        do i = 1, npx-idx
          if (tx(idx+i) .eq. tx(idx)) icnt = icnt + 1
        end do
        !
        !       replace first data by mean of duplicates points
        !       -----------------------------------------------
        ave = sum(x(idx:idx+icnt)) / real(icnt+1, dp)
        x(idx) = ave
        !
        !       shift remaining points
        !       ----------------------
        do j = 1, icnt
          do i = idx+1, npx-1
            tx(i) = tx(i+1)
            x(i) = x(i+1)
          end do
          npx  = npx - 1
        end do
      end if
      idx = idx + 1
      if (idx .eq. npx-1) exit
    end do
    !
    ! save averaged data set
    ! ----------------------
    if (err_dupl .eqv. .true.) then
      errflagx = .true.
      open(90, file = "TimeSeriesx.avg")
      do i = 1, npx
        write(90,*) tx(i), x(i)
      end do
      close(90)
    end if
    !
  end subroutine check1

  !--------------------------------------------------------------------------
  subroutine check2()
    !--------------------------------------------------------------------------
    use timeser
    use error
    implicit none
    !
    real(dp)    :: ave
    integer :: i, j, idx, icnt
    logical :: err_order = .false.
    logical :: err_dupl = .false.
    !
    ! check for descending time axis
    ! ------------------------------
    do i = 1, npy-1
      if (ty(i) .gt. ty(i+1)) then
        err_order = .true.
        write(errio,'(1x,a,1x,e12.6)') 'Descending time axis at ty =', ty(i)
      end if
    end do
    if (err_order .eqv. .true.) stop
    !
    ! check for duplicates time and replace values of time-dependent variable by their mean
    ! -------------------------------------------------------------------------------------
    idx = 1
    do while (.true.)
      if (ty(idx+1) .eq. ty(idx)) then
        err_dupl = .true.
        write(errio,'(1x,a,1x,e12.6,1x,a)') 'Duplicate time at ty =', ty(idx), '...averaging data'
        !
        !       search for multiple occurrences
        !       -------------------------------
        icnt = 0
        do i = 1, npy-idx
          if (ty(idx+i) .eq. ty(idx)) icnt = icnt + 1
        end do
        !
        !       replace first data by mean of duplicates points
        !       -----------------------------------------------
        ave = sum(y(idx:idx+icnt)) / real(icnt+1, dp)
        y(idx) = ave
        !
        !       shift remaining points
        !       ----------------------
        do j = 1, icnt
          do i = idx+1, npy-1
            ty(i) = ty(i+1)
            y(i) = y(i+1)
          end do
          npy  = npy - 1
        end do
      end if
      idx = idx + 1
      if (idx .eq. npy-1) exit
    end do
    !
    ! save averaged data set
    ! ----------------------
    if (err_dupl .eqv. .true.) then
      errflagy = .true.
      open(90, file = "TimeSeriesy.avg")
      do i = 1, npy
        write(90,*) ty(i), y(i)
      end do
      close(90)
    end if
    !
  end subroutine check2

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
    complex(dp),dimension(nout) :: cpxy  
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
    allocate(txwk(nsegx), xwk(nsegx), tywk(nsegy), ywk(nsegy), stat = ialloc)
    if (ialloc .ne. 0) call allocerr("a")
    allocate(ftrx(lfreq), ftix(lfreq), ftry(lfreq), ftiy(lfreq), stat = ialloc)
    if (ialloc .ne. 0) call allocerr("a")
    if (ini .eqv. .true.) then
      allocate(txcos(nsegx,nfreq,n50),tycos(nsegy,nfreq,n50), stat = ialloc)
      if (ialloc .ne. 0) call allocerr("a")
      allocate(txsin(nsegx,nfreq,n50),tysin(nsegy,nfreq,n50), stat = ialloc)
      if (ialloc .ne. 0) call allocerr("a")
      allocate(wxtau(nfreq,n50),wytau(nfreq,n50), stat = ialloc)
      if (ialloc .ne. 0) call allocerr("a")
      allocate(wwx(nsegx,n50), wwy(nsegy,n50), stat = ialloc)
      if (ialloc .ne. 0) call allocerr("a")
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
      do j = 1, nout
        gxx(j) = gxx(j) + (ftrx(j)*ftrx(j) + ftix(j)*ftix(j))
        gyy(j) = gyy(j) + (ftry(j)*ftry(j) + ftiy(j)*ftiy(j))
      end do
      !
      !    cross and phase spectra
      !    --------------------------
      do j = 1,nout
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
    do i = 1, nout
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
      do i=1,nout
        if(gxy(i)==0 .or. gxx(i)==0 .or. gyy(i)==0) then
          cxy(i) = 0
        else 
          cxy(i)= (gxy(i) * gxy(i))/(gxx(i) * gyy(i))
        end if
      end do
    end if
    !
    deallocate(txwk, xwk, tywk, ywk, stat = ialloc)
    if (ialloc .ne. 0) call allocerr("d")
    deallocate(ftrx, ftix, ftry, ftiy, stat = ialloc)
    if (ialloc .ne. 0) call allocerr("d")
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
    use error
    use gammad_module
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
      if (ierr .eq. 1) return
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
          write(errio,'(a)') "Error in GETCHI2: Iter > ItMax"
          ierr = 1
          return
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
    use error
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
          write(errio,'(a)') "Error in GETZ: Iter > ItMax"
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
  subroutine gettau(rhopre,tx,x,npx,tau)
    !----------------------------------------------------------------------
    use const
    use param, only : n50
    use error
    implicit none
    real(dp), intent(in) :: rhopre
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
    !if prescribed tau
    if(rhopre .ge.0) then
      tau = -avgdt / log(rhopre)
    else if (rhopre.lt.0)then  
      rhosum = 0
      nseg = int(2 * npx / (n50 + 1))         ! points per segment
      allocate(twk(nseg), xwk(nseg), stat = ialloc)
      if (ialloc .ne. 0) call allocerr("a")
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
        call tauest(twk(1:nseg), xwk(1:nseg), nseg, tau, rho)
        if (ierr .eq. 1) then
          write(errio,*) ' Error in TAUEST'
          return
        end if
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
        ierr = 2
        write (errio,*) 'Warning: GETTAU returned tau =', tau
        write (errio,*) '         Negative tau is forced to zero.'
        tau = 0
      end if
      !
      deallocate(twk, xwk, stat = ialloc)
      if (ialloc .ne. 0) call allocerr("d")
      ! 
    end if
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
  subroutine tauest(t, x, np, tau, rhoavg)
    !
    use const
    use error
    use meanvar_module
    !
    implicit none
    !
    integer, intent(in) :: np
    real(dp), dimension(np), intent(in) :: t, x
    real(dp), intent (out)  :: tau
    real(dp), dimension(np) :: tscal, xscal
    real(dp) :: fac, avg, var, dt, rho, scalt, amin, rhoavg
    real(dp) :: damin
    integer :: i, mult
    !
    ! Correct time direction; assume that ages are input
    ! --------------------------------------------------
    do i = 1, np
      tscal(i) = -t(np+1-i)
      xscal(i) = x(np+1-i)
    end do
    !
    ! Scaling of x
    ! ------------
    call meanvar(xscal(1:np), avg, var)
    fac = sqrt(var)
    xscal(1:np) = xscal(1:np) / fac
    !
    ! Scaling of t (=> start value of a = 1/e)
    ! ---------------------------------------
    dt = (tscal(np)-tscal(1)) / real(np-1, dp)
    call rhoest(np, xscal(1:np), rho)
    if (rho .le. 0) then
      rho = 0.05_dp
      write(errio,*) 'Warning: rho estimation: < 0'
      ierr = 2
    else if (rho .gt. 1) then
      rho = 0.95_dp
      write(errio,*) 'Warning: rho estimation: > 1'
      ierr = 2
    end if
    scalt = -log(rho)/dt
    tscal(1:np) = tscal(1:np) * scalt
    !
    ! Estimation
    ! ----------
    call minls(np, tscal(1:np), xscal(1:np), damin, mult)
    if (ierr .eq. 1) then
      write(errio,*) ' Error in MNILS'
      return
    end if
    amin = sngl(damin)
    if (mult .eq. 1) then
      write(errio,*) ' Estimation problem: LS function has > 1 minima'
      return
    end if
    if (amin .le. 0) then
      write(errio,*) ' Estimation problem: a_min =< 0'
      return
    else if (amin .ge. 1) then
      write(errio,*) ' Estimation problem: a_min >= 1'
      return
    end if
    !
    ! determine tau
    ! -------------
    tau = -1 /(scalt*log(amin))
    !
    ! determine rho, corresponding to tau
    ! -----------------------------------
    rhoavg = exp(-dt / tau)
    !
  end subroutine tauest

  ! Numerical Recipes (modified): Brent's search in one direction:
  ! function brent(ax,bx,cx,f,tol,xmin,xfunc,yfunc,nfunc)

  !----------------------------------------------------------------------
  ! Least-squares function
  !----------------------------------------------------------------------
  real(dp) function ls(a,t,x,n)
    implicit none
    real(dp), intent(in) :: a, t(:), x(:)
    integer, intent(in) :: n
    integer i
    ls=0.0_dp
    do i=2,n
      ls=ls+(x(i)-x(i-1)*sign(1.0_dp,a)* abs(a)**(t(i)-t(i-1)))**2
    end do
    return
  end function ls

  !----------------------------------------------------------------------
  ! Minimization of least-squares function ls.
  !----------------------------------------------------------------------
  subroutine minls(n, t, x, amin, nmu_)
    !
    use error
    use newbrent_module
    !
    implicit none
    !
    real(dp), parameter :: a_ar1 = 0.367879441_dp ! 1/e
    real(dp), parameter :: tol = 3e-8_dp           ! Brent's search, precision
    real(dp), parameter :: tol2 = 1e-6_dp          ! multiple solutions, precision
    integer n
    real(dp) t(1:n),x(1:n)
    real(dp) amin
    integer nmu_
    real(dp) dum1,dum2,dum3,dum4,a_ar11,a_ar12,a_ar13
    !
    nmu_=0
    dum1=newbrent(ls, -2.0_dp, a_ar1, +2.0_dp, tol, a_ar11, t, x, n)
    dum2=newbrent(ls, a_ar1, 0.5_dp*(a_ar1+1.0_dp), +2.0_dp, tol, a_ar12, t, x, n)
    dum3=newbrent(ls,-2.0_dp, 0.5_dp*(a_ar1-1.0_dp),  a_ar1, tol, a_ar13, t, x, n)
    if (ierr .eq. 1) then
      write(errio, *) ' Error in MINLS (call to brent)'
      return
    end if
    if  ((abs(a_ar12-a_ar11).gt.tol2.and.abs(a_ar12-a_ar1).gt.tol2) &
      .or.(abs(a_ar13-a_ar11).gt.tol2.and.abs(a_ar13-a_ar1).gt.tol2)) &
      nmu_=1
    dum4=min(dum1,dum2,dum3)
    if (dum4.eq.dum2) then
      amin=a_ar12
    else if (dum4.eq.dum3) then
      amin=a_ar13
    else
      amin=a_ar11
    end if
    return
  end subroutine minls

  !----------------------------------------------------------------------
  ! Autocorrelation coefficient estimation (equidistant data).
  !----------------------------------------------------------------------
  subroutine rhoest(n,x,rho)
    !
    implicit none
    !
    integer n
    real(dp) x(1:n)
    real(dp) rho
    integer i
    real(dp) sum1,sum2
    !
    sum1=0
    sum2=0
    do i=2,n
      sum1=sum1+x(i)*x(i-1)
      sum2=sum2+x(i)**2
    end do
    rho=sum1/sum2
    return
  end subroutine rhoest

end module redfit_x_module

! The following functions downto and including gammad are copied directly from:
! https://people.sc.fsu.edu/~jburkardt/f_src/asa239/asa239.f90
  
function alngam ( xvalue, ifault )

!*****************************************************************************80
!
!! alngam() computes the logarithm of the gamma function.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by Allan Macleod.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Allan Macleod,
!    Algorithm AS 245,
!    A Robust and Reliable Algorithm for the Logarithm of the Gamma Function,
!    Applied Statistics,
!    Volume 38, Number 2, 1989, pages 397-402.
!
!  Parameters:
!
!    Input, real ( kind = rk ) XVALUE, the argument of the Gamma function.
!
!    Output, integer IFAULT, error flag.
!    0, no error occurred.
!    1, XVALUE is less than or equal to 0.
!    2, XVALUE is too big.
!
!    Output, real ( kind = rk ) ALNGAM, the logarithm of the gamma function of X
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) alngam
  real ( kind = rk ), parameter :: alr2pi = 0.918938533204673D+00
  integer ifault
  real ( kind = rk ), dimension ( 9 ) :: r1 = (/ &
    -2.66685511495D+00, &
    -24.4387534237D+00, &
    -21.9698958928D+00, &
     11.1667541262D+00, &
     3.13060547623D+00, &
     0.607771387771D+00, &
     11.9400905721D+00, &
     31.4690115749D+00, &
     15.2346874070D+00 /)
  real ( kind = rk ), dimension ( 9 ) :: r2 = (/ &
    -78.3359299449D+00, &
    -142.046296688D+00, &
     137.519416416D+00, &
     78.6994924154D+00, &
     4.16438922228D+00, &
     47.0668766060D+00, &
     313.399215894D+00, &
     263.505074721D+00, &
     43.3400022514D+00 /)
  real ( kind = rk ), dimension ( 9 ) :: r3 = (/ &
    -2.12159572323D+05, &
     2.30661510616D+05, &
     2.74647644705D+04, &
    -4.02621119975D+04, &
    -2.29660729780D+03, &
    -1.16328495004D+05, &
    -1.46025937511D+05, &
    -2.42357409629D+04, &
    -5.70691009324D+02 /)
  real ( kind = rk ), dimension ( 5 ) :: r4 = (/ &
     0.279195317918525D+00, &
     0.4917317610505968D+00, &
     0.0692910599291889D+00, &
     3.350343815022304D+00, &
     6.012459259764103D+00 /)
  real ( kind = rk ) x
  real ( kind = rk ) x1
  real ( kind = rk ) x2
  real ( kind = rk ), parameter :: xlge = 5.10D+05
  real ( kind = rk ), parameter :: xlgst = 1.0D+30
  real ( kind = rk ) xvalue
  real ( kind = rk ) y

  x = xvalue
  alngam = 0.0D+00
!
!  Check the input.
!
  if ( xlgst <= x ) then
    ifault = 2
    return
  end if

  if ( x <= 0.0D+00 ) then
    ifault = 1
    return
  end if

  ifault = 0
!
!  Calculation for 0 < X < 0.5 and 0.5 <= X < 1.5 combined.
!
  if ( x < 1.5D+00 ) then

    if ( x < 0.5D+00 ) then

      alngam = - log ( x )
      y = x + 1.0D+00
!
!  Test whether X < machine epsilon.
!
      if ( y == 1.0D+00 ) then
        return
      end if

    else

      alngam = 0.0D+00
      y = x
      x = ( x - 0.5D+00 ) - 0.5D+00

    end if

    alngam = alngam + x * (((( &
        r1(5)   * y &
      + r1(4) ) * y &
      + r1(3) ) * y &
      + r1(2) ) * y &
      + r1(1) ) / (((( &
                  y &
      + r1(9) ) * y &
      + r1(8) ) * y &
      + r1(7) ) * y &
      + r1(6) )

    return

  end if
!
!  Calculation for 1.5 <= X < 4.0.
!
  if ( x < 4.0D+00 ) then

    y = ( x - 1.0D+00 ) - 1.0D+00

    alngam = y * (((( &
        r2(5)   * x &
      + r2(4) ) * x &
      + r2(3) ) * x &
      + r2(2) ) * x &
      + r2(1) ) / (((( &
                  x &
      + r2(9) ) * x &
      + r2(8) ) * x &
      + r2(7) ) * x &
      + r2(6) )
!
!  Calculation for 4.0 <= X < 12.0.
!
  else if ( x < 12.0D+00 ) then

    alngam = (((( &
        r3(5)   * x &
      + r3(4) ) * x &
      + r3(3) ) * x &
      + r3(2) ) * x &
      + r3(1) ) / (((( &
                  x &
      + r3(9) ) * x &
      + r3(8) ) * x &
      + r3(7) ) * x &
      + r3(6) )
!
!  Calculation for 12.0 <= X.
!
  else

    y = log ( x )
    alngam = x * ( y - 1.0D+00 ) - 0.5D+00 * y + alr2pi

    if ( x <= xlge ) then

      x1 = 1.0D+00 / x
      x2 = x1 * x1

      alngam = alngam + x1 * ( ( &
             r4(3)   * &
        x2 + r4(2) ) * &
        x2 + r4(1) ) / ( ( &
        x2 + r4(5) ) * &
        x2 + r4(4) )

    end if

  end if

  return
end
function alnorm ( x, upper )

!*****************************************************************************80
!
!! alnorm() computes the cumulative density of the standard normal distribution.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by David Hill.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Hill,
!    Algorithm AS 66:
!    The Normal Integral,
!    Applied Statistics,
!    Volume 22, Number 3, 1973, pages 424-427.
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, is one endpoint of the semi-infinite interval
!    over which the integration takes place.
!
!    Input, logical UPPER, determines whether the upper or lower
!    interval is to be integrated:
!    .TRUE.  => integrate from X to + Infinity;
!    .FALSE. => integrate from - Infinity to X.
!
!    Output, real ( kind = rk ) ALNORM, the integral of the standard normal
!    distribution over the desired interval.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), parameter :: a1 = 5.75885480458D+00
  real ( kind = rk ), parameter :: a2 = 2.62433121679D+00
  real ( kind = rk ), parameter :: a3 = 5.92885724438D+00
  real ( kind = rk ) alnorm
  real ( kind = rk ), parameter :: b1 = -29.8213557807D+00
  real ( kind = rk ), parameter :: b2 = 48.6959930692D+00
  real ( kind = rk ), parameter :: c1 = -0.000000038052D+00
  real ( kind = rk ), parameter :: c2 = 0.000398064794D+00
  real ( kind = rk ), parameter :: c3 = -0.151679116635D+00
  real ( kind = rk ), parameter :: c4 = 4.8385912808D+00
  real ( kind = rk ), parameter :: c5 = 0.742380924027D+00
  real ( kind = rk ), parameter :: c6 = 3.99019417011D+00
  real ( kind = rk ), parameter :: con = 1.28D+00
  real ( kind = rk ), parameter :: d1 = 1.00000615302D+00
  real ( kind = rk ), parameter :: d2 = 1.98615381364D+00
  real ( kind = rk ), parameter :: d3 = 5.29330324926D+00
  real ( kind = rk ), parameter :: d4 = -15.1508972451D+00
  real ( kind = rk ), parameter :: d5 = 30.789933034D+00
  real ( kind = rk ), parameter :: ltone = 7.0D+00
  real ( kind = rk ), parameter :: p = 0.398942280444D+00
  real ( kind = rk ), parameter :: q = 0.39990348504D+00
  real ( kind = rk ), parameter :: r = 0.398942280385D+00
  logical up
  logical upper
  real ( kind = rk ), parameter :: utzero = 18.66D+00
  real ( kind = rk ) x
  real ( kind = rk ) y
  real ( kind = rk ) z

  up = upper
  z = x

  if ( z < 0.0D+00 ) then
    up = .not. up
    z = - z
  end if

  if ( ltone < z .and. ( ( .not. up ) .or. utzero < z ) ) then

    if ( up ) then
      alnorm = 0.0D+00
    else
      alnorm = 1.0D+00
    end if

    return

  end if

  y = 0.5D+00 * z * z

  if ( z <= con ) then

    alnorm = 0.5D+00 - z * ( p - q * y &
      / ( y + a1 + b1 &
      / ( y + a2 + b2 & 
      / ( y + a3 ))))

  else

    alnorm = r * exp ( - y ) &
      / ( z + c1 + d1 &
      / ( z + c2 + d2 &
      / ( z + c3 + d3 &
      / ( z + c4 + d4 &
      / ( z + c5 + d5 &
      / ( z + c6 ))))))

  end if

  if ( .not. up ) then
    alnorm = 1.0D+00 - alnorm
  end if

  return
end

subroutine gamma_inc_values ( n_data, a, x, fx )

!*****************************************************************************80
!
!! gamma_inc_values() returns some values of the incomplete Gamma function.
!
!  Discussion:
!
!    The (normalized) incomplete Gamma function P(A,X) is defined as:
!
!      PN(A,X) = 1/Gamma(A) * Integral ( 0 <= T <= X ) T^(A-1) * exp(-T) dT.
!
!    With this definition, for all A and X,
!
!      0 <= PN(A,X) <= 1
!
!    and
!
!      PN(A,INFINITY) = 1.0
!
!    In Mathematica, the function can be evaluated by:
!
!      1 - GammaRegularized[A,X]
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = rk ) A, the parameter of the function.
!
!    Output, real ( kind = rk ) X, the argument of the function.
!
!    Output, real ( kind = rk ) FX, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 20

  real ( kind = rk ) a
  real ( kind = rk ), save, dimension ( n_max ) :: a_vec = (/ &
    0.10D+00, &
    0.10D+00, &
    0.10D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.10D+01, &
    0.10D+01, &
    0.10D+01, &
    0.11D+01, &
    0.11D+01, &
    0.11D+01, &
    0.20D+01, &
    0.20D+01, &
    0.20D+01, &
    0.60D+01, &
    0.60D+01, &
    0.11D+02, &
    0.26D+02, &
    0.41D+02 /)
  real ( kind = rk ) fx
  real ( kind = rk ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.7382350532339351D+00, &
    0.9083579897300343D+00, &
    0.9886559833621947D+00, &
    0.3014646416966613D+00, &
    0.7793286380801532D+00, &
    0.9918490284064973D+00, &
    0.9516258196404043D-01, &
    0.6321205588285577D+00, &
    0.9932620530009145D+00, &
    0.7205974576054322D-01, &
    0.5891809618706485D+00, &
    0.9915368159845525D+00, &
    0.1018582711118352D-01, &
    0.4421745996289254D+00, &
    0.9927049442755639D+00, &
    0.4202103819530612D-01, &
    0.9796589705830716D+00, &
    0.9226039842296429D+00, &
    0.4470785799755852D+00, &
    0.7444549220718699D+00 /)
  integer n_data
  real ( kind = rk ) x
  real ( kind = rk ), save, dimension ( n_max ) :: x_vec = (/ &
    0.30D-01, &
    0.30D+00, &
    0.15D+01, &
    0.75D-01, &
    0.75D+00, &
    0.35D+01, &
    0.10D+00, &
    0.10D+01, &
    0.50D+01, &
    0.10D+00, & 
    0.10D+01, &
    0.50D+01, &
    0.15D+00, &
    0.15D+01, &
    0.70D+01, &
    0.25D+01, &
    0.12D+02, &
    0.16D+02, &
    0.25D+02, &
    0.45D+02 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function gammad ( x, p, ifault )

!*****************************************************************************80
!
!! gammad() computes the Incomplete Gamma Integral
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by B Shea.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    B Shea,
!    Algorithm AS 239:
!    Chi-squared and Incomplete Gamma Integral,
!    Applied Statistics,
!    Volume 37, Number 3, 1988, pages 466-473.
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, P, the parameters of the incomplete 
!    gamma ratio.  0 <= X, and 0 < P.
!
!    Output, integer IFAULT, error flag.
!    0, no error.
!    1, X < 0 or P <= 0.
!
!    Output, real ( kind = rk ) GAMMAD, the value of the incomplete 
!    Gamma integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) alnorm
  real ( kind = rk ) alngam
  real ( kind = rk ) an
  real ( kind = rk ) arg
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ), parameter :: elimit = - 88.0D+00
  real ( kind = rk ) gammad
  integer, intent(out) :: ifault  ! KJ added intent
  real ( kind = rk ), parameter :: oflo = 1.0D+37
  real ( kind = rk ), intent(in) :: p  ! KJ added intent
  real ( kind = rk ), parameter :: plimit = 1000.0D+00
  real ( kind = rk ) pn1
  real ( kind = rk ) pn2
  real ( kind = rk ) pn3
  real ( kind = rk ) pn4
  real ( kind = rk ) pn5
  real ( kind = rk ) pn6
  real ( kind = rk ) rn
  real ( kind = rk ), parameter :: tol = 1.0D-14
  logical upper
  real ( kind = rk ), intent(in) :: x  ! KJ added intent
  real ( kind = rk ), parameter :: xbig = 1.0D+08

  gammad = 0.0D+00
!
!  Check the input.
!
  if ( x < 0.0D+00 ) then
    ifault = 1
    return
  end if

  if ( p <= 0.0D+00 ) then
    ifault = 1
    return
  end if

  ifault = 0

  if ( x == 0.0D+00 ) then
    gammad = 0.0D+00
    return
  end if
!
!  If P is large, use a normal approximation.
!
  if ( plimit < p ) then

    pn1 = 3.0D+00 * sqrt ( p ) * ( ( x / p )**( 1.0D+00 / 3.0D+00 ) &
    + 1.0D+00 / ( 9.0D+00 * p ) - 1.0D+00 )

    upper = .false.
    gammad = alnorm ( pn1, upper )
    return

  end if
!
!  If X is large set GAMMAD = 1.
!
  if ( xbig < x ) then
    gammad = 1.0D+00
    return
  end if
!
!  Use Pearson's series expansion.
!  (Note that P is not large enough to force overflow in ALOGAM).
!  No need to test IFAULT on exit since P > 0.
!
  if ( x <= 1.0D+00 .or. x < p ) then

    arg = p * log ( x ) - x - alngam ( p + 1.0D+00, ifault )
    c = 1.0D+00
    gammad = 1.0D+00
    a = p

    do

      a = a + 1.0D+00
      c = c * x / a
      gammad = gammad + c

      if ( c <= tol ) then
        exit
      end if

    end do

    arg = arg + log ( gammad )

    if ( elimit <= arg ) then
      gammad = exp ( arg )
    else
      gammad = 0.0D+00
    end if
!
!  Use a continued fraction expansion.
!
  else 

    arg = p * log ( x ) - x - alngam ( p, ifault )
    a = 1.0D+00 - p
    b = a + x + 1.0D+00
    c = 0.0D+00
    pn1 = 1.0D+00
    pn2 = x
    pn3 = x + 1.0D+00
    pn4 = x * b
    gammad = pn3 / pn4

    do

      a = a + 1.0D+00
      b = b + 2.0D+00
      c = c + 1.0D+00
      an = a * c
      pn5 = b * pn3 - an * pn1
      pn6 = b * pn4 - an * pn2

      if ( pn6 /= 0.0D+00 ) then

        rn = pn5 / pn6

        if ( abs ( gammad - rn ) <= min ( tol, tol * rn ) ) then
          exit
        end if

        gammad = rn

      end if

      pn1 = pn3
      pn2 = pn4
      pn3 = pn5
      pn4 = pn6
!
!  Re-scale terms in continued fraction if terms are large.
!
      if ( oflo <= abs ( pn5 ) ) then
        pn1 = pn1 / oflo
        pn2 = pn2 / oflo
        pn3 = pn3 / oflo
        pn4 = pn4 / oflo
      end if

    end do

    arg = arg + log ( gammad )

    if ( elimit <= arg ) then
      gammad = 1.0D+00 - exp ( arg )
    else
      gammad = 1.0D+00
    end if

  end if

  return
end
!-------- End of functions copied from J.Burkardt's home page -------

module sort_module
  use precision
  implicit none

contains

  logical function lsame(ca, cb)
    ! A simplified version of the netlib lsame function that assumes ASCII
    character, intent(in) :: ca, cb
    integer ia, ib
    lsame = ca == cb 
    if (.not. lsame) then
      ia = ichar(ca)
      ib = ichar(cb)
      if (ia > 96 .and. ia < 123) ia = ia - 32
      if (ib > 96 .and. ib < 123) ib = ib - 32
      lsame = ia == ib
    end if
  end function lsame

  ! The following subroutine was converted from the netlib dlasrt
  ! using f77_to_f90 obtained from people.math.sc.edu/Burkardt.
  ! The only other change is commenting out lsame declaration
  ! and all references to xerbla, as well as changing double
  ! precision to real(dp)
  SUBROUTINE DLASRT (ID, N, D, INFO) 
    !                                                                       
    !  -- LAPACK routine (version 3.1) --                                   
    !     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..    
    !     November 2006                                                     
    !                                                                       
    !     .. Scalar Arguments ..                                            
    CHARACTER ID 
    INTEGER INFO, N 
    !     ..                                                                
    !     .. Array Arguments ..                                             
    REAL(DP) D ( * ) 
    !     ..                                                                
    !                                                                       
    !  Purpose                                                              
    !  =======                                                              
    !                                                                       
    !  Sort the numbers in D in increasing order (if ID = 'I') or           
    !  in decreasing order (if ID = 'D' ).                                  
    !                                                                       
    !  Use Quick Sort, reverting to Insertion sort on arrays of             
    !  size <= 20. Dimension of STACK limits N to about 2**32.              
    !                                                                       
    !  Arguments                                                            
    !  =========                                                            
    !                                                                       
    !  ID      (input) CHARACTER*1                                          
    !          = 'I': sort D in increasing order;                           
    !          = 'D': sort D in decreasing order.                           
    !                                                                       
    !  N       (input) INTEGER                                              
    !          The length of the array D.                                   
    !                                                                       
    !  D       (input/output) DOUBLE PRECISION array, dimension (N)         
    !          On entry, the array to be sorted.                            
    !          On exit, D has been sorted into increasing order             
    !          (D(1) <= ... <= D(N) ) or into decreasing order              
    !          (D(1) >= ... >= D(N) ), depending on ID.                     
    !                                                                       
    !  INFO    (output) INTEGER                                             
    !          = 0:  successful exit                                        
    !          < 0:  if INFO = -i, the i-th argument had an illegal value   
    !                                                                       
    !  =====================================================================
    !                                                                       
    !     .. Parameters ..                                                  
    INTEGER SELECT 
    PARAMETER (SELECT = 20) 
    !     ..                                                                
    !     .. Local Scalars ..                                               
    INTEGER DIR, ENDD, I, J, START, STKPNT 
    REAL(DP) D1, D2, D3, DMNMX, TMP 
    !     ..                                                                
    !     .. Local Arrays ..                                                
    INTEGER STACK (2, 32) 
    !     ..                                                                
    !     .. External Functions ..                                          
    ! LOGICAL LSAME 
    ! EXTERNAL LSAME 
    !     ..                                                                
    !     .. External Subroutines ..                                        
    ! EXTERNAL XERBLA 
    !     ..                                                                
    !     .. Executable Statements ..                                       
    !                                                                       
    !     Test the input paramters.                                         
    !                                                                       
    INFO = 0 
    DIR = - 1 
    IF (LSAME (ID, 'D') ) THEN 
      DIR = 0 
    ELSEIF (LSAME (ID, 'I') ) THEN 
      DIR = 1 
    ENDIF
    IF (DIR.EQ. - 1) THEN 
      INFO = - 1 
    ELSEIF (N.LT.0) THEN 
      INFO = - 2 
    ENDIF
    IF (INFO.NE.0) THEN 
      ! CALL XERBLA ('DLASRT', - INFO) 
      RETURN 
    ENDIF
    !                                                                       
    !     Quick return if possible                                          
    !                                                                       
    IF (N.LE.1) RETURN 
    !                                                                       
    STKPNT = 1 
    STACK (1, 1) = 1 
    STACK (2, 1) = N 
10  CONTINUE 
    START = STACK (1, STKPNT) 
    ENDD = STACK (2, STKPNT) 
    STKPNT = STKPNT - 1 
    IF (ENDD-START.LE.SELECT.AND.ENDD-START.GT.0) THEN 
      !                                                                       
      !        Do Insertion sort on D( START:ENDD )                           
      !                                                                       
      IF (DIR.EQ.0) THEN 
        !                                                                       
        !           Sort into decreasing order                                  
        !                                                                       
        DO 30 I = START + 1, ENDD 
          DO 20 J = I, START + 1, - 1 
            IF (D (J) .GT.D (J - 1) ) THEN 
              DMNMX = D (J) 
              D (J) = D (J - 1) 
              D (J - 1) = DMNMX 
            ELSE 
              GOTO 30 
            ENDIF
20        END DO
30      END DO
        !                                                                       
      ELSE 
        !                                                                       
        !           Sort into increasing order                                  
        !                                                                       
        DO 50 I = START + 1, ENDD 
          DO 40 J = I, START + 1, - 1 
            IF (D (J) .LT.D (J - 1) ) THEN 
              DMNMX = D (J) 
              D (J) = D (J - 1) 
              D (J - 1) = DMNMX 
            ELSE 
              GOTO 50 
            ENDIF
40        END DO
50      END DO
        !                                                                       
      ENDIF
      !                                                                       
    ELSEIF (ENDD-START.GT.SELECT) THEN 
      !                                                                       
      !        Partition D( START:ENDD ) and stack parts, largest one first   
      !                                                                       
      !        Choose partition entry as median of 3                          
      !                                                                       
      D1 = D (START) 
      D2 = D (ENDD) 
      I = (START + ENDD) / 2 
      D3 = D (I) 
      IF (D1.LT.D2) THEN 
        IF (D3.LT.D1) THEN 
          DMNMX = D1 
        ELSEIF (D3.LT.D2) THEN 
          DMNMX = D3 
        ELSE 
          DMNMX = D2 
        ENDIF
      ELSE 
        IF (D3.LT.D2) THEN 
          DMNMX = D2 
        ELSEIF (D3.LT.D1) THEN 
          DMNMX = D3 
        ELSE 
          DMNMX = D1 
        ENDIF
      ENDIF
      !                                                                       
      IF (DIR.EQ.0) THEN 
        !                                                                       
        !           Sort into decreasing order                                  
        !                                                                       
        I = START - 1 
        J = ENDD+1 
60      CONTINUE 
70      CONTINUE 
        J = J - 1 
        IF (D (J) .LT.DMNMX) GOTO 70 
80      CONTINUE 
        I = I + 1 
        IF (D (I) .GT.DMNMX) GOTO 80 
        IF (I.LT.J) THEN 
          TMP = D (I) 
          D (I) = D (J) 
          D (J) = TMP 
          GOTO 60 
        ENDIF
        IF (J - START.GT.ENDD-J - 1) THEN 
          STKPNT = STKPNT + 1 
          STACK (1, STKPNT) = START 
          STACK (2, STKPNT) = J 
          STKPNT = STKPNT + 1 
          STACK (1, STKPNT) = J + 1 
          STACK (2, STKPNT) = ENDD 
        ELSE 
          STKPNT = STKPNT + 1 
          STACK (1, STKPNT) = J + 1 
          STACK (2, STKPNT) = ENDD 
          STKPNT = STKPNT + 1 
          STACK (1, STKPNT) = START 
          STACK (2, STKPNT) = J 
        ENDIF
      ELSE 
        !                                                                       
        !           Sort into increasing order                                  
        !                                                                       
        I = START - 1 
        J = ENDD+1 
90      CONTINUE 
100     CONTINUE 
        J = J - 1 
        IF (D (J) .GT.DMNMX) GOTO 100 
110     CONTINUE 
        I = I + 1 
        IF (D (I) .LT.DMNMX) GOTO 110 
        IF (I.LT.J) THEN 
          TMP = D (I) 
          D (I) = D (J) 
          D (J) = TMP 
          GOTO 90 
        ENDIF
        IF (J - START.GT.ENDD-J - 1) THEN 
          STKPNT = STKPNT + 1 
          STACK (1, STKPNT) = START 
          STACK (2, STKPNT) = J 
          STKPNT = STKPNT + 1 
          STACK (1, STKPNT) = J + 1 
          STACK (2, STKPNT) = ENDD 
        ELSE 
          STKPNT = STKPNT + 1 
          STACK (1, STKPNT) = J + 1 
          STACK (2, STKPNT) = ENDD 
          STKPNT = STKPNT + 1 
          STACK (1, STKPNT) = START 
          STACK (2, STKPNT) = J 
        ENDIF
      ENDIF
    ENDIF
    IF (STKPNT.GT.0) GOTO 10 
    RETURN 
    !                                                                       
    !     End of DLASRT                                                     
    !                                                                       
  END SUBROUTINE DLASRT

  subroutine SORT(D)
    double precision, dimension(:), intent(inout) :: D
    integer :: INFO

    ! Call DLASRT with 'I' for increasing order.
    call DLASRT('I', size(D), D, INFO)

    if (INFO /= 0) then
      print *, 'Error in DLASRT: INFO=', INFO
      stop
    endif
  end subroutine SORT

end module sort_module
