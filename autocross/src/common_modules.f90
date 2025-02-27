! This file contains modules introduced in 2024–2025 to replace Numerical
! recipes subroutines (see comments at the top of refit-x.f90 and pearsont3.f90
! for details).
!
! It also has minls_module with routines minls, tauest_x, rhoest, and lstau take
! from pearsont3, but also used by redfitx. The corresponding routines in
! redfitx have been deleted (tauest_x was tauest and lstau was ls in redfitx).

! --------------------------------------------------------------------------
! M O D U L E   D E F I N I T O N S
! --------------------------------------------------------------------------
! global array dimensions and constants
! -------------------------------------

module precision  
  integer, parameter :: dp = kind(1d0)
end module precision

module time
  real :: start_time
contains
  subroutine tic
    call cpu_time(start_time)
  end subroutine tic

  subroutine tocprint(header)
    character(len=*), intent(in), optional :: header
    if (present(header)) then
      print '(A,": ",F0.3," s")', header, toc()
    else
      print '(F0.3, " s")', toc()
    end if
  end subroutine tocprint
  
  function toc() result(elapsed)
    real :: elapsed
    call cpu_time(elapsed)
    call tic()
  end function toc
end module time

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
    integer :: n, i, clock
    call random_seed(size = n)
    allocate(seedvector(n))
    if (present(seed)) then
      seedvector = 0
      seedvector(1) = seed
    else
      call system_clock(count=clock)
      seedvector = clock + 37*[(i, i=1,n)]
    end if
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
  !     * Richard Brent, "Algorithms for minimization without derivatives",
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

module asa239

  contains
  
  ! The functions in this module are copied directly from:
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
  end function alngam
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
  end function alnorm

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
  end subroutine gamma_inc_values
  
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
    ! real ( kind = rk ) alnorm  ! Commented out when the asa239 module was introduced
    ! real ( kind = rk ) alngam  ! do
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
  end function gammad
  !-------- End of functions copied from J.Burkardt's home page -------
end module asa239

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

module minls_module
  ! The routines in this module come from pearsont3
contains
  function lstau(a,t,x,n)
    use precision
    implicit none
    integer, intent(in) :: n
    real(dp), dimension(:), intent(in) :: t
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(in) :: a
    !       Least-squares function for tau estimation.
    real(dp) :: lstau
    integer :: i
    lstau=0.0_dp
    do i=2,n
      lstau=lstau+(x(i)-x(i-1)*sign(1.0_dp,a)*abs(a)**(t(i)-t(i-1)))**2
    end do
  end function lstau
  !
  !=============================================================================
  !
  subroutine minls(n,t,x,amin,nmu_)
    use precision
    use newbrent_module
    implicit none
    integer, intent(in) :: n
    real(dp), dimension(:), intent(in) :: t
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(out) :: amin
    integer, intent(out)  :: nmu_
    !       Minimizes least-squares function s using Brent's method.
    !
    real(dp), parameter :: a_ar1=0.3678794_dp   ! 1/e
    real(dp), parameter :: tol =3.0e-08_dp      ! Brent's search, precision
    real(dp), parameter :: tol2 = 1.0e-06_dp    ! Multiple solutions, precision
    real(dp) :: dum1=-999.0_dp
    real(dp) :: dum2=-999.0_dp
    real(dp) :: dum3=-999.0_dp
    real(dp) :: dum4=-999.0_dp
    real(dp) :: a_ar11=-999.0_dp
    real(dp) :: a_ar12=-999.0_dp
    real(dp) :: a_ar13=-999.0_dp
    real(dp) :: a, bm, b, bp, c
    nmu_ = 0
    a = -2
    b = a_ar1
    bm = (b - 1)/2
    bp = (b + 1)/2
    c = 2
    dum1 = newbrent(lstau, a, b,  c, tol, a_ar11, t, x, n)
    dum2 = newbrent(lstau, b, bp, c, tol, a_ar12, t, x, n)
    dum3 = newbrent(lstau, a, bm, b, tol, a_ar13, t, x, n) 
    if ((abs(a_ar12-a_ar11) > tol2 .and. abs(a_ar12-a_ar1) > tol2) .or. &
      (abs(a_ar13-a_ar11) > tol2 .and. abs(a_ar13-a_ar1) > tol2))     &
      nmu_=1
    dum4 = min(dum1, dum2, dum3)
    if (dum4 == dum2) then
      amin = a_ar12
    else if (dum4 == dum3) then
      amin = a_ar13
    else
      amin = a_ar11
    end if
  end subroutine minls

  subroutine tauest_x(t_in, x_in, n, tau, rhoout)
    use precision
    use meanvar_module
    implicit none
    integer, intent(in) :: n        ! number of data points
    real(dp), dimension(:), intent(in) :: t_in  ! time
    real(dp), dimension(:), intent(in) :: x_in  ! time-series values
    real(dp), intent(out) :: tau    ! result: persistence time tau
    real(dp), intent(out) :: rhoout ! result: equivalent autocorrelation
    ! coefficient
    !
    !       Estimates AR(1) persistence time (tau) and
    !       equivalent autocorrelation coefficient (rho),using
    !       the least-squares algorithm of Mudelsee (2002).
    !
    !       Notes: (1) assumes t = age (see Point 1)
    !
    !              (2) automatic bias correction (see Point 5)
    !
    real(dp), dimension(n) :: x      ! x renamed
    real(dp), dimension(n) :: t      ! t renamed
    real(dp) :: avex=-999.0_dp      ! average x
    real(dp) :: varx=-999.0_dp      ! variance x
    real(dp) :: delta=-999.0_dp      ! average spacing
    real(dp) :: scalt=-999.0_dp      !  factor
    real(dp) :: rho=-999.0_dp      ! rho x, used in scaling t
    real(dp) :: rho_non=-999.0_dp      ! equivalent autocorrelation coefficient
    real(dp) :: amin=-999.0_dp      ! output from minls. Estimated value of
    ! a=exp(-scalt/tau)
    integer :: mult=-999        ! flag (multiple solution) mult=0 if 1
    ! solution,mult=1 if more than one minimum
    integer :: i=-999
    real(dp) :: rho_max = 0.99_dp
    real(dp) :: rho_min = 0.01_dp
    !
    ! 1.    Rename and change time direction
    !       ===============================
    do i=n,1,-1
      t(i)=-1.0_dp*t_in(n+1-i)
      x(i)=+1.0_dp*x_in(n+1-i)
    end do
    !
    ! 2.    Scaling (x)
    !       ==========
    call meanvar(x,avex,varx)
    x=x/sqrt(varx)
    !
    ! 3.    Scaling (t)
    !       ==========
    !       => start value of a = 1/e
    !
    delta=abs(t(n)-t(1))/(n-1)    ! average sampling
    call rhoest(n,x,rho)     ! rho estimated for x from eq 2.4 in Mudelsee(2010)
    if (rho <= 0.0_dp) then
      rho=0.05_dp
    else if (rho >= 1.0_dp) then ! this rho setting equal to 0.05 or 0.95 is not
      rho=0.95_dp        ! problematic since we require rho only for t-scaling
    end if
    scalt=-1.0_dp*log(rho)/delta    ! scaling factor
    t=t*scalt        ! t-scaled: t*(-log(a)/delta)
    !
    ! 4.    Estimation (tau)
    !       ===============
    !ls - least square function
    !minls - minimization of least-squares function ls
    !brent - Brent's search - minimum ls value
    !
    call minls(n,t,x,amin,mult)
    !
    ! 5. Result ! checked if amin (the value for a) is equal or less than zero,
    !    ====== ! or if it is equal or bigger than 1.0. That can cause
    !           ! problems in the estimation
    !
    if (mult == 1 .or. amin <= 0.0_dp .or. amin >= 1.0_dp) then
      if (amin <= 0.0_dp) then
        rho_non=rho_min
        tau=-delta/log(rho_min)
      else if (amin >= 1.0_dp) then
        rho_non=rho_max
        tau=-delta/log(rho_max)
      end if
    else
      tau = -1.0_dp/(scalt*log(amin))    ! tau - rescaled, without bias
      ! correction tau=-1/log(a)
      rho_non = exp(-delta/tau)      ! equivalent autocorrelation coefficient
      ! for tau, used in the bias correction

      ! Bias correction(unknown mean) for rho (Kendall, 1954) Solved equation
      ! (this equation: a_est =a_est_bc+(1.0_dp + 3.0_dp * a_est_bc)/(n-1)
      ! solved for estimated rho bias corrected = a_est_bc
      rho_non = (rho_non * (n-1.0_dp) + 1.0_dp) / (n-4.0_dp)

      !     print*,'rho_non',rho_non
      !     print*,'delta', delta

      if (rho_non >= 1.0_dp)then
        rho_non = exp(-delta/tau)  ! If the bias corrected equivalent
        ! autocorrelation coefficient becomes >1 then
        ! the bias correction is not performed
      end if          ! it can occur if n is small and the autocorrelation
      ! coefficent is large.
      !
      !  print*,'bias corrected rho_non',rho_non

      tau=-delta/log(rho_non)      ! tau-calculated from the bias corrected rho 
      ! a=exp(delta/tau)  gives tau=-delta/log(rho)
    end if

    rhoout = rho_non
    !
  end subroutine tauest_x

  subroutine rhoest(n,x,rho)
    use precision
    implicit none
    integer, intent(in) :: n
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(out) :: rho
    !       Estimates autocorrelation coefficient (equidistant data).
    integer i
    real(dp) :: sum1
    real(dp) :: sum2  
    sum1=0.0_dp
    sum2=0.0_dp    
    do i=2,n
      sum1=sum1+x(i)*x(i-1)
      sum2=sum2+x(i)**2
    end do
    rho=sum1/sum2
  end subroutine rhoest

end module minls_module

