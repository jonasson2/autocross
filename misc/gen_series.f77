! Spearman Data

! Program to generate test data for the Spearman program
!
! It generates:
!  o Normal bivariate AR(1) processes x and y
!    with uneven time spacing t.
!    File : normal_data.dat
!  o Lognormal bivariate AR(1) processes x_log and y_log
!    with uneven time spacing t.
!    File : lognormal_data.dat
!
! Persistence time (for x and y) and correlation coefficient
! (between x and y) are prescribed.
!
!
! Date:      May 2011
!
! Author: Kristin B. Olafsdottir
!
!==================================================================================================
!
include 'c:\nr\nr\double\nrtype.f90'
include 'c:\nr\nr\double\nr.f90'
include 'c:\nr\nr\double\nrutil.f90'
!
!==================================================================================================
!
module setting
      use nrtype
      implicit none
        save
!        
        integer,parameter  :: n = 100           !number of points
       
        real(sp) :: p =16.0_sp                  !time spacing
       
       
!
        real(sp) ::  tau_x= 2.0_sp        !persistence time (X)
        real(sp) ::  tau_y= 3.0_sp        !persistence time (Y)
      real(sp) ::  rho_epsilon= 0.60_sp   !correlation coefficient
!     
 end module setting
!==================================================================================================
!
 module data
      use nrtype
      use setting, only : n
      implicit none
      save
!     
      real(sp), dimension(n):: t, x, y, x_log, y_log
!     
 end module data
!
!==================================================================================================
!
 program random_data
      use nrtype
      implicit none
      
      call ranseed
!     
      call generate_times
      call generate
!     
 end program random_data
!
!==================================================================================================
!
 FUNCTION gamdev(ia)
        USE nrtype; USE nrutil, ONLY : assert
        USE nr, ONLY : ran
        IMPLICIT NONE
!       Numerical Recipes, modified: uses ran.
!
        INTEGER(I4B), INTENT(IN) :: ia
        REAL(SP) :: gamdev
        REAL(SP) :: am,e,h,s,x,y,v(2),arr(5)
        integer :: i
        integer :: idum=1
        call assert(ia >= 1, 'gamdev arg')
        if (ia < 6) then
                do i=1,ia
                   arr(i)=ran(idum)
                end do
                x=-log(product(arr(1:ia)))
        else
                do
                        v(1)=ran(idum)
                        v(2)=ran(idum)
                        v(2)=2.0_sp*v(2)-1.0_sp
                        if (dot_product(v,v) > 1.0) cycle
                        y=v(2)/v(1)
                        am=ia-1
                        s=sqrt(2.0_sp*am+1.0_sp)
                        x=s*y+am
                        if (x <= 0.0) cycle
                        e=(1.0_sp+y**2)*exp(am*log(x/am)-s*y)
                        h=ran(idum)
                        if (h <= e) exit
                end do
        end if
        gamdev=x
!
 END FUNCTION gamdev
!
!==================================================================================================
!
 SUBROUTINE gasdev(harvest)
        USE nrtype
        USE nr, ONLY : ran
        IMPLICIT NONE
!       Numerical Recipes, modified: uses ran.
!
        REAL(SP), INTENT(OUT) :: harvest
        REAL(SP) :: rsq,v1,v2
        REAL(SP), SAVE :: g
        LOGICAL, SAVE :: gaus_stored=.false.
        integer :: idum=1
        if (gaus_stored) then
                harvest=g
                gaus_stored=.false.
        else
                do
!                        call ran(v1)
!                        call ran(v2)
                        v1=ran(idum)
                        v2=ran(idum)
                        v1=2.0_sp*v1-1.0_sp
                        v2=2.0_sp*v2-1.0_sp
                        rsq=v1**2+v2**2
                        if (rsq > 0.0 .and. rsq < 1.0) exit
                end do
                rsq=sqrt(-2.0_sp*log(rsq)/rsq)
                harvest=v1*rsq
                g=v2*rsq
                gaus_stored=.true.
        end if
!
 END SUBROUTINE gasdev
 
 
 
!==================================================================================================
!
subroutine generate
        use nrtype
        use nr, only : ran
        use data
        use setting
        implicit none
!
!       Generates (x, y).
!
        integer :: i=0
        integer :: j=0
        integer :: idum=1
        real(sp), parameter :: zero = 0.0_sp
        real(sp), parameter :: half = 0.5_sp
        real(sp), parameter ::  one = 1.0_sp
        real(sp), parameter ::  two = 2.0_sp
        real(sp), parameter :: tiny = 1.0e-10_sp
        real(sp) :: z1=-999.0_sp
!

!
      real(sp), allocatable, dimension (:) :: ex ! epsilon_X
        real(sp), allocatable, dimension (:) :: ey ! epsilon_Y
      real(sp) :: a_x ! equivalent autocorrelation coefficient, X
      real(sp) :: a_y ! equivalent autocorrelation coefficient, Y
      real(sp) :: rho_epsilon_factor ! time-dependent white-noise correlation
!
      allocate (ex(n), ey(n))
      ex=zero
      ey=zero
!
! 1.        Bivariate white noise, (ex, ey)
!       ===============================
!
        call gasdev(z1)
        ex(1)=z1
        call gasdev(z1)
      ey(1)=rho_epsilon*ex(1)+sqrt(one-rho_epsilon**two)*z1
!
      do i=2,n
         a_x=exp(-(t(i)-t(i-1))/tau_x)
         a_y=exp(-(t(i)-t(i-1))/tau_y)
         call gasdev(z1)
         ex(i)=z1
         call gasdev(z1)
         rho_epsilon_factor=rho_epsilon*(one-a_x*a_y)/((one-a_x**two)*(one-a_y**two))**half
         ey(i)=rho_epsilon_factor*ex(i)+sqrt(one-rho_epsilon_factor**two)*z1
!
!        print *, ' subroutine generate ---      rho_epsilon        = ',rho_epsilon
!        print *, ' subroutine generate ---      rho_epsilon_factor = ',rho_epsilon_factor
!        pause
      end do
!
! 2.        Gaussian AR(1)
!       ==============
!
      x(1)=ex(1)
      y(1)=ey(1)
!
        do i=2,n
           x(i)=a_x*x(i-1)+sqrt(one-a_x**two)*ex(i)
           y(i)=a_y*y(i-1)+sqrt(one-a_y**two)*ey(i)
        end do
!
! 3.    Transform Gaussian to lognormal
!       ============================================
!
!
         do i=1,n
            x_log(i)=exp(x(i))
            y_log(i)=exp(y(i))
         end do


! 4.    Write results into two files
!       ============================================
!        
      
!     open (unit=70, file='normal_data.dat', status='unknown',    &
      open (unit=70, file='1.tmp', status='unknown',        &
            form='formatted', action='write')
         do i = 1,n
            write (70,'(f25.13,2(1x,f25.13))') t(i),x(i),y(i)
         end do
      close (unit=70, status='keep')
!       
      open (unit=60, file='lognormal_data.dat', status='unknown',       &
            form='formatted', action='write')
         do i = 1,n
            write (60,'(f25.13,2(1x,f25.13))') t(i),x_log(i),y_log(i)
         end do  
      close (unit=60, status='keep')  
!
      deallocate(ex,ey)
!
      
!
 end subroutine generate
!
!==================================================================================================
!
subroutine generate_times
        use nrtype
        use nr, only : ran
        use data, only: t
        use setting, only: n,p
        implicit none
!
!       Generates t(1:n) with average spacing = 1.
!
!       if p(1) <  0.5 : jittered: t(i) = 1.0*i + U(-p(1),+p(1)),
!                                  where U(a,b) is the uniform distribution
!                                  over the interval [a; b],
!                                  which includes equidistance (p(1) = 0.0)
!                                  coefficient of variation = p(1) * SQRT(2/3)
!
!       if p(1) >= 1.0 : gamma(INT(p(1)),
!                        coefficient of variation = SQRT(1/p(1))
!
!       if p(1) == 0.0 : equidistant
!
        integer :: i=0
        integer :: idum=1
        real(sp), parameter :: zero=0.0_sp
        real(sp), parameter :: half=0.5_sp
        real(sp), parameter ::  one=1.0_sp
        real(sp), parameter ::  two=2.0_sp
        real(sp) :: gamdev
        real(sp), dimension(n-1) :: d ! spacing
        real(sp) :: sum_d=0.0_sp
        real(sp) :: u=-999.0_sp
        d=zero
        if (p == zero) then
           do i=1,n
              t(i)=one*i
           end do
        else if (p > zero .and. p < half) then
           do i=1,n
              u=ran(idum)
              t(i)=one*i+two*(u-half)*p
           end do
        else if (p >= one) then
           do i=1,n-1
              d(i)=gamdev(int(p))
           end do
           t(1)=one
           do i=1,n-1
              t(i+1)=t(i)+d(i)
           end do
        else
           print *, ' p value not permitted - program terminates.'
           stop
        end if
!
! set average d = one
!
        do i=1,n-1
           d(i)=t(i+1)-t(i)
        end do
        sum_d=sum(d)
        d=d*(n-1)/sum_d
        t(1)=one
        do i=1,n-1
           t(i+1)=t(i)+d(i)
        end do
 end subroutine generate_times
!
!===================================================================================================
!
 subroutine ranseed
        use nrtype
        use nr, only: ran
        implicit none
!       Seeds random number generator.
        integer :: seed=0
        real(sp) :: dummy=-999.0_sp
!        print *,'seed random number generator (negative integer): '
!        read (5,*) seed
      seed=-286761062
        dummy=ran(seed)
 end subroutine ranseed
!
!====================================================================================================
!
 include 'c:\nr\nr\ran.f90'
 
 
 
 
 
 
 
 
 
 
 

Kær kveðja | Kind regards
Kristín Björg

Outlook-vl1yydtg.png

Kristín Björg Ólafsdóttir
Sérfræðingur  |  Scientist
Veðurfarsrannsóknir  |  Climate Research
Þjónustu- og rannsóknasvið  |  Service and Research Division

Veðurstofa Íslands  |  Icelandic Met Office
Bústaðavegur 7-9  |  105 Reykjavík   |  Iceland
Sími  |  Tel.  +354 522 6000







