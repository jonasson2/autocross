!
! PearsonT:
!
! Estimating Pearson's correlation coefficient with
! bootstrap confidence interval from
! serially dependent time series
!
! Author
! ======
!
! Manfred Mudelsee
! Climate Risk Analysis
! Schneiderberg 26
! 30167 Hannover
! Germany
! Email: mudelsee@mudelsee.com
! URL: http://www.mudelsee.com
!
!===================================================================================================
!
! Change log
! ==========
!
! Version  Date                Comments
!
! 1.00     July 2002            o original version
! 1.10     November 2007     	o adapted to Gnuplot 4.2
!                               o mean detrending method added
! 1.20     March 2010          	o subroutine bootstrap: 'p = 1.0 commented out
!                               o subroutine plot, point 3.1: simplified (but: has to adapted in future)
! 1.30     June 2011            o compiler switch to gfortran
!
!===================================================================================================
!
! kbo-change log (May 2013)
! =============
!
! PearsonT3:
!
! Estimating Pearson's correlation coefficient with
! calibrated bootstrap confidence interval from
! serially dependent time series!
!
! Author: 
! ========
! Kristin B. Olafsdottir (kbo)
! Climate Risk Analysis
! E-mail: olafsdottir@climate-risk-analysis.com
!
! Main changes:
!	o Bootstrap method changed from stationary bootstrap to pairwise MBB and block length selector changed
!       o Confidence intervals changed from BCa to calibrated bootstrap Student's t CI
!
!===================================================================================================
!
include 'nrtype.f90'
include 'nr.f90'
include 'nrutil.f90'
!
!===================================================================================================
!
module data1
        use nrtype
        implicit none
!       Data file (original data): n1 rows of (t1, x1, y1).
        save
        real(sp), allocatable, dimension (:) :: t1          ! time (original)
        real(sp), allocatable, dimension (:) :: x1          ! data
        real(sp), allocatable, dimension (:) :: y1          ! data
        real(sp), allocatable, dimension (:) :: t2          ! time (interval)
        real(sp), allocatable, dimension (:) :: x2          ! data
        real(sp), allocatable, dimension (:) :: y2          ! data
        real(sp), allocatable, dimension (:) :: x3          ! data (detrended)
        real(sp), allocatable, dimension (:) :: y3          ! data (detrended)
       	real(sp), allocatable, dimension (:,:) :: x3_resample1 ! resampled x3 - formed in the 1st bootstrap loop
       	real(sp), allocatable, dimension (:,:) :: y3_resample1 ! resampled y3 - formed in the 1st bootstrap loop
       	real(sp), allocatable, dimension (:) :: x3_resample2   ! resampled x3 - formed in the 2nd bootstrap loop
       	real(sp), allocatable, dimension (:) :: y3_resample2   ! resampled y3 - formed in the 2nd bootstrap loop
end module data1
!
!=====================================================================================================
!
module parameters
        use nrtype
        implicit none
!       Constants.
        save
!
! 1.    Minimum number of points
!       ========================
!
        integer, parameter :: nmin=10
!
! 2.    Number of bootstrap simulations
!       ===============================
!
        integer, parameter :: b1=2000
        integer, parameter :: b2=1000
!
! 3.    Confidence level (1 - 2 * alpha)
!       ================================
!
        real(sp) :: alpha=0.025_sp
!
! 4.    Output file names
!       =================
!
        character(len=13), parameter :: outputfile='pearsont3.dat'
!
! 5.    Some numerical parameters
!       =========================
!
        integer, parameter :: imax=10000 ! big integer
!
! 6.    Miscelleanous
!       =============
!
        integer, parameter :: ntry=5     ! user input: maximum number of errors
        
!
! 7.	Number of grid points lambda (Calibrated CI)
!	============================================
!
	integer, parameter :: n_lambda=500 	! n_lambda=500, gives lambda=0.001,0.002,,,,,,0.499,0.5
!
end module parameters
!
!==================================================================================================
!
module inv_student
	use nrtype
	implicit none
!	Inverse Student's t-distribution
	save
	real(sp) :: t_ave_alpha=-999.0_sp	! percentage point of Student's t distribution
						! with v-degrees of freedom
end module inv_student
!
!===================================================================================================
!
module inv_tlambda
	use nrtype
	use parameters, only : n_lambda
	implicit none
!	Inverse Student's t distribution over lambda range (2nd bootstrap loop)
	save
	real(sp), dimension(n_lambda) :: t_inv_lambda=-999.0_sp  ! percentage points of Student's t distribution
								 ! with v-degrees of freedom
end module inv_tlambda
!
!===================================================================================================
!
module resample_data
        use nrtype
	use parameters, only: b1,b2,n_lambda
        implicit none
!       Resampling data.
        save
	integer:: i_lambda
	real(sp), dimension(b1) :: r_resample1=-999.0_sp   		! r_xy  for resamples formed in loop1
	real(sp), dimension(b2) :: r_resample2=-999.0_sp  		! r_xy for resamples formed in loop2
	real(sp):: se_r_resample1=-999.0_sp     			! standard error for r_resample1
	real(sp),dimension(b1) :: se_r_resample2=-999.0_sp 		! standard error for r_resample2 (b1 times)
	real(sp),dimension(b1,n_lambda) :: r_low_resample1=-999.0_sp 	! b1 times lower bound (r) for the resamples1 over grid of lambda
	real(sp),dimension(b1,n_lambda) :: r_upp_resample1=-999.0_sp 	! b1 times lower bound (r) for the resamples1 over grid of lambda
!	
	real(sp),dimension(n_lambda),parameter :: 				&
	         lambda=(/ (0.5_sp*i_lambda/n_lambda, i_lambda=1,n_lambda) /)	! lambda over grid of values from 0-0.5. 
										! how tight the values are, is decided with n_lambda	
!
end module resample_data
!
!===================================================================================================
!
module result1
        use nrtype
        implicit none
        save
        real (sp) :: r=-999.0_sp        	! r_XY (Pearson's) 
        real (sp) :: r_low=-999.0_sp    	! CI lower bound 
        real (sp) :: r_upp=-999.0_sp    	! CI upper bound   
        real (sp) :: taux3=-999.0_sp    	! persistence time (x3)
        real (sp) :: tauy3=-999.0_sp    	! persistence time (y3)
	real (sp) :: rhox3=-999.0_sp		! equivalent autocorrelation coefficient for tau (x3)
	real (sp) :: rhoy3=-999.0_sp		! equivalent autocorrelation coefficient for tau (y3)
!
end module result1
!
!===================================================================================================
!
module setting
        use nrtype
        implicit none
        save
        character (len=79) :: datafile
!         character (len=1) :: dtrtype    	! detrending type ('l', linear; 'm', mean)
	character (len=1) :: dtrtype='m'	! detrending type put to 'm', mean    
        integer :: n1=-999              	! data size (original)
        integer :: n2=-999              	! data size (time interval)      
	integer :: l_mbb = -999.0_sp		! block length for MBB-moving block bootstrap   
end module setting
!
!===================================================================================================
!
module own_interfaces
!
interface
subroutine bootstrap(n,x,y,l,x_resample,y_resample)
        use nrtype
        use nr, only: ran
        implicit none
        integer, intent(in) :: n
        real(sp), dimension(n), intent(in)  :: x
        real(sp), dimension(n), intent(in)  :: y
        integer,intent(in):: l
        real(sp), dimension(n), intent(out) :: x_resample
        real(sp), dimension(n), intent(out) :: y_resample
end subroutine bootstrap
end interface
!
interface
function brent(ax,bx,cx,lstau,tol,xmin,xfunc,yfunc,nfunc)
        use nrtype; use nrutil, only : nrerror
        implicit none
        integer, intent(in) :: nfunc
        real(sp), dimension(nfunc) :: xfunc
        real(sp), dimension(nfunc) :: yfunc
        real(sp), intent(in) :: ax,bx,cx,tol
        real(sp), intent(out) :: xmin
        real(sp) :: brent
        interface
        function lstau(a,t,x,n)
               use nrtype
               implicit none
               integer, intent(in) :: n
               real(sp), dimension(n), intent(in) :: t
               real(sp), dimension(n), intent(in) :: x
               real(sp), intent(in) :: a
               real(sp) :: lstau
        end function lstau
        end interface
end function brent
end interface
!
interface
subroutine decis1(c)
        use parameters, only: ntry
        implicit none
        character(len=1), intent(out) :: c
end subroutine decis1
end interface
!
interface
subroutine fit(x,y,ndata,a,b)
        use nrtype
        use nrutil
        implicit none
        integer, intent(in) :: ndata
        real(sp), dimension(ndata), intent(in) :: x,y
        real(sp), intent(out) :: a,b
end subroutine fit
end interface
!
interface
subroutine info1(c,filec)
        use setting, only: datafile,n1,n2
        use data1, only: t1,t2
        implicit none
        character(len=1), intent(in) :: c
        character(len=12), intent(in), optional :: filec
end subroutine info1
end interface
!
interface
subroutine info2(c,filec)
        use result1, only: taux3,tauy3
        implicit none
        character(len=1), intent(in) :: c
        character(len=12), intent(in), optional :: filec
end subroutine info2
end interface
!
interface
subroutine info3(c,filec)
        use setting, only: dtrtype
        implicit none
        character(len=1), intent(in) :: c
        character(len=12), intent(in), optional :: filec
end subroutine info3
end interface
!
interface
subroutine info4(c,filec)
        use nrtype
        use data1, only: t2,x2,y2,x3,y3
        use parameters, only: alpha
        use result1, only: r,r_low,r_upp
        use setting, only: n2
        implicit none
        character(len=1), intent(in) :: c
        character(len=12), intent(in), optional :: filec
end subroutine info4
end interface
!
interface
subroutine interval(c)
        use nrtype
        use data1, only: t1,x1,y1,t2,x2,y2,x3,y3,x3_resample1,y3_resample1,	&
                         x3_resample2,y3_resample2
        use parameters, only: nmin,ntry
        use setting, only: n1,n2
        implicit none
        character(len=1), intent(in) :: c
end subroutine interval
end interface
!
interface
function invtanh(x)
	use nrtype
	implicit none
	real(sp),intent(in) :: x
	real(sp) :: invtanh
end function
end interface
!
interface
function lstau(a,t,x,n)
        use nrtype
        implicit none
        integer, intent(in) :: n
        real(sp), dimension(n), intent(in) :: t
        real(sp), dimension(n), intent(in) :: x
        real(sp), intent(in) :: a
        real(sp) :: lstau
end function lstau
end interface
!
interface
 subroutine minls(n,t,x,amin,nmu_)
        use nrtype
        implicit none
        integer, intent(in) :: n
        real(sp), dimension(n), intent(in) :: t
        real(sp), dimension(n), intent(in) :: x
        real(sp), intent(out) :: amin
        integer, intent(out)  :: nmu_
        interface
        function lstau(a,t,x,n)
               use nrtype
               implicit none
               integer, intent(in) :: n
               real(sp), dimension(n), intent(in) :: t
               real(sp), dimension(n), intent(in) :: x
               real(sp), intent(in) :: a
               real(sp) :: lstau
        end function lstau
        end interface
end subroutine minls
end interface
!
interface
subroutine pearsn(x,y,n,r)
        use nrtype
        use nrutil
        implicit none
        integer, intent(in) :: n
        real(sp), dimension(n), intent(in) :: x,y
        real(sp), intent(out) :: r
end subroutine pearsn
end interface
!
interface
function phi(x)
        use nrtype
        use nr, only: erfcc
        implicit none
        real(sp) :: phi
        real(sp), intent(in) :: x
end function phi
end interface
!
interface
function phi_inv(prob)
        use nrtype
        implicit none
        real(sp) :: phi_inv
        real(sp), intent(in) :: prob
end function phi_inv
end interface
!
interface
subroutine plot(type)
        use nrtype
        use data1, only: t1,t2,x1,x2,x3,y1,y2,y3
        use result1, only: r,r_low,r_upp
        use setting, only: n1,n2
        implicit none
        integer, intent(in) :: type
end subroutine plot
end interface
!
interface
subroutine rhoest(n,x,rho)
        use nrtype
        implicit none
        integer, intent(in) :: n
        real(sp), dimension(n), intent(in) :: x
        real(sp), intent(out) :: rho
end subroutine rhoest
end interface
!
interface
 subroutine tauest_x(t_in,x_in,n,c,tau,rhoout)
        use nrtype
        use nr, only: avevar
        implicit none
        character (len=1), intent(in) :: c
        integer, intent(in) :: n
        real(sp), dimension(n), intent(in) :: t_in
        real(sp), dimension(n), intent(in) :: x_in
        real(sp), intent(out) :: tau
	real(sp), intent(out) :: rhoout
end subroutine tauest_x
end interface
!
interface
function t_inv(alpha,dof)
	use nrtype
	use inv_student
	implicit none
	real(sp),intent(in) :: alpha
	integer,intent(in) :: dof
	real(sp) :: t_inv
end function t_inv
end interface
!
end module own_interfaces
!
!===================================================================================================
!
program pear
	use own_interfaces
        use parameters, only: imax
        implicit none
        character(len=1) :: c1
        integer :: i1
!
! 1.    Welcome
!       =======
        call info0      		! welcome message
        call ranseed
!
! 2.    Data
!       =====
        call chsett1   		! changes setting: data file name
        call chsett2    		! changes setting: n1
        call allocate0  		! t1, x1, y1
        call init0        		! t1, x1, y1
        call read1      		! reads data
!
! 3.    Time interval extraction and calculation
!       ======================================
        call init1a     			! n2=n1
	call calc_inv_student  	! calculates percentage point tv(alpha) of the Student's t distribution (after alpha and n have been defined)
	call calc_t_inv_lambda	! calculates percentage point tv(lambda) over a lambda grid (Calibrated CI)      
        call allocate1  			! t2, x2, y2, x3, y3, x3_resample1, y3_resample1, x3_resample2, y3_resample2
        call init1b     			! t2, x2, y2, x3, y3, x3_resample1, y3_resample1, x3_resample2, y3_resample2
outer4: do i1=1,imax
!                call chsett3           ! changes setting: dtrtype
              call r_est               	! detrends (x2->x3, y2->y3)
					! estimates r(x3, y3)
              call tauest             	! estimates persistence times taux3, tauy3 and rhox3 and rhoy3)
              call chsett4            	! changes setting: l_mbb , block length
              call confidence        	! estimates [r_low; r_upp]
              call plot(1)             	! t2, x2, y2
              call plot(2)             	! x2, y2
              call info1('p')          	! extracted time interval
              call info2('p')          	! persistence times
              call info3('p')          	! detrending method 
              call info4('p')          	! r [ r_low; r_upp ]
              call decis1(c1)          	! decision tree: Part 1
           if (c1 == 'o') then      	! original time interval
              call interval(c1)
           else if (c1 == 'n') then ! new time interval
              call interval(c1)
           else if (c1 == 'x') then
              exit outer4
           end if
        end do outer4
!
! 4.    Output and exit
!       ===============
        call output
        call deallocate1
end program pear
!
!===================================================================================================
!
subroutine allocate0
        use nrtype
        use setting, only: n1
        use data1, only: t1,x1,y1
        implicit none
!       Allocates t1, x1, y1.
        integer :: error=0
        allocate(t1(n1), x1(n1), y1(n1), stat=error)
        if (error /= 0) then
           print *,'Subroutine allocate0: space requested not possible - PearsonT terminates.'
           stop
        end if
end subroutine allocate0
!
!===================================================================================================
!
subroutine allocate1
        use setting, only: n2
	use data1, only: t2,x2,y2,x3,y3,x3_resample1,y3_resample1,      &
	                 x3_resample2,y3_resample2
	use parameters, only: b1
        implicit none
!       Allocates t2, x2, y2, x3, y3, x3_resample1,y3_resample1, x3_resample2, y3_resample2
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
!===================================================================================================
!
subroutine bootstrap(n,x,y,l,x_resample,y_resample)
        use nrtype
        use nr, only: ran
        implicit none
        integer, intent(in) :: n			    ! number of data points
        real(sp), dimension(n), intent(in)  :: x	    ! original x time series
        real(sp), dimension(n), intent(in)  :: y	    ! original y time series
        integer, intent(in) :: l			    ! block length
        real(sp), dimension(n), intent(out) :: x_resample   ! bootstrap resampled x time series	
        real(sp), dimension(n), intent(out) :: y_resample   ! bootstrap resampled y time series
!
!       Pairwise Moving Block bootstrap resampling
!	Algorithm 7.2 and 3.1 in (Mudelsee, 2010)
!	Block length = l = l_mbb calculated in subroutine chsett4.
!	Blocks are randomly selected and concatenated until n data
!	are resampled (resize last block by removing last observation in it).
!	By using the same random bootstrap index for x_resample and y_resample,
!	(x(i),y(i)) pairs are resampled. 
!
        integer, dimension(n) :: indxx
        integer :: i=0
        integer :: j=0
        integer :: k=0
        integer :: idum=1
        integer :: n_block_start
        n_block_start=n-l+1
        if (l == 1) then 				! block length less or equal to average spacing, ordinary bootstrap is used
           do i=1,n
              indxx(i)=int(n*ran(idum))+1
           end do
        else
           k=1						! set counter
outer:     do i=1,n
              indxx(k)=int(n_block_start*ran(idum))+1	! draw random block start
              k=k+1
              if (k > n) exit outer
              do j=1,l-1				! fill up the block of length l
                 indxx(k)=indxx(k-1)+1
                 k=k+1
                 if (k > n) exit outer
              end do					! the loop keeps on going until k>n or k==n
           end do outer
        end if
        x_resample(:)=x(indxx(:))
        y_resample(:)=y(indxx(:))
!
end subroutine bootstrap
!
!===================================================================================================
!
function brent(ax,bx,cx,lstau,tol,xmin,xfunc,yfunc,nfunc)
        use nrtype; use nrutil, only : nrerror
        implicit none
!       brents minimization (numerical recipes, modified: last three arguments).
        integer, intent(in) :: nfunc
        real(sp), dimension(nfunc) :: xfunc
        real(sp), dimension(nfunc) :: yfunc
        real(sp), intent(in) :: ax,bx,cx,tol
        real(sp), intent(out) :: xmin
        real(sp) :: brent
        interface
        function lstau(a,t,x,n)
               use nrtype
               implicit none
               integer, intent(in) :: n
               real(sp), dimension(n), intent(in) :: t
               real(sp), dimension(n), intent(in) :: x
               real(sp), intent(in) :: a
               real(sp) :: lstau
        end function lstau
        end interface
        integer(i4b), parameter :: itmax=100
        real(sp), parameter :: cgold=0.3819660_sp,zeps=1.0e-3_sp*epsilon(ax)
        integer(i4b) :: iter
        real(sp) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
        a=min(ax,cx)
        b=max(ax,cx)
        v=bx
        w=v
        x=v
        e=0.0
        fx=lstau(x,xfunc,yfunc,nfunc)
        fv=fx
        fw=fx
        do iter=1,itmax
                xm=0.5_sp*(a+b)
                tol1=tol*abs(x)+zeps
                tol2=2.0_sp*tol1
                if (abs(x-xm) <= (tol2-0.5_sp*(b-a))) then
                        xmin=x
                        brent=fx
                        return
                end if
                if (abs(e) > tol1) then
                        r=(x-w)*(fx-fv)
                        q=(x-v)*(fx-fw)
                        p=(x-v)*q-(x-w)*r
                        q=2.0_sp*(q-r)
                        if (q > 0.0) p=-p
                        q=abs(q)
                        etemp=e
                        e=d
                        if (abs(p) >= abs(0.5_sp*q*etemp) .or. &
                                p <= q*(a-x) .or. p >= q*(b-x)) then
                                e=merge(a-x,b-x, x >= xm )
                                d=cgold*e
                        else
                                d=p/q
                                u=x+d
                                if (u-a < tol2 .or. b-u < tol2) d=sign(tol1,xm-x)
                        end if
                else
                        e=merge(a-x,b-x, x >= xm )
                        d=cgold*e
                end if
                u=merge(x+d,x+sign(tol1,d), abs(d) >= tol1 )
                fu=lstau(u,xfunc,yfunc,nfunc)
                if (fu <= fx) then
                        if (u >= x) then
                                a=x
                        else
                                b=x
                        end if
                        call shft(v,w,x,u)
                        call shft(fv,fw,fx,fu)
                else
                        if (u < x) then
                                a=u
                        else
                                b=u
                        end if
                        if (fu <= fw .or. w == x) then
                                v=w
                                fv=fw
                                w=u
                                fw=fu
                        else if (fu <= fv .or. v == x .or. v == w) then
                                v=u
                                fv=fu
                        end if
                end if
        end do
        call nrerror('brent: exceed maximum iterations')
        contains
!
        subroutine shft(a,b,c,d)
        real(sp), intent(out) :: a
        real(sp), intent(inout) :: b,c
        real(sp), intent(in) :: d
        a=b
        b=c
        c=d
        end subroutine shft
end function brent
!
!==================================================================================================
!
subroutine calc_inv_student
        use nrtype
        use inv_student
        use parameters, only: alpha
        use setting, only: n2
        implicit none
!       Calculates inverse of t(alpha, dof),
!       dof = n - 2.	degrees of freedom
        real(sp) :: t_inv
!
        t_ave_alpha = t_inv(alpha,n2-2)
!        
 end subroutine calc_inv_student
!
!===================================================================================================
!
subroutine calc_t_inv_lambda
	use nrtype
	use inv_tlambda
	use parameters, only: n_lambda
	use setting, only: n2
	use resample_data, only: lambda
!       Calculates inverse of t(alpha, dof).
!       dof = n -1 (mean estimation) (von Storch and Zwiers 1999: p. 92)
!
! 	calculates for lambda grid; dof = 2n-5 (x, y; E_x, E_y, S_x,S_y, r_xy estimation).
!
	real(sp):: t_inv
	integer:: i
!
	do i=1,n_lambda
	   t_inv_lambda(i)=t_inv(lambda(i),2*n2-5)
	end do
!		
end subroutine calc_t_inv_lambda
!
!===================================================================================================
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
           datafile = 'test_data.txt'
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
!===================================================================================================
!
subroutine chsett2
        use nrtype
        use parameters, only: nmin
        use setting, only: datafile,n1
        implicit none
!       Determines n, see REDFIT35.F90 subroutine setdim
!       (Schulz and Mudelsee 2002).
        character (len = 1) :: flag
        integer :: i,iocheck,open_error
        real(sp) :: tdum, xdum, ydum  ! dummy
        open (unit=1, file=datafile, status='old',                      &
              form='formatted', action='read', iostat=open_error)
        if (open_error /= 0 ) then
           print *,'Error during data file opening - PearsonT terminates.'
           stop
        end if
!
! 1.    Skip header
!       ===========
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
!       ==========
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
!===================================================================================================
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
!===================================================================================================
!
subroutine chsett4
 	use nrtype
 	use data1, only: t2
 	use result1, only: rhox3,rhoy3 	
 	use setting, only: n2, l_mbb
 	implicit none
! 	
!	Changes setting: 
!	Calculates block length (formula from Carlstein et al, 1986. Sherman et al., 1998 adapted it to MBB)
!	Equations 7.31 and 7.32 in Mudelsee, 2010  	
!
	real(sp) :: axy			
	real(sp), parameter :: onethird=0.3333333333333333333333333333333333333333_sp
	real(sp), parameter :: twothird=0.6666666666666666666666666666666666666666_sp
	real(sp), parameter :: half=0.5_sp
	real(sp), parameter :: one=1.0_sp
	real(sp), parameter :: two=2.0_sp
!
	axy = sqrt(rhox3*rhoy3) !Eq 7.31 in (Mudelsee, 2010)
!	
	l_mbb = anint(((6.0_sp**half*axy)/(one-axy**two))**twothird*n2**onethird)        !Eq 7.32 in (Mudelsee,2010)
!	
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
!===================================================================================================
!
subroutine confidence
	use nrtype
	use nr, only: avevar
	use data1, only: x3,x3_resample1,x3_resample2,                  &
	                 y3,y3_resample1,y3_resample2
	use inv_student
	use inv_tlambda
	use parameters, only: alpha,b1,b2,n_lambda
	use resample_data
	use result1, only: r,r_low,r_upp
	use setting, only: n2,l_mbb
	implicit none
!
!       Calibrated Student's t confidence interval for r_xy 
!
	integer :: j=0			! for j=1,b1
	integer :: k=0			! for k=1,b2
	integer :: l=0			! for l=1,n_lambda
!	
	real(sp) :: ave_se_1=-999.0_sp
	real(sp) :: var_se_1=-999.0_sp
	real(sp) :: ave_se_2=-999.0_sp
	real(sp) :: var_se_2=-999.0_sp
!	
	real(sp),dimension(n_lambda)::plambda=-999.0_sp
	real(sp):: alphatest=-999.0_sp				! test alpha (1 - 2 * alpha)
	integer :: l_loc=0 					! gives the location of lambda value, within the lambda,						
!								! which gives the coverage closest to 1-2*alpha  array
	real(sp) :: invtanh
	character(len=1),parameter::flag_fisher='y'		! 'y' = Fisher's transform
!								! 'n' = without Fisher's transform
!
! 	1. First bootstrap loop (Student's t confidence interval)
! 	=========================================================
! 
!	1.1 Bootstrap
!	=============
!	Form bootstrap resamples x3_resample1 and y3_resample1 and store it
!	Estimate Pearson's correlation coefficient, r_resample1 and store it
!
	do j=1,b1
	   call bootstrap(n2,x3,y3,l_mbb,x3_resample1(j,1:n2),y3_resample1(j,1:n2))
	   call pearsn(x3_resample1(j,1:n2),y3_resample1(j,1:n2),n2,r_resample1(j))
	end do
    print *, 'n2,b1=', n2, b1
    print *, 'x3=', x3_resample1(1,1:5), x3_resample1(1,n2-4:n2)
    print *, 'x3=', x3_resample1(b1,1:5), x3_resample1(b1,n2-4:n2)
    print *, 'y3=', y3_resample1(1,1:5), y3_resample1(1,n2-4:n2)
    print *, 'y3=', y3_resample1(b1,1:5), y3_resample1(b1,n2-4:n2)
    print *, 'r_=', r_resample1(1:5), r_resample1(b1-4:b1)
!
!	1.2 Fisher's z-transformation
!	=============================
!	z-transform the r_resamples z=invtanh(r_resample1)
!
	if (flag_fisher == 'y')then
	      do j=1,b1
	         r_resample1(j)=invtanh(r_resample1(j))
	      end do
	  
	end if
!
!	1.3 Bootstrap standard error
!	============================
!	Estimate the bootstrap se from all the replications
!	
	call avevar(r_resample1(1:b1),ave_se_1,var_se_1)
	se_r_resample1 = sqrt(var_se_1)	
!

!	1.4 Bootstrap student's t confidence interval			
!	=============================================
!	Calculate CI for r and re-transform z
!
	
!	if (flag_fisher == 'y')then	   
!	      r_low= tanh(invtanh(r)+ t_ave_alpha * se_r_resample1)
!	      r_upp = tanh(invtanh(r)- t_ave_alpha * se_r_resample1)
!	      print*,'r_low=',r_low
!	      print*,'r_upp=',r_upp	   
!	else if (flag_fisher == 'n')then
!	      r_low = r + t_ave_alpha * se_r_resample1
!	      r_upp = r- t_ave_alpha * se_r_resample1
!	      print*,'r_low=',r_low
!	      print*,'r_upp=',r_upp
!	end if

	print*, '1st bootstrap loop finished '
!
! 2. Second bootstrap loop (Forms CI for each resample1)
! =========================================================
!
!	2.1 Big bootstrap loop
!	======================
!	Goes b1=2000 times through the 2nd bootstrap loop
        
!	
boot:	do j=1,b1

	   If (j ==100) print*, '2nd bootstrap loop working...' 
	   If (mod(j,500) == 0) print*, '2nd bootstrap loop working...' ,j,'/',b1
	
!
!	2.2 2nd Bootstrap
!	=================
!	Forms bootstrap resamples x3_resample2 and y3_resample2 from the resamples
!	from bootstrap loop 1. Forms b1*b2 resamples2 (no need to store them between loops).
!	The block length is overtaken from the original samples x3 and y3.
!	Estimates Pearson's correlation coefficient r_resample2 and stores it b2 times (just within this loop,
!	not within the big loop).
!
	   do k=1,b2
	      call bootstrap(n2,x3_resample1(j,1:n2),y3_resample1(j,1:n2),l_mbb,x3_resample2(1:n2),y3_resample2(1:n2))
	      call pearsn(x3_resample2(1:n2),y3_resample2(1:n2),n2,r_resample2(k))
	   end do
!
!	2.3 Fisher's z-transform
!	========================
!
	   if (flag_fisher == 'y')then
	      do k=1,b2
	            r_resample2(k) = invtanh(r_resample2(k))
	      end do
	   end if
!
!	2.4 Confidence intervals for resamples 1
!	========================================
!
!	2.4.1 Bootstrap standard error
!	==============================
!	Estimates the bootstrap standard error from the replications, r_resample2.
!	Results with b1 times se, one for each big loop or one for each of the 2000 CI for resample1
!	
	   call avevar(r_resample2(1:b2),ave_se_2,var_se_2)
	   se_r_resample2(j) = sqrt(var_se_2)  	
!	
!
!	2.4.2 Student's t CI for each resample1 (b1 times r_upp and r_low)
!	==================================================================
!	r_upp and r_low estimated for each bootstrap resample1 (b1 times)
!	over a grid of confidence levels (lambda)
!	
	   if (flag_fisher == 'y')then
	         do l=1,n_lambda
	            r_low_resample1(j,l) = tanh(r_resample1(j)+t_inv_lambda(l)*se_r_resample2(j))
	            r_upp_resample1(j,l) = tanh(r_resample1(j)-t_inv_lambda(l)*se_r_resample2(j))
	         end do
	   else if (flag_fisher == 'n')then
	         do l=1,n_lambda
	            r_low_resample1(j,l) = r_resample1(j)+t_inv_lambda(l)*se_r_resample2(j)
	            r_upp_resample1(j,l) = r_resample1(j)-t_inv_lambda(l)*se_r_resample2(j)
	         end do
	   end if
!		
!
!	2.5 End of big Bootstrap loop
!	=============================
	end do boot
	
!	
!
! 3. Determination of two-sides plambda
! =====================================
! 	Equation 3.47 in Mudelsee(2010)
!
	do l = 1,n_lambda
	   plambda(l)= 1.0_sp *count(r .ge. r_low_resample1(1:b1,l) .and.	&
	                                         r .le. r_upp_resample1(1:b1,l))/b1
	end do
!	
! 4. Calibration - calibrated Student's t confidence interval
! ===========================================================
!	
!	4.1 Find new confidence points
!	==============================
!	Find the lambda value which gives coverage closest to the nominal value (1-2*alpha)
!		
	alphatest=1.0_sp-(2.0_sp*alpha)
!	 
	l_loc = maxval(minloc(abs(alphatest-plambda(:))))
!
!	4.2 Calculate calibrated Student's t CI
!	=======================================
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
!===================================================================================================
!
subroutine deallocate1
        use data1, only: t2,x2,y2,x3,y3,x3_resample1,y3_resample1,	&
                         x3_resample2,y3_resample2
        implicit none
!       Deallocates t2, x2, y2, x3, y3, x3_resample1, y3_resample1, x3_resample2, y3_remsample2
        integer :: error=0
!	
        deallocate(t2, x2, y2, x3, y3, x3_resample1, y3_resample1,      &
                   x3_resample2,y3_resample2, stat=error)
!
        if (error /= 0) then
           print *,'Subroutine deallocate: Unexpected deallocation error - PearsonT terminates.'
           stop
        end if
end subroutine deallocate1
!
!===================================================================================================
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
!===================================================================================================
!
subroutine fit(x,y,ndata,a,b)
        use nrtype
        use nrutil
        implicit none
!       Straight line fit (Numerical Recipes, modified).
        integer, intent(in) :: ndata
        real(sp), dimension(ndata), intent(in) :: x,y
        real(sp), intent(out) :: a,b
        real(sp) :: ss,sx,sxoss,sy,st2
        real(sp), dimension(ndata) :: t
        ss=real(ndata,sp)
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
!===================================================================================================
!
subroutine info0
        implicit none
!       Welcomes.
        integer :: i
        do i=1,10
           print '(a)',' '
        end do
        print '(a)',     '==============================================================================='
        print '(a)',     '                                                                               			'
        print '(a,a,a)', 'PearsonT3 (Version 1.0',',',' May 2013)                                   		'
        print '(a)',     '                                                                               			'
        print '(a)',     '==============================================================================='
        print '(a)',     '                                                                               			'
        print '(a)',     '  Reference:                                                                   		'
        print '(a)',     '                                                                               			'
        print '(a)',     '  Olafsdottir and Mudelsee. More accurate, calibrated bootstrap 	'
        print '(a)',     '                            confidence intervals for correlating two          	'
        print '(a)',     '                            serially dependent time series.                      	'
        print '(a)',     '                            (submitted to Mathematical Geosciences)                             	'
        print '(a)',     '                                                                               			'
        print '(a)',     '                                                                               			'
        print '(a)',     '  Download:                                                                    		'
        print '(a)',     '                                                                               			'
        print '(a)',     '  http://www.climate-risk-analysis.com                                                      	'
        print '(a)',     '                                                                               			'
        print '(a)',     '==============================================================================='
        print '(a)',     '                                                                               			'
end subroutine info0
!
!==================================================================================================
!
subroutine info1(c,filec)
        use setting, only: datafile,n1,n2
        use data1, only: t1,t2
        implicit none
        character(len=1), intent(in) :: c
        character(len=13), intent(in), optional :: filec
!       Informs about: filename, datatype, number of points:
!       observed and extracted time interval.
        integer :: open_error
        if (c == 'p') then
           print '(a,a)',                                              			 	&
           'Data file name:              ',datafile
           print '(a,f11.3,a,f11.3,a,i8,a)',                            			&
           'Time interval - original:    [ ',t1(1),'; ',t1(n1),' ]     - ',n1,' points'
           if (n1 /= n2)                                                			&
              print '(a,f11.3,a,f11.3,a,i8,a)',                         			&
              'Time interval - extracted:   [ ',t2(1),'; ',t2(n2),' ]     - ',n2,' points'
        else if (c == 'w') then
           open (unit=1, file=filec, status='replace',                   			&
                 form='formatted', action='write', iostat=open_error)
           if (open_error /= 0 ) then
              print *,'Error during data file opening - terminate'
              stop
           end if
           write (unit=1, fmt='(a)')                                    			&
           'PearsonT3 (www.climate-risk-analysis.com):  Output'
           write (unit=1, fmt='(a)')                                    			&
           ' '
           write (unit=1, fmt='(a,a)')                                  			&
           'Data file name:              ',datafile
           write (unit=1, fmt='(a,f11.3,a,f11.3,a,i8,a)')               	&
           'Time interval - original:    [ ',t1(1),'; ',t1(n1),' ]     - ',n1,' points'
           if (n1 /= n2)                                                			&
              write (unit=1, fmt='(a,f11.3,a,f11.3,a,i8,a)')            	&
             'Time interval - extracted:   [ ',t2(1),'; ',t2(n2),' ]     - ',n2,' points'
           close (unit=1, status='keep')
        end if
end subroutine info1
!
!==================================================================================================
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
           print '(a,23x,f11.3)',                                       		&
           'Persistence time (x):',taux3
           print '(a,23x,f11.3)',                                       		&
           'Persistence time (y):',tauy3
        else if (c == 'w') then
           open (unit=1, file=filec, status='old', position='append',    		&
                 form='formatted', action='write', iostat=open_error)
           if (open_error /= 0 ) then
              print *,'Error during data file opening - PearsonT terminates.'
              stop
           end if
           write (unit=1, fmt='(a)')                                    		&
           ' '
           write (unit=1, fmt='(a,23x,f11.3)')                          		&
           'Persistence time (x):',taux3
           write (unit=1, fmt='(a,23x,f11.3)')                          		&
           'Persistence time (y):',tauy3
           close (unit=1, status='keep')
        end if
end subroutine info2
!
!==================================================================================================
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
              print '(a,29x,a)',                                        		&
              'Detrending type:     ','linear'
           else if (dtrtype .eq. 'm') then
              print '(a,29x,a)',                                        		&
              'Detrending type:     ','mean'
           end if
        else if (c == 'w') then
           open (unit=1, file=filec, status='old', position='append',  			&
                 form='formatted', action='write', iostat=open_error)
           if (open_error /= 0 ) then
              print *,'Error during data file opening - PearsonT terminates.'
              stop
           end if
           write (unit=1, fmt='(a)')                                    		&
           ' '
           if (dtrtype .eq. 'l') then
              write (unit=1, fmt='(a,29x,a)')                           		&
              'Detrending type:     ','linear'
           else if (dtrtype .eq. 'm') then
              write (unit=1, fmt='(a,29x,a)')                           		&
              'Detrending type:     ','mean'
           end if
           close (unit=1, status='keep')
        end if
end subroutine info3
!
!==================================================================================================
!
subroutine info4(c,filec)
        use nrtype
        use data1, only: t2,x2,y2,x3,y3
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
           print '(a,a,i2,a,11x,4(f6.3,a))',                            		&
                 'Pearson''s r',', ',            					&
                 nint(100*(1-2.0_sp*alpha)),' % confidence interval:',   		&
                 r, '  [ ',r_low,'; ',r_upp,' ]'
        else if (c == 'w') then
	   open (unit=1, file=filec, status='old', position='append',  			&
                    form='formatted', action='write', iostat=open_error)
           if (open_error /= 0 ) then
	      print *,'Error during data file opening - PearsonT terminates.'
              stop
           end if
           write (unit=1, fmt='(a)') ' '
           write (unit=1, fmt='(a,a,i2,a,11x,4(f6.3,a))')              	 		&
                 'Pearson''s r',', ',                                   		&
                 nint(100*(1-2.0_sp*alpha)),' % confidence interval:',   	 	&
                  r, '  [ ',r_low,'; ',r_upp,' ]'
           write (unit=1, fmt='(a)') ' '
           write (unit=1,fmt='(2x,7(a))')                               		&
                 ' time          ',                                     		&
                 ' x             ',                                     		&
                 ' y             ',                                    	 		&
                 ' trend(x)      ',                                     		&
                 ' trend(y)      ',                                     		&
                 ' x-detrended   ',                                     		&
                 ' y-detrended   '
           do i=1,n2
              write (unit=1,fmt='(7(1x,es14.6))') t2(i),x2(i),y2(i),    		&
                                                  x2(i)-x3(i),          		&
                                                  y2(i)-y3(i),          		&
                                                  x3(i),y3(i)
           end do
           close (unit=1, status='keep')
        end if
end subroutine info4
!
!===================================================================================================
!
subroutine init0
        use nrtype
        use data1, only: t1, x1, y1
        implicit none
!          Initializes t1, x1, y1.
        t1=-999.0_sp
        x1=-999.0_sp
        y1=-999.0_sp
end subroutine init0
!
!===================================================================================================
!
subroutine init1a
        use setting, only: n1,n2
        implicit none
!         Initializes: n2=n1.
        n2=n1
end subroutine init1a
!
!===================================================================================================
!
subroutine init1b
        use nrtype
        use data1, only: t1,t2,x1,x2,x3,y1,y2,y3,x3_resample1,y3_resample1,	&
                         x3_resample2,y3_resample2    
        implicit none
!       Initializes t2, x2, y2, x3, y3, x3_resample1, y3_resample1, x3_resample2, y3_resample2
        t2=t1
        x2=x1
        y2=y1
        x3=x2
        y3=y2
        x3_resample1=-999.0_sp
        y3_resample1=-999.0_sp
        x3_resample2=-999.0_sp
        y3_resample2=-999.0_sp
!
end subroutine init1b
!
!===================================================================================================
!
subroutine interval(c)
        use nrtype
        use data1, only: t1,x1,y1,t2,x2,y2,x3,y3,x3_resample1,y3_resample1,	&
                         x3_resample2,y3_resample2
        use parameters, only: nmin,ntry
        use setting, only: n1,n2
        implicit none
        character(len=1), intent(in) :: c
!       Defines new: n2, t2, x2, y2.
        integer :: i1,i2,j,k1
        real(sp) :: tl,tr
        if (c == 'n') then
           do i1=1,ntry
              do i2=1,ntry
                 print '(a,i6,a)','New interval must contain at least',nmin,' points'
                 print '(3(a))',  '[left boundary',',',' right boundary]: '
                 read (5,*) tl, tr
                 tl=max(tl,t1(1))
                 tr=min(tr,t1(n1))
                 if (tl < t1(n1) .and. tr > t1(1) .and. tl < tr) exit
              end do
              if (i2 > ntry) then
                 print *,'OK - PearsonT terminates.'
                 stop
              end if
              if (tl == t1(1)) j=1
              if (tr == t1(n1)) k1=n1
              do i2=2,n1
                 if (tl >  t1(i2-1) .and. tl <= t1(i2)) j=i2
                 if (tr >= t1(i2-1) .and. tr <  t1(i2)) k1=i2-1
              end do
              n2=k1-j+1
              if (n2 >= nmin) exit
              print '(a,i6)','Time interval contains too few points: ',n2
           end do
           if (i1 > ntry) then
              print *,'OK - PearsonT terminates.'
              stop
           end if
           call deallocate1	! t2, x2, y2, x3, y3, x3_resample1, y3_resample1, x3_resample2, y3_resample2
           call allocate1   	! t2, x2, y2, x3, y3, x3_resample1, y3_resample1, x3_resample2, y3_resample2
           do i2=j,k1
              t2(i2-j+1)=t1(i2)
              x2(i2-j+1)=x1(i2)
              y2(i2-j+1)=y1(i2)
           end do
           x3=x2
           y3=y2
           x3_resample1=-999.0_sp
           y3_resample1=-999.0_sp
           x3_resample2=-999.0_sp
           y3_resample2=-999.0_sp
        else if (c == 'o') then
           call init1a      	! n2=n1
           call deallocate1	! t2, x2, y2, x3, y3, x3_resample1, y3_resample1, x3_resample2, y3_resample2
           call allocate1   	! t2, x2, y2, x3, y3, x3_resample1, y3_resample1, x3_resample2, y3_resample2
           call init1b      	! t2, x2, y2, x3, y3, x3_resample1, y3_resample1, x3_resample2, y3_resample2
        end if
end subroutine interval
!
!===================================================================================================
!
 function invtanh(x)
        use nrtype
        implicit none
        real(sp) :: invtanh
!       Inverse hyperbolic tangens.
!
        real(sp),  intent(in) :: x
	real(sp), parameter :: tinyx=1.0e-13
	real(sp), parameter :: half=0.5_sp
	real(sp), parameter :: one=1.0_sp
	real(sp), parameter :: big=16.0_sp
	if (x.ge.one-tinyx) then
	   invtanh=big
	else if (x.le.-one+tinyx)then
	   invtanh=-big	
	else
	   invtanh=half*log((one+x)/(one-x))
	end if
 end function invtanh
!
!===================================================================================================
!
function lstau(a,t,x,n)
        use nrtype
        implicit none
        integer, intent(in) :: n
        real(sp), dimension(n), intent(in) :: t
        real(sp), dimension(n), intent(in) :: x
        real(sp), intent(in) :: a
!       Least-squares function for tau estimation.
        real(sp) :: lstau
        integer :: i
        lstau=0.0_sp
        do i=2,n
           lstau=lstau+(x(i)-x(i-1)*sign(1.0_sp,a)*abs(a)**(t(i)-t(i-1)))**2.0_sp
        end do
end function lstau
!
!===================================================================================================
!
 subroutine minls(n,t,x,amin,nmu_)
        use nrtype
        implicit none
        integer, intent(in) :: n
        real(sp), dimension(n), intent(in) :: t
        real(sp), dimension(n), intent(in) :: x
        real(sp), intent(out) :: amin
        integer, intent(out)  :: nmu_
!       Minimizes least-squares function s using Brent's method.
!
        real(sp), parameter :: a_ar1=0.3678794_sp	! 1/e
        real(sp), parameter :: tol =3.0e-08_sp		!Brent's search, precision
        real(sp), parameter :: tol2=1.0e-06_sp		!Multiple solutions, precision
        real(sp) :: dum1=-999.0_sp
        real(sp) :: dum2=-999.0_sp
        real(sp) :: dum3=-999.0_sp
        real(sp) :: dum4=-999.0_sp
        real(sp) :: a_ar11=-999.0_sp
        real(sp) :: a_ar12=-999.0_sp
        real(sp) :: a_ar13=-999.0_sp
        real(sp) :: brent
        interface
        function lstau(a,t,x,n)
               use nrtype
               implicit none
               integer, intent(in) :: n
               real(sp), dimension(n), intent(in) :: t
               real(sp), dimension(n), intent(in) :: x
               real(sp), intent(in) :: a
               real(sp) :: lstau
        end function lstau
        end interface
        nmu_=0
        dum1=brent(-2.0_sp,                a_ar1, +2.0_sp, lstau, tol, a_ar11, t,x,n)
        dum2=brent( a_ar1, 0.5_sp*(a_ar1+1.0_sp), +2.0_sp, lstau, tol, a_ar12, t,x,n)
        dum3=brent(-2.0_sp, 0.5_sp*(a_ar1-1.0_sp),  a_ar1, lstau, tol, a_ar13, t,x,n)
        if ((abs(a_ar12-a_ar11) > tol2 .and. abs(a_ar12-a_ar1) > tol2) .or. &
            (abs(a_ar13-a_ar11) > tol2 .and. abs(a_ar13-a_ar1) > tol2))     &
           nmu_=1
        dum4=min(dum1,dum2,dum3)
        if (dum4 == dum2) then
           amin=a_ar12
        else if (dum4 == dum3) then
           amin=a_ar13
        else
           amin=a_ar11
        end if
end subroutine minls
!
!===================================================================================================
!
subroutine output
	use own_interfaces
        use result1
        use parameters, only: outputfile
        implicit none
!       Writes output to file.
        character(len=1) :: dummy
        print *
        print '(a,a)','Output written to file ',outputfile
        call info1('w',outputfile)
        call info2('w',outputfile)
        call info3('w',outputfile)
        call info4('w',outputfile)
        print '(a)','Press Enter to exit'
        !read (5,'(a)') dummy
end subroutine output
!
!===================================================================================================
!
subroutine pearsn(x,y,n,r)
        use nrtype
        use nrutil
        implicit none
!       Pearson's correlation coefficient (Numerical Recipes, modified).
        integer, intent(in) :: n
        real(sp), dimension(n), intent(in) :: x,y
        real(sp), intent(out) :: r
        real(sp), parameter :: tiny=1.0e-08_sp
        real(sp), dimension(n) :: xt,yt
        real(sp) :: ax,ay,sxx,sxy,syy
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
!===================================================================================================
!
function phi(x)
        use nrtype
        use nr, only: erfcc
        implicit none
        real(sp) :: phi
        real(sp), intent(in) :: x
!       Cumulative normal density distribution,
!       calculated using NR approximation erfcc.
!       NR: Fractional error [phi] < 1.2e-7.
        phi=1.0_sp-0.5_sp*erfcc(x/sqrt2)
end function phi
!
!===================================================================================================
!
function phi_inv(prob)
        use nrtype
        implicit none
        real(sp) :: phi_inv
!       Inverse cumulative normal density approximation after
!       Odeh and Evans (1974, Appl. Stat. 23, 96-97).
!       Error [phi_inv] < 1.5e-8.
        real(sp), intent(in) :: prob
        real(sp), parameter ::  big=1.0e+38_sp
        real(sp), parameter :: plim=1.0e-20_sp
        real(sp), parameter :: p0=-0.322232431088e0_sp
        real(sp), parameter :: p1=-1.0e0_sp
        real(sp), parameter :: p2=-0.342242088547e0_sp
        real(sp), parameter :: p3=-0.0204231210245e0_sp
        real(sp), parameter :: p4=-0.453642210148e-04_sp
        real(sp), parameter :: q0= 0.0993484626060e0_sp
        real(sp), parameter :: q1= 0.588581570495e0_sp
        real(sp), parameter :: q2= 0.531103462366e0_sp
        real(sp), parameter :: q3= 0.103537752850e0_sp
        real(sp), parameter :: q4= 0.38560700634e-02_sp
        real(sp) :: p=-999.0_sp
        real(sp) :: y=-999.0_sp
        real(sp), parameter :: half=0.5_sp
        real(sp), parameter ::  one=1.0_sp
        real(sp), parameter :: zero=0.0_sp
        phi_inv=-big
        p=prob
!
!       Ranges for p
!       ============
!
        if (p < plim) then
           phi_inv=-big
        else if (p > one-plim) then
           phi_inv=big
        else if (p == half) then
           phi_inv=zero
        else if (p >= plim .and. p < half) then
           y=sqrt(log(1.0_sp/(p*p)))
           phi_inv=y+((((y*p4+p3)*y+p2)*y+p1)*y+p0)/((((y*q4+q3)*y+q2)*y+q1)*y+q0)
           phi_inv=-1.0_sp*phi_inv
        else if (p > half .and. p <= one-plim) then
           p=one-p
           y=sqrt(log(1.0_sp/(p*p)))
           phi_inv=y+((((y*p4+p3)*y+p2)*y+p1)*y+p0)/((((y*q4+q3)*y+q2)*y+q1)*y+q0)
        end if
end function phi_inv
!
!===================================================================================================
!
subroutine plot(type)
        use nrtype
        use data1, only: t1,t2,x1,x2,x3,y1,y2,y3
        use result1, only: r,r_low,r_upp
        use setting, only: n1,n2
        implicit none
        integer, intent(in) :: type             ! Plot type
!
!       Plots data.
!
!       Plot type                   Description
!       =======================================
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
!       ==========
!
        if (type /=1 .and. type /= 2) then
           print *,'Plot-type value out of range - PearsonT terminates.'
           stop
        end if
!
! 2.    Write plot data
!       ===============
!
        open (unit=1, file=file1, status='unknown',                     &
              form='formatted', action='write', iostat=open_error)
        if (open_error /= 0 ) then
           print *,'Error during plot data file opening - PearsonT terminates.'
           stop
        end if
        if (type == 1 .or. type == 2) then
           write (unit=1,fmt='(2x,7(a))')                               	&
                 ' time          ',                                     		&
                 ' x             ',                                     		&
                 ' y             ',                                     		&
                 ' trend(x)      ',                                     		&
                 ' trend(y)      ',                                     		&
                 ' x-detrended   ',                                     		&
                 ' y-detrended   '
           do i=1,n2
              write (unit=1,fmt='(7(1x,es14.6))') t2(i),x2(i),y2(i), &
                                                  x2(i)-x3(i),          		&
                                                  y2(i)-y3(i),          		&
                                                  x3(i),y3(i)
           end do
        end if
        close (unit=1, status='keep')
!
! 3.    Write gnuplot scriptfile
!       ========================
!
!
! 3.1   Key and grid
!       ============
!
        open (unit=2, file=file2, status='unknown',            	&
              form='formatted', action='write', iostat=open_error)
!
!***************************************************************************************************
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
!***************************************************************************************************
!
!       See notes to subroutine plotgn.
!
!***************************************************************************************************
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
!***************************************************************************************************
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
!       ===========
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
!       ==============
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
!       ============
!
        if (type == 1)                                                  		&
           then
           write (unit=2, fmt='(a,a,a,a,a,a,a,a,a)')                            &
                 ' plot ''',file1,''' u 1:2 axes x1y1 title ''x''         w lp 12 ',', ',   &
                 '''1.tmp'' u 1:4 axes x1y1 title ''trend(x) '' w l  11 ',', ',         &
                 '''1.tmp'' u 1:3 axes x1y2 title ''y''         w lp  3 ',', ',           &
                 '''1.tmp'' u 1:5 axes x1y2 title ''trend(y) '' w l   4 '
        else if (type == 2)                                             		 &
           then
           write (unit=2, fmt='(a,a,a)')                                		 &
                 ' plot ''',file1,''' u 6:7 axes x1y1 notitle w p 3'
        end if
!
!       See notes to subroutine plotgn.
!***************************************************************************************************
!       Sun workstation, Unix, X11                    EPCF77 2.7
!       <=> following lines uncommented
 !      write (unit=2, fmt='(1x,a)') 'pause -1'
!***************************************************************************************************
!
!***************************************************************************************************
!       PC, Suse Linux 8.0                            Portland PGF90 3.2-4
!       <=> following lines uncommented
!     write (unit=2, fmt='(1x,a)') 'pause -1'
!***************************************************************************************************
!
!       IBM workstation, Windows XPP 64 bit           Absoft 10.0
!       <=> following lines uncommented
       write (unit=2, fmt='(1x,a)') 'pause mouse key'
       close (unit=2, status ='keep')
!
! 4.    Run gnuplot
!       ===========
!
        call plotgn
end subroutine plot
!
!===================================================================================================
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
!***************************************************************************************************
!       IBM-PC, Windows 9x/NT, DOS window             AbSoft 6.2, AbSoft 10.0
!       <=> following lines uncommented
!         integer*4 plot, system
!         external system
!         plot=system( "gnuplot.exe 2.tmp" )
!***************************************************************************************************
!       IBM-PC, Linux                                 Absoft 6.0
!       <=> following lines uncommented
!       integer*4 plot, system
!       external system
!       plot=system( ' gnuplot 2.tmp ' )
!***************************************************************************************************
!        Sun workstation, Unix, X11                   EPCF77 2.7
!       <=> following lines uncommented
!        call system('/u2/users/mm/gnuplot/gnuplot
!       +             -geometry -5-335 -background white
!       +             /u2/users/mm/xtrend/2.tmp')
! ****************************************************************************
!
!       IBM/Samsung/PC, Windows XP                                    gfortran
!
!       <=> following lines uncommented
!
	!call execute_command_line("wgnuplot.exe 2.tmp")
!***************************************************************************************************
!       PC, Suse Linux 8.0                            Portland PGF90 3.2-4
!       <=> following lines uncommented
!       integer*4 plot, system
!       external system
!       plot=system( ' gnuplot 2.tmp ' )
!***************************************************************************************************
!
end subroutine plotgn
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
        seed=-286761062
        dummy=ran(seed)
end subroutine ranseed
!
!===================================================================================================
!
subroutine read1
        use nrtype
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
!===================================================================================================
!
subroutine r_est
        use nrtype
        use data1, only: t2,x2,y2,x3,y3
        use result1, only: r
        use setting, only: dtrtype,n2
        implicit none
!       Estimates Pearson's correlation coefficient r(x2, y2).
!       x2 and y2 are detrended (linearly or mean), renamed x3, y3, prior to estimation.
        integer :: i
        real(sp) :: a=-999.0_sp
        real(sp) :: b=-999.0_sp
        real(sp) :: mux=-999.0_sp        
        real(sp) :: muy=-999.0_sp 
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
!
end subroutine r_est
!
!===================================================================================================
!
subroutine rhoest(n,x,rho)
        use nrtype
        implicit none
        integer, intent(in) :: n
        real(sp), dimension(n), intent(in) :: x
        real(sp), intent(out) :: rho
!       Estimates autocorrelation coefficient (equidistant data).
        integer i
        real(sp) :: sum1
        real(sp) :: sum2	
        sum1=0.0_sp
        sum2=0.0_sp    
        do i=2,n
           sum1=sum1+x(i)*x(i-1)
           sum2=sum2+x(i)**2.0_sp
        end do
        rho=sum1/sum2
end subroutine rhoest
!
!===================================================================================================
!
subroutine tauest
        use data1, only: t2,x3,y3
        use result1, only: taux3,tauy3,rhox3,rhoy3
        use setting, only: n2
        implicit none
!
        call tauest_x(t2,x3,n2,'x',taux3,rhox3)
        call tauest_x(t2,y3,n2,'y',tauy3,rhoy3)
!
end subroutine tauest
!
!===================================================================================================
!
 subroutine tauest_x(t_in,x_in,n,c,tau,rhoout)
        use nrtype
        use nr, only: avevar
        implicit none
        character (len=1), intent(in) :: c		! x or y
        integer, intent(in) :: n				! number of data points
        real(sp), dimension(n), intent(in) :: t_in	! time
        real(sp), dimension(n), intent(in) :: x_in	! time-series values
        real(sp), intent(out) :: tau			! result: persistence time tau
	real(sp), intent(out) :: rhoout			! result: equivalent autocorrelation coefficient
!
!       Estimates AR(1) persistence time (tau) and
!       equivalent autocorrelation coefficient (rho),using
!       the least-squares algorithm of Mudelsee (2002).
!
!       Notes: (1) assumes t = age (see Point 1)
!
!              (2) automatic bias correction (see Point 5)
!
        real(sp), dimension(n) :: x			! x renamed
        real(sp), dimension(n) :: t			! t renamed
        real(sp) :: avex=-999.0_sp			! average x
        real(sp) :: varx=-999.0_sp			! variance x
        real(sp) :: delta=-999.0_sp			! average spacing
        real(sp) :: scalt=-999.0_sp			!  factor
        real(sp) :: rho=-999.0_sp			! rho x, used in scaling t
        real(sp) :: rho_non=-999.0_sp			! equivalent autocorrelation coefficient
        real(sp) :: amin=-999.0_sp			! output from minls. Estimated value of a=exp(-scalt/tau)
        integer :: mult=-999				! flag (multiple solution) mult=0 if 1 solution,mult=1 if more than one minimum
        integer :: i=-999
	real(sp) :: rho_max = 0.99_sp
	real(sp) :: rho_min = 0.01_sp
!
! 1.    Rename and change time direction
!       ================================
        do i=n,1,-1
           t(i)=-1.0_sp*t_in(n+1-i)
           x(i)=+1.0_sp*x_in(n+1-i)
        end do
!
! 2.    Scaling (x)
!       ===========
        call avevar(x,avex,varx)
        x=x/sqrt(varx)
!
! 3.    Scaling (t)
!       ===========
!       => start value of a = 1/e
!
        delta=abs(t(n)-t(1))/(n-1)		! average sampling		
        call rhoest(n,x,rho)			! rho estimated for x from eq 2.4 in Mudelsee(2010)
        if (rho <= 0.0_sp) then
           rho=0.05_sp
        else if (rho >= 1.0_sp) then		! this rho setting equal to 0.05 or 0.95 is not
           rho=0.95_sp				! problematic since we require rho only for t-scaling
        end if
        scalt=-1.0_sp*log(rho)/delta		! scaling factor
        t=t*scalt				! t-scaled: t*(-log(a)/delta)
!
! 4.    Estimation (tau)
!       ================
	!ls - least square function
	!minls - minimization of least-squares function ls
	!brent - Brent's search - minimum ls value
!	
	call minls(n,t,x,amin,mult)
        print *,'amin=', amin
!
! 5.    Result					! checked if amin (the value for a) is equal or less than zero,
!       ======					! or if it is equal or bigger than 1.0. That can cause problems in the estimation
!
        if (mult == 1 .or. amin <= 0.0_sp .or. amin >= 1.0_sp) then
           if (amin <= 0.0_sp) then
              rho_non=rho_min
	       tau=-delta/log(rho_min)
           else if (amin >= 1.0_sp) then
              rho_non=rho_max
	       tau=-delta/log(rho_max)
           end if
        else
	   tau = -1.0_sp/(scalt*log(amin))		! tau - rescaled, without bias correction tau=-1/log(a)
	   rho_non = exp(-delta/tau)			! equivalent autocorrelation coefficient for tau, used in the bias correction		
	   
	   ! Bias correction(unknown mean) for rho (Kendall, 1954)
	   ! Solved equation (this equation: a_est =a_est_bc+(1.0_sp + 3.0_sp * a_est_bc)/(n-1)
	   ! solved for estimated rho bias corrected = a_est_bc
	   rho_non = (rho_non * (n-1.0_sp) + 1.0_sp) / (n-4.0_sp)
	   
!	   print*,'rho_non',rho_non
!	   print*,'delta', delta
	   
	   if (rho_non >= 1.0_sp)then
	      rho_non = exp(-delta/tau)	! If the bias corrected equivalent autocorrelation coefficient becomes >1 then the bias correction is not performed
	   end if					! it can occur if n is small and the autocorrelation coefficent is large. 
!
!	print*,'bias corrected rho_non',rho_non

	   tau=-delta/log(rho_non)			! tau-calculated from the bias corrected rho 
							! a=exp(delta/tau)  gives tau=-delta/log(rho)
	end if
	
	rhoout = rho_non
!
end subroutine tauest_x
!
!===================================================================================================
!
function t_inv(alpha,dof)
        use nrtype
        use inv_student
        implicit none
        real(sp) :: t_inv
        real(sp) :: phi_inv
!       	Inverse Student's t distribution function, approximation
!       	after Abramowitz and Stegun (1964).
!
!       	Highest inaccuracy occurs when:         o dof is small,
!                                               		    o alpha is high.
!       	For example: dof = 4, alpha = 0.025.
!       	There results: t = 2.7756 (instead of 2.776),
!      	 that is, a negligible error.
!
        real(sp),  intent(in) :: alpha
        integer,  intent(in) :: dof
        real(sp) :: u=-999.0_sp
        integer  :: d=0
        u=phi_inv(alpha)
        d=dof
        t_inv=u                                                            + &
        (1.0_sp/    4)*(                                u**3+    u) / d    + &
        (1.0_sp/   96)*(                    5*u**5+  16*u**3+  3*u) / d**2 + &
        (1.0_sp/  384)*(          3*u**7+  19*u**5+  17*u**3- 15*u) / d**3 + &
        (1.0_sp/92160)*(79*u**9+776*u**7+1482*u**5-1920*u**3-945*u) / d**4
 end function t_inv
!
!===================================================================================================
!
! Used Numerical Recipes files on the computer:
!
include 'avevar.f90'
include 'erfcc.f90'
include 'ran.f90'


