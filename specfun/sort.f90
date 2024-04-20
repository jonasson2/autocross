module sort_module
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
  ! The only other change is commenting out of lsame declaration
  ! and xerbla
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
    DOUBLEPRECISION D ( * ) 
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
    DOUBLEPRECISION D1, D2, D3, DMNMX, TMP 
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
