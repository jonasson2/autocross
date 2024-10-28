      LOGICAL FUNCTION LSAME(CA, CB)
      CHARACTER CA, CB

C     LSAME returns .TRUE. if CA is the same letter as CB regardless of case.
C     It returns .FALSE. otherwise. The comparison is case-insensitive.

      LSAME = CA .EQ. CB
      IF (LSAME) RETURN

C     Check if characters are letters and convert to uppercase if they are lowercase.
      IF (CA .GE. 'a' .AND. CA .LE. 'z') CA = CHAR(IACHAR(CA) - 32)
      IF (CB .GE. 'a' .AND. CB .LE. 'z') CB = CHAR(IACHAR(CB) - 32)

C     Perform case-insensitive comparison
      LSAME = CA .EQ. CB

      RETURN
      END
