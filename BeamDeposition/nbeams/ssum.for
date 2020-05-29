C  Purpose:    Sum the values of a single precision vector.
C
C  Usage:      SSUM(N, SX, INCX)
C
C  Arguments:
C     N      - Length of vectors X.  (Input)
C     SX     - Real vector of length N*INCX.  (Input)
C     INCX   - Displacement between elements of SX.  (Input)
C              X(I) is defined to be SX(1+(I-1)*INCX). INCX must be
C              greater than 0.
C     SSUM   - Single precision sum from I=1 to N of X(I).  (Output)
C              X(I) refers to a specific element of SX.
C
C-----------------------------------------------------------------------
C
      REAL FUNCTION SSUM (N, SX, INCX)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    N, INCX
      REAL       SX(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, M, MP1, NINCX
C                                  SPECIFICATIONS FOR SPECIAL CASES
C     INTRINSIC  MOD
      INTRINSIC  MOD
      INTEGER    MOD
C
      SSUM = 0.0E0
      IF (N .GT. 0) THEN
         IF (INCX .NE. 1) THEN
C                                  CODE FOR INCREMENT NOT EQUAL TO 1
            NINCX = N*INCX
            DO 10  I=1, NINCX, INCX
               SSUM = SSUM + SX(I)
   10       CONTINUE
         ELSE
C                                  CODE FOR INCREMENT EQUAL TO 1
            M = MOD(N,6)
C                                  CLEAN-UP LOOP
            DO 30  I=1, M
               SSUM = SSUM + SX(I)
   30       CONTINUE
            MP1 = M + 1
            DO 40  I=MP1, N, 6
               SSUM = SSUM + SX(I) + SX(I+1) + SX(I+2) + SX(I+3) +
     &                SX(I+4) + SX(I+5)
   40       CONTINUE
         END IF
      END IF
      RETURN
      END
