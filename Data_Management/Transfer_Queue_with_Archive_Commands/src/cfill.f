      SUBROUTINE CFILL( N, DVALUE, A, INCX )
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      COMPLEX            DVALUE
*     ..
*     .. Array Arguments ..
      COMPLEX            A( * )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IX
*     ..
*     .. Executable Statements ..
      IF( N.LE.0 ) THEN
         RETURN
      ENDIF
      IF( INCX.EQ.1 ) THEN
* special case
         DO 10 I = 1, N
            A( I ) = DVALUE
   10    CONTINUE
   20    CONTINUE
      ELSE
* more general case
         IX = 1
         DO 30 I = 1, N
            A( IX ) = DVALUE
            IX = IX + INCX
   30    CONTINUE
   40    CONTINUE
      ENDIF
      RETURN
      END
