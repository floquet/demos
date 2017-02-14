      SUBROUTINE ICOPY( N, X, INCX, Y, INCY )
*     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
*     ..
*     .. Array Arguments ..
      INTEGER            X( * ), Y( * )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IX, IY
*     ..
*     .. Executable Statements ..
      IF( N.LE.0 ) THEN
         RETURN
      ENDIF
      IF( ( INCX.EQ.1 ) .AND. ( INCY.EQ.1 ) ) THEN
         DO 10 I = 1, N
            Y( I ) = X( I )
   10    CONTINUE
   20    CONTINUE
      ELSE
         IX = 1
         IY = 1
         DO 30 I = 1, N
            Y( IY ) = X( IX )
            IX = IX + INCX
            IY = IY + INCY
   30    CONTINUE
   40    CONTINUE
      ENDIF
      RETURN
      END
