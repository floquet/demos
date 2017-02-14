      INTEGER          FUNCTION NUMROCINV( NROW, NB, IPROC, ISRCPROC,
     $                 NPROCS )
*
*
*  -- ScaLAPACK auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     Oct 10, 1996
*
*
* Purpose:
* =======
*
* Inverse operation to NUMROC
*
* Return the largest global (gnrow) such that
* nrow == numroc( gnrow...)
* and numroc( gnrow+1...) > nrow
*
*
*     .. Parameters ..
      INTEGER            IHUGE
      PARAMETER          ( IHUGE = 536870912 )
*     ..
*     .. Scalar Arguments ..
      INTEGER            IPROC, ISRCPROC, NB, NPROCS, NROW
*     ..
*     .. Local Scalars ..
      LOGICAL            ISFOUND, ISOK
      INTEGER            GNROW, GNROWMAX, GNROWMIN, I, NROW1, NROW2
*     ..
*     .. External Functions ..
      INTEGER            INDXL2G, NUMROC
      EXTERNAL           INDXL2G, NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           ASSERT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*
* Use a simple minded algorithm.
* Just try gnrow between gnrowmin (lower bound)
* and gnrowmax (upper bound).
*
*
      IF( NROW.GE.1 ) THEN
         GNROWMIN = INDXL2G( NROW, NB, IPROC, ISRCPROC, NPROCS )
      ELSE
         GNROWMIN = 0
      ENDIF
      GNROWMAX = IHUGE
*
* Use binary search to find gnrow.
*
      ISFOUND = .false.
      DO 10 I = 1, 64
         GNROW = ( GNROWMIN+GNROWMAX ) / 2
         GNROW = MAX( 1, MIN( IHUGE, GNROW ) )
         NROW1 = NUMROC( GNROW, NB, IPROC, ISRCPROC, NPROCS )
         NROW2 = NUMROC( GNROW+1, NB, IPROC, ISRCPROC, NPROCS )
         ISFOUND = ( NROW1.EQ.NROW ) .AND. ( NROW2.GT.NROW )
         IF( ISFOUND ) THEN
            GOTO 20
         ENDIF
         IF( NROW1.GT.NROW ) THEN
            GNROWMAX = GNROW
         ELSE
            GNROWMIN = GNROW
         ENDIF
   10 CONTINUE
   20 CONTINUE
      IF( .NOT.ISFOUND ) THEN
         DO 30 GNROW = GNROWMIN, GNROWMAX
            NROW1 = NUMROC( GNROW, NB, IPROC, ISRCPROC, NPROCS )
            NROW2 = NUMROC( GNROW+1, NB, IPROC, ISRCPROC, NPROCS )
            ISFOUND = ( NROW1.EQ.NROW ) .AND. ( NROW2.GT.NROW )
            IF( ISFOUND ) THEN
               GOTO 40
            ENDIF
   30    CONTINUE
   40    CONTINUE
      ENDIF
      CALL ASSERT( ISFOUND, ' ** numrocinv, gnrow ', GNROW )
*
* Double check again.
*
      ISOK = ( NUMROC( GNROW, NB, IPROC, ISRCPROC, NPROCS ).EQ.NROW )
     $        .AND. ( NUMROC( GNROW+1, NB, IPROC, ISRCPROC, NPROCS ).GT.
     $       NROW )
      CALL ASSERT( ISOK, '** numrocinv ** ', GNROW )
      NUMROCINV = ( GNROW )
      RETURN
      END
