      SUBROUTINE PFMAX2SIZE( ASIZE, M, N, MM, NN, IODEV )
*
*
*  -- ScaLAPACK auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     Oct 10, 1996
*
*
*
* Purpose:
* ========
* Compute the appropriate size of matrix patch to work with
* out of core algorithms given the size of problem (m by n)
* and amount of available storage (Asize)
*
* intent(out) :: mm,nn
*
*
*     .. Scalar Arguments ..
      INTEGER            ASIZE, IODEV, M, MM, N, NN
*     ..
*     .. Local Scalars ..
      LOGICAL            ASIZEQUERY, ISVALID
      CHARACTER          ROWORCOL
      INTEGER            CSRC, ICONTXT, INFO, J, MB, MMB, MYPCOL,
     $                   MYPROW, NB, NNB, NPCOL, NPROW, NTRY, RSRC
*     ..
*     .. External Functions ..
      INTEGER            ICEIL
      EXTERNAL           ICEIL
*     ..
*     .. External Subroutines ..
      EXTERNAL           ASSERT, BLACS_GRIDINFO, LAIO_INFO, PFMAXSIZE
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
      CALL LAIO_INFO( IODEV, MM, NN, MMB, NNB, MB, NB, CSRC, RSRC,
     $                ICONTXT )
      CALL BLACS_GRIDINFO( ICONTXT, NPROW, NPCOL, MYPROW, MYPCOL )
      NTRY = MAX( 1, ICEIL( M, MMB ) )
      ASIZEQUERY = ( ASIZE.EQ.-1 )
      ROWORCOL = 'C'
      DO 10 J = NTRY, 1, -1
         MM = MIN( M, J*MMB )
         CALL PFMAXSIZE( ROWORCOL, ASIZE, MM, NN, MB, NB, RSRC, CSRC,
     $                   ICONTXT, INFO )
         ISVALID = ( INFO.EQ.0 ) .OR. ( ( INFO.EQ.-2 ) .AND.
     $             ( ASIZEQUERY ) )
         CALL ASSERT( ISVALID, '** pfmax2size: pfmaxsize returns info ',
     $                INFO )
         IF( NN.GE.NNB ) THEN
            NN = NN - MOD( NN, NNB )
            RETURN
         ENDIF
   10 CONTINUE
   20 CONTINUE
* try again
      MM = MIN( M, MB*NPROW )
      ROWORCOL = 'C'
      CALL PFMAXSIZE( ROWORCOL, ASIZE, MM, NN, MB, NB, RSRC, CSRC,
     $                ICONTXT, INFO )
      ISVALID = ( INFO.EQ.0 ) .OR. ( ( INFO.EQ.-2 ) .AND.
     $          ( ASIZEQUERY ) )
      CALL ASSERT( ISVALID, '** pfmax2size: pfmaxsize returns info ',
     $             INFO )
      IF( NN.GE.1 ) THEN
         RETURN
      ENDIF
* last try
      NN = MIN( N, NB*NPCOL )
      ROWORCOL = 'R'
      CALL PFMAXSIZE( ROWORCOL, ASIZE, MM, NN, MB, NB, RSRC, CSRC,
     $                ICONTXT, INFO )
      ISVALID = ( INFO.EQ.0 ) .OR. ( ( INFO.EQ.-2 ) .AND.
     $          ( ASIZEQUERY ) )
      CALL ASSERT( INFO.EQ.0, '** pfmax2size: pfmaxsize returns info ',
     $             INFO )
      IF( MM.GE.1 ) THEN
         RETURN
      ENDIF
* Treat as error condition
      MM = -1
      NN = -1
      RETURN
      END
