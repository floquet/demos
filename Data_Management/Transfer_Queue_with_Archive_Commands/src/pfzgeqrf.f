      SUBROUTINE PFZGEQRF( M, N, A, IA, JA, DESCA, TAU, WORK, LWORK,
     $                     INFO )
*
*
*  -- ScaLAPACK routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     Oct 10, 1996
*
*
*
* Purpose:
* ========
*
* Perform out-of-core QR factorization.
* Argument list is similar to ScaLAPACK PxGEQRF routine.
*
* Note descriptor for matrix A has extra fields to identify it
* as an out-of-core matrix.
*
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DT_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1 )
      INTEGER            CTXT_, MB_, NB_
      PARAMETER          ( CTXT_ = 2, MB_ = 5, NB_ = 6 )
      INTEGER            RSRC_, CSRC_
      PARAMETER          ( RSRC_ = 7, CSRC_ = 8 )
      INTEGER            IODEV_, SIZE_, DISK_BLOCK_CYCLIC_2D
      PARAMETER          ( IODEV_ = 10, SIZE_ = 11,
     $                   DISK_BLOCK_CYCLIC_2D = 601 )
      INTEGER            FLEN_
      PARAMETER          ( FLEN_ = 12 )
*     ..
*     .. Scalar Arguments ..
      INTEGER            IA, INFO, JA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( FLEN_ )
      COMPLEX*16         A( * ), TAU( * ), WORK( * )
*     ..
*     .. Local Scalars ..
      LOGICAL            ASIZEQUERY, ISOK, LWORKQUERY
      INTEGER            ANEED, ASIZE, CSRC, ICONTXT, ICTXT, IDUMMY,
     $                   INEED, IODEV, LOCP, LOCQ, LOCQ1, LOCQ2, LWORK0,
     $                   LWORK1, MB, MM, MMB, MP0, MYPCOL, MYPROW, NB,
     $                   NCOL, NN, NNB, NPCOL, NPROW, NQ0, P0, Q0, RSRC,
     $                   SAVEDT, SIZE1, SIZE2
*     ..
*     .. Local Arrays ..
      INTEGER            DESC1( DLEN_ ), DESC2( DLEN_ )
*     ..
*     .. External Functions ..
      INTEGER            NUMROC
      EXTERNAL           NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           ASSERT, BLACS_GRIDINFO, CHK1MAT, DESCINIT,
     $                   LAIO_INFO, PCHK1MAT, PFMAXSIZE, PFZQRFACT2,
     $                   PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
      INFO = 0
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYPROW,
     $                     MYPCOL )
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 6*100+CTXT_ )
         RETURN
      ENDIF
      SAVEDT = DESCA( DT_ )
      DESCA( DT_ ) = BLOCK_CYCLIC_2D
      CALL CHK1MAT( M, 1, N, 2, IA, JA, DESCA, 6, INFO )
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PFZGEQRF', 6 )
         RETURN
      ENDIF
      IDUMMY = 1
      CALL PCHK1MAT( M, 1, N, 2, IA, JA, DESCA, 6, IDUMMY, IDUMMY,
     $               IDUMMY, INFO )
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PFZGEQRF', 6 )
         RETURN
      ENDIF
      DESCA( DT_ ) = SAVEDT
      IF( DESCA( DT_ ).NE.DISK_BLOCK_CYCLIC_2D ) THEN
         INFO = -( 6*100+DT_ )
         CALL PXERBLA( DESCA( CTXT_ ), 'PFZGEQRF', 6 )
         RETURN
      ENDIF
      INFO = 0
      CALL LAIO_INFO( DESCA( IODEV_ ), MM, NN, MMB, NNB, MB, NB, CSRC,
     $                RSRC, ICONTXT )
      IF( MB.NE.DESCA( MB_ ) ) THEN
         INFO = -( 6*100+MB_ )
      ENDIF
      IF( NB.NE.DESCA( NB_ ) ) THEN
         INFO = -( 6*100+NB_ )
      ENDIF
      IF( ICONTXT.NE.DESCA( CTXT_ ) ) THEN
         INFO = -( 6*100+CTXT_ )
      ENDIF
      IF( DESCA( RSRC_ ).NE.RSRC ) THEN
         INFO = -( 6*100+RSRC_ )
      ENDIF
      IF( DESCA( CSRC_ ).NE.CSRC ) THEN
         INFO = -( 6*100+CSRC_ )
      ENDIF
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PFZGEQRF', 6 )
         RETURN
      ENDIF
      IODEV = DESCA( IODEV_ )
      ASIZE = DESCA( SIZE_ )
      CALL LAIO_INFO( IODEV, MM, NN, MMB, NNB, MB, NB, CSRC, RSRC,
     $                ICTXT )
      LWORKQUERY = ( LWORK.EQ.-1 )
      ASIZEQUERY = ( ASIZE.EQ.-1 )
*
*  Calculate storage for buffer.
*  At least same amount of in core panel.
*
      P0 = MYPROW
      Q0 = MYPCOL
      LOCP = NUMROC( M, MB, MYPROW, P0, NPROW )
      LOCQ1 = NUMROC( NNB, NB, MYPCOL, Q0, NPCOL )
      SIZE1 = LOCP*LOCQ1
      ANEED = 2*SIZE1
*
*  Sufficient storage to accomodate two panels.
*  Try to do better by letting Panel 2 takes the rest of storage.
*
      NCOL = NNB
      IF( ANEED.LE.ASIZE ) THEN
         NCOL = -1
         CALL PFMAXSIZE( 'C', ASIZE-SIZE1, M, NCOL, MB, NB, RSRC, CSRC,
     $                   ICTXT, INFO )
         CALL ASSERT( INFO.EQ.0, '** pfmaxsize returns info ', INFO )
         ISOK = ( NCOL.GE.NNB )
         IF( .NOT.ISOK ) THEN
            WRITE( *, FMT = * )'** PFGEQRF: Asize,size1, ncol,nnb ',
     $         ASIZE, SIZE1, NCOL, NNB
         ENDIF
         CALL ASSERT( NCOL.GE.NNB,
     $                '**  PFGEQRF: internal error ncol < nnb, ncol ',
     $                NCOL )
         NCOL = MAX( 1, MIN( NCOL, N ) )
         LOCQ2 = NUMROC( NCOL, NB, MYPCOL, CSRC, NPCOL )
         SIZE2 = LOCP*LOCQ2
      ENDIF
*
*  Check amount of temp storage required.
*  Overestimate storage slightly for safety margin.
*
*
*    lwork0 storage for PxGEQRF
*
      MP0 = NUMROC( M+MB-1, MB, MYPROW, P0, NPROW ) + MB
      NQ0 = NUMROC( NCOL+NB-1, NB, MYPCOL, Q0, NPCOL ) + NB
      LWORK0 = NB*( MP0+NQ0+NB )
*
*     extra storage for intermediate "tau" vectors.
*
      LOCQ = NUMROC( N, NB, MYPCOL, CSRC, NPCOL )
      LWORK1 = 2*LOCQ + MAX( NNB, NCOL )
      INEED = LWORK0 + LWORK1
      IF( ASIZE.LT.ANEED ) THEN
* need more storage to proceed effectively
         A( 1 ) = ANEED
         WORK( 1 ) = INEED
         INFO = -( 6*100+SIZE_ )
         IF( .NOT.ASIZEQUERY ) THEN
            CALL PXERBLA( DESCA( CTXT_ ), 'PFZGEQRF', 6 )
         ENDIF
         RETURN
      ENDIF
      IF( LWORK.LT.INEED ) THEN
         A( 1 ) = ANEED
         WORK( 1 ) = INEED
         INFO = -9
         IF( .NOT.LWORKQUERY ) THEN
            CALL PXERBLA( DESCA( CTXT_ ), 'PFZGEQRF', 9 )
         ENDIF
         RETURN
      ENDIF
      CALL DESCINIT( DESC1, M, NNB, MB, NB, RSRC, CSRC, ICTXT,
     $               MAX( 1, LOCP ), INFO )
      CALL ASSERT( INFO.EQ.0, '** PFGEQRF: descinit returns info ',
     $             INFO )
      CALL DESCINIT( DESC2, M, NCOL, MB, NB, RSRC, CSRC, ICTXT,
     $               MAX( 1, LOCP ), INFO )
      CALL ASSERT( INFO.EQ.0, '** PFGEQRF: descinit returns info ',
     $             INFO )
      CALL PFZQRFACT2( IODEV, M, N, A( 1 ), IA, JA, DESCA, TAU, WORK,
     $                 LWORK )
      WORK( 1 ) = INEED
      A( 1 ) = ANEED
      INFO = 0
      RETURN
      END
