      SUBROUTINE PFDTRSM( SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A, IA,
     $                    JA, DESCA, B, IB, JB, DESCB, WORK, LWORK,
     $                    INFO )
*
*
*  -- ScaLAPACK auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     Oct 10, 1996
*
*
*   Purpose:
*   ========
*
*   PFxTRSM performs a triangular solve, similar to
*   PBLAS PxTRSM  routine.
*
*   A is assumed to be an out-of-core matrix.
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
      CHARACTER          DIAG, SIDE, TRANS, UPLO
      INTEGER            IA, IB, INFO, JA, JB, LWORK, M, N
      DOUBLE PRECISION   ALPHA
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( FLEN_ ), DESCB( DLEN_ )
      DOUBLE PRECISION   A( * ), B( * ), WORK( * )
*     ..
*     .. Local Scalars ..
      LOGICAL            ASIZEQUERY, ISLEFT, ISVALID
      INTEGER            ANEED, ASIZE, CSRC, ICONTXT, ICTXT, IDUMMY,
     $                   LOCP, LOCQ, MB, MM, MMB, MYPCOL, MYPROW, NB,
     $                   NN, NNB, NPCOL, NPROW, P0, Q0, RSRC, SAVEDT
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            NUMROC
      EXTERNAL           LSAME, NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, LAIO_INFO, PCHK1MAT,
     $                   PFDLATRSM, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
* ===========================================================
      INFO = 0
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYPROW,
     $                     MYPCOL )
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 11*100+CTXT_ )
         RETURN
      ENDIF
      ICTXT = DESCB( CTXT_ )
      CALL BLACS_GRIDINFO( DESCB( CTXT_ ), NPROW, NPCOL, MYPROW,
     $                     MYPCOL )
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 15*100+CTXT_ )
         RETURN
      ENDIF
      ICONTXT = DESCA( CTXT_ )
      ISVALID = ( LSAME( SIDE, 'R' ) .OR. LSAME( SIDE, 'L' ) )
      IF( .NOT.ISVALID ) THEN
         INFO = -1
         CALL PXERBLA( ICONTXT, 'PFDTRSM', 1 )
         RETURN
      ENDIF
      ISVALID = ( LSAME( TRANS, 'T' ) .OR. LSAME( TRANS, 'C' ) .OR.
     $          LSAME( TRANS, 'N' ) )
      IF( .NOT.ISVALID ) THEN
         INFO = -3
         CALL PXERBLA( ICONTXT, 'PFDTRSM', 3 )
         RETURN
      ENDIF
      ISVALID = ( LSAME( UPLO, 'U' ) .OR. LSAME( UPLO, 'L' ) )
      IF( .NOT.ISVALID ) THEN
         INFO = -2
         CALL PXERBLA( ICONTXT, 'PFDTRSM', 2 )
         RETURN
      ENDIF
      ISVALID = ( LSAME( DIAG, 'U' ) .OR. LSAME( DIAG, 'N' ) )
      IF( .NOT.ISVALID ) THEN
         INFO = -4
         CALL PXERBLA( ICONTXT, 'PFDTRSM', 4 )
         RETURN
      ENDIF
      WORK( 1 ) = 0
*
*
*   sub(A)  is a square m x m matrix if side = 'L'
*   sub(A)  is a square n x n matrix if side = 'R'
*
*   sub(B)  is a m x n matrix
*
*
      ISLEFT = LSAME( SIDE, 'L' )
      IF( ISLEFT ) THEN
         SAVEDT = DESCA( DT_ )
         DESCA( DT_ ) = BLOCK_CYCLIC_2D
         CALL CHK1MAT( M, 5, M, 5, IA, JA, DESCA, 11, INFO )
         IF( INFO.NE.0 ) THEN
            CALL PXERBLA( DESCA( CTXT_ ), 'PFDTRSM', 11 )
            RETURN
         ENDIF
         IDUMMY = 1
         CALL PCHK1MAT( M, 5, M, 5, IA, JA, DESCA, 11, IDUMMY, IDUMMY,
     $                  IDUMMY, INFO )
         IF( INFO.NE.0 ) THEN
            CALL PXERBLA( DESCA( CTXT_ ), 'PFDTRSM', 11 )
            RETURN
         ENDIF
         DESCA( DT_ ) = SAVEDT
      ELSE
         SAVEDT = DESCA( DT_ )
         DESCA( DT_ ) = BLOCK_CYCLIC_2D
         CALL CHK1MAT( N, 6, N, 6, IA, JA, DESCA, 11, INFO )
         IF( INFO.NE.0 ) THEN
            CALL PXERBLA( DESCA( CTXT_ ), 'PFDTRSM', 11 )
            RETURN
         ENDIF
         IDUMMY = 1
         CALL PCHK1MAT( N, 6, N, 6, IA, JA, DESCA, 11, IDUMMY, IDUMMY,
     $                  IDUMMY, INFO )
         IF( INFO.NE.0 ) THEN
            CALL PXERBLA( DESCA( CTXT_ ), 'PFDTRSM', 11 )
            RETURN
         ENDIF
         DESCA( DT_ ) = SAVEDT
      ENDIF
      SAVEDT = DESCB( DT_ )
      DESCB( DT_ ) = BLOCK_CYCLIC_2D
      CALL CHK1MAT( M, 5, N, 6, IB, JB, DESCB, 15, INFO )
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCB( CTXT_ ), 'PFDTRSM', 15 )
         RETURN
      ENDIF
      IDUMMY = 1
      CALL PCHK1MAT( M, 5, N, 6, IB, JB, DESCB, 15, IDUMMY, IDUMMY,
     $               IDUMMY, INFO )
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCB( CTXT_ ), 'PFDTRSM', 15 )
         RETURN
      ENDIF
      DESCB( DT_ ) = SAVEDT
      IF( DESCA( DT_ ).NE.DISK_BLOCK_CYCLIC_2D ) THEN
         INFO = -( 11*100+DT_ )
         CALL PXERBLA( DESCA( CTXT_ ), 'PFDTRSM', 11 )
         RETURN
      ENDIF
      INFO = 0
      CALL LAIO_INFO( DESCA( IODEV_ ), MM, NN, MMB, NNB, MB, NB, CSRC,
     $                RSRC, ICONTXT )
      IF( MB.NE.DESCA( MB_ ) ) THEN
         INFO = -( 11*100+MB_ )
      ENDIF
      IF( NB.NE.DESCA( NB_ ) ) THEN
         INFO = -( 11*100+NB_ )
      ENDIF
      IF( ICONTXT.NE.DESCA( CTXT_ ) ) THEN
         INFO = -( 11*100+CTXT_ )
      ENDIF
      IF( DESCA( RSRC_ ).NE.RSRC ) THEN
         INFO = -( 11*100+RSRC_ )
      ENDIF
      IF( DESCA( CSRC_ ).NE.CSRC ) THEN
         INFO = -( 11*100+CSRC_ )
      ENDIF
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PFDTRSM', 11 )
         RETURN
      ENDIF
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYPROW,
     $                     MYPCOL )
*
*       Need at least one panel.
*
      IF( LSAME( SIDE, 'L' ) ) THEN
         MM = M
         NN = MIN( NNB, M )
      ELSE
         MM = N
         NN = MIN( NNB, N )
      ENDIF
*
*       Note the use of p0,q0 to overestimate storage.
*
      P0 = MYPROW
      Q0 = MYPCOL
      LOCP = NUMROC( MM, MB, MYPROW, P0, NPROW )
      LOCQ = NUMROC( NN, NB, MYPCOL, Q0, NPCOL )
      ANEED = LOCP*LOCQ
      ASIZE = DESCA( SIZE_ )
      ASIZEQUERY = ( ASIZE.EQ.-1 )
      IF( ANEED.GT.ASIZE ) THEN
         INFO = -( 11*100+SIZE_ )
         A( 1 ) = ANEED
         IF( .NOT.ASIZEQUERY ) THEN
            CALL PXERBLA( DESCA( CTXT_ ), 'PFDTRSM', 11 )
         ENDIF
         RETURN
      ENDIF
      CALL PFDLATRSM( SIDE, UPLO, TRANS, DIAG, DESCA( IODEV_ ), M, N,
     $                IA, JA, ALPHA, B, IB, JB, DESCB, A, ASIZE, INFO )
      IF( INFO.EQ.-16 ) THEN
         INFO = -( 11*100+SIZE_ )
         IF( .NOT.ASIZEQUERY ) THEN
            CALL PXERBLA( DESCA( CTXT_ ), 'PFDTRSM', 11 )
         ENDIF
         RETURN
      ENDIF
      END
