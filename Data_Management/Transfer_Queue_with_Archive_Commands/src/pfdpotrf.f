      SUBROUTINE PFDPOTRF( UPLO, N, A, IA, JA, DESCA, INFO )
*
*
*  -- ScaLAPACK routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     Oct 10, 1996
*
*
*
*
*  Purpose
*  =======
*
*  PFPOTRF computes the Cholesky factorization of an N-by-N
*  symmetric positive definite distributed matrix sub( A )
*  hermitian positive definite distributed matrix sub( A )
*  denoting A(IA:IA+N-1, JA:JA+N-1).
*
*  The factorization has the form
*
*            sub( A ) = U' * U ,  if UPLO = 'U', or
*
*            sub( A ) = L  * L',  if UPLO = 'L',
*
*  where U is an upper triangular matrix and L is lower triangular.
*
*
*
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DT_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1 )
      INTEGER            CTXT_, M_, N_, MB_, NB_
      PARAMETER          ( CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6 )
      INTEGER            RSRC_, CSRC_
      PARAMETER          ( RSRC_ = 7, CSRC_ = 8 )
      INTEGER            IODEV_, SIZE_, DISK_BLOCK_CYCLIC_2D
      PARAMETER          ( IODEV_ = 10, SIZE_ = 11,
     $                   DISK_BLOCK_CYCLIC_2D = 601 )
      INTEGER            FLEN_
      PARAMETER          ( FLEN_ = 12 )
*     ..
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            IA, INFO, JA, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( FLEN_ )
      DOUBLE PRECISION   A( * )
*     ..
*     .. Local Scalars ..
      LOGICAL            ASIZEQUERY, ISVALID
      INTEGER            ANEED, ANEED0, ANEED1, ASIZE, CSRC, ICONTXT,
     $                   ICTXT, IDUMMY, IODEV, LOCP, LOCQ, M, MB, MM,
     $                   MMB, MYPCOL, MYPROW, NB, NN, NNB, NPCOL, NPROW,
     $                   P0, Q0, RSRC, SAVEDT
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            NUMROC
      EXTERNAL           LSAME, NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, LAIO_INFO, PCHK1MAT,
     $                   PFDLCHFACT, PFDUCHFACT, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
* Check parameters.
*
      INFO = 0
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYPROW,
     $                     MYPCOL )
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 6*100+CTXT_ )
         RETURN
      ENDIF
      ISVALID = ( LSAME( UPLO, 'U' ) .OR. LSAME( UPLO, 'L' ) )
      IF( .NOT.ISVALID ) THEN
         INFO = -1
         CALL PXERBLA( DESCA( CTXT_ ), 'PFDPOTRF', 1 )
         RETURN
      ENDIF
      SAVEDT = DESCA( DT_ )
      DESCA( DT_ ) = BLOCK_CYCLIC_2D
      CALL CHK1MAT( N, 2, N, 2, IA, JA, DESCA, 6, INFO )
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PFDPOTRF', 6 )
         RETURN
      ENDIF
      IDUMMY = 1
      CALL PCHK1MAT( N, 2, N, 2, IA, JA, DESCA, 6, IDUMMY, IDUMMY,
     $               IDUMMY, INFO )
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PFDPOTRF', 6 )
         RETURN
      ENDIF
      DESCA( DT_ ) = SAVEDT
      IF( DESCA( DT_ ).NE.DISK_BLOCK_CYCLIC_2D ) THEN
         INFO = -( 6*100+DT_ )
         CALL PXERBLA( DESCA( CTXT_ ), 'PFDPOTRF', 6 )
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
         CALL PXERBLA( DESCA( CTXT_ ), 'PFDPOTRF', 6 )
         RETURN
      ENDIF
*
*       Need at least 2 panels.
*       Note the use of p0 = myprow, q0 = mypcol
*       to overestimate storage use.
*
      M = N
      P0 = MYPROW
      Q0 = MYPCOL
      IF( LSAME( UPLO, 'L' ) ) THEN
*
*         Lower Triangular part.
*         Panel is m x nnb.
*
         LOCP = NUMROC( M, DESCA( MB_ ), MYPROW, P0, NPROW )
         LOCQ = NUMROC( NNB, DESCA( NB_ ), MYPCOL, Q0, NPCOL )
      ELSE
*
*       Upper triangular part.
*       Panel is mmb x n.
*
         LOCP = NUMROC( MMB, DESCA( MB_ ), MYPROW, P0, NPROW )
         LOCQ = NUMROC( N, DESCA( NB_ ), MYPCOL, Q0, NPCOL )
      ENDIF
      ANEED0 = 2*MAX( 1, LOCP )*MAX( 1, LOCQ )
*
*   Total storage to hold the entire matrix.
*
      P0 = MYPROW
      Q0 = MYPCOL
      LOCP = NUMROC( DESCA( M_ ), DESCA( MB_ ), MYPROW, P0, NPROW )
      LOCQ = NUMROC( DESCA( N_ ), DESCA( NB_ ), MYPCOL, Q0, NPCOL )
      ANEED1 = MAX( LOCP, 1 )*MAX( LOCQ, 1 )
      ANEED = MIN( ANEED0, ANEED1 )
      ASIZE = DESCA( SIZE_ )
      ASIZEQUERY = ( ASIZE.EQ.-1 )
      IF( ASIZE.LT.ANEED ) THEN
         A( 1 ) = ANEED
         INFO = -( 6*100+SIZE_ )
         IF( .NOT.ASIZEQUERY ) THEN
            CALL PXERBLA( DESCA( CTXT_ ), 'PFDPOTRF', 6 )
         ENDIF
         RETURN
*
*       Check storage.
*
      ENDIF
      IODEV = DESCA( IODEV_ )
      ASIZE = DESCA( SIZE_ )
      IF( LSAME( UPLO, 'L' ) ) THEN
         CALL PFDLCHFACT( IODEV, N, IA, JA, A, ASIZE, INFO )
      ELSE
         CALL PFDUCHFACT( IODEV, N, IA, JA, A, ASIZE, INFO )
      ENDIF
      A( 1 ) = ANEED
      RETURN
      END
