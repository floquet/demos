      SUBROUTINE PFDGETRF( M, N, A, IA, JA, DESCA, IPIV, INFO )
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
*   Purpose
*   =======
*
*   PZGETRF computes an LU factorization of a general M-by-N distributed
*   matrix sub( A ) = A(IA:IA+M-1,JA:JA+N-1) using partial pivoting with
*   row interchanges.
*
*   The factorization has the form sub( A ) = P * L * U, where P is a
*   permutation matrix, L is lower triangular with unit diagonal ele-
*   ments (lower trapezoidal if m > n), and U is upper triangular
*   (upper trapezoidal if m < n). L and U are stored in sub( A ).
*
*   This is the left-looking Parallel out-of-core version of the
*   algorithm.
*
*   For details, see routine PxGETRF.
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
      INTEGER            LIWORK
      PARAMETER          ( LIWORK = 100*1000 )
*     ..
*     .. Scalar Arguments ..
      INTEGER            IA, INFO, JA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( FLEN_ ), IPIV( * )
      DOUBLE PRECISION   A( * )
*     ..
*     .. Local Scalars ..
      INTEGER            ANEED, CSRC, ICONTXT, ICTXT, IDUMMY, LOCP,
     $                   LOCQ, MB, MM, MMB, MYPCOL, MYPROW, NB, NCOPY,
     $                   NN, NNB, NPCOL, NPROW, RSRC, SAVEDT
*     ..
*     .. External Functions ..
      INTEGER            NUMROC
      EXTERNAL           NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, LAIO_INFO, PCHK1MAT,
     $                   PFDGETF2, PXERBLA
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
      SAVEDT = DESCA( DT_ )
      DESCA( DT_ ) = BLOCK_CYCLIC_2D
      CALL CHK1MAT( M, 1, N, 2, IA, JA, DESCA, 6, INFO )
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PFDGETRF', 6 )
         RETURN
      ENDIF
      IDUMMY = 1
      CALL PCHK1MAT( M, 1, N, 2, IA, JA, DESCA, 6, IDUMMY, IDUMMY,
     $               IDUMMY, INFO )
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PFDGETRF', 6 )
         RETURN
      ENDIF
      DESCA( DT_ ) = SAVEDT
      IF( DESCA( DT_ ).NE.DISK_BLOCK_CYCLIC_2D ) THEN
         INFO = -( 6*100+DT_ )
         CALL PXERBLA( DESCA( CTXT_ ), 'PFDGETRF', 6 )
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
         CALL PXERBLA( DESCA( CTXT_ ), 'PFDGETRF', 6 )
         RETURN
      ENDIF
      NCOPY = 2
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYPROW,
     $                     MYPCOL )
      LOCP = NUMROC( M, DESCA( MB_ ), MYPROW, MYPROW, NPROW )
      LOCQ = NUMROC( NNB, DESCA( NB_ ), MYPCOL, MYPCOL, NPCOL )
      ANEED = NCOPY*( LOCP*LOCQ )
      IF( DESCA( SIZE_ ).LT.ANEED ) THEN
         A( 1 ) = ANEED
         INFO = -( 6*100+SIZE_ )
         IF( DESCA( SIZE_ ).NE.-1 ) THEN
            CALL PXERBLA( DESCA( CTXT_ ), 'PFDGETRF', 6 )
         ENDIF
         RETURN
      ENDIF
      CALL PFDGETF2( M, N, A, IA, JA, DESCA, IPIV, INFO )
      A( 1 ) = ANEED
      RETURN
      END
