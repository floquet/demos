      SUBROUTINE PFCLUSOLVE( IODEV, M, NRHS, IA, JA, DESCA, IPIV, B, IB,
     $                       JB, DESCB, WORK, LWORK )
*
*
*  -- ScaLAPACK auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     Oct 10, 1996
*
*
*
*   Purpose
*   =======
*
*   PFLUSOLVE is an auxiliary routine called by PFGETRF
*   that computes an LU factorization of a general M-by-N distributed
*   matrix sub( A ) = (IA:IA+M-1,JA:JA+N-1) using partial pivoting with
*   row interchanges.
*
*   The factorization has the form sub( A ) = P * L * U, where P is a
*   permutation matrix, L is lower triangular with unit diagonal ele-
*   ments (lower trapezoidal if m > n), and U is upper triangular
*   (upper trapezoidal if m < n). L and U are stored in sub( A ).
*
*   This is the left-looking out-of-core version of the
*   algorithm.
*
*
*     .. Parameters ..
      INTEGER            DLEN_
      PARAMETER          ( DLEN_ = 9 )
      INTEGER            CTXT_, M_, N_, MB_
      PARAMETER          ( CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5 )
      INTEGER            RSRC_
      PARAMETER          ( RSRC_ = 7 )
*     ..
*     .. Scalar Arguments ..
      INTEGER            IA, IB, IODEV, JA, JB, LWORK, M, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( DLEN_ ), DESCB( DLEN_ ), IPIV( * )
      COMPLEX            B( * ), WORK( * )
*     ..
*     .. Local Scalars ..
      LOGICAL            ISVALID
      CHARACTER          DIAG, SIDE, TRANS, UPLO
      INTEGER            IAONTEXT, IDUM1, IFREE, INFO, LDIPIV, MM, MYID,
     $                   MYPCOL, MYPROW, NFREE, NN, NPCOL, NPROC, NPROW,
     $                   NROW, P0, Q0
      COMPLEX            DALPHA, ONE, ZERO
*     ..
*     .. Local Arrays ..
      INTEGER            DESCIPIV( DLEN_ )
*     ..
*     .. External Functions ..
      INTEGER            NUMROC
      EXTERNAL           NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           ASSERT, BLACS_GRIDINFO, BLACS_PINFO, DESCINIT,
     $                   PCLAPIV, PFCLATRSM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, REAL
*     ..
*     .. Executable Statements ..
*
*       Solve sub( C ) *X = sub( B ) where sub( C ) = L*U.
*       L is unit diagonal.
*
      ONE = CMPLX( REAL( 1 ) )
      ZERO = CMPLX( REAL( 0 ) )
      CALL BLACS_PINFO( MYID, NPROC )
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYPROW,
     $                     MYPCOL )
      IFREE = 1
      NFREE = LWORK
      ISVALID = ( 1.LE.IA ) .AND. ( 1.LE.JA ) .AND. ( 1.LE.M ) .AND.
     $          ( 1.LE.NRHS ) .AND. ( M.LE.DESCB( M_ ) ) .AND.
     $          ( NRHS.LE.DESCB( N_ ) )
      CALL ASSERT( ISVALID, ' ** PFLUSOLVE: invalid parameters ', MYID )
*
*       Solve L*X = sub( B ), overwriting sub( B ) with X.
*
      NROW = NUMROC( DESCA( M_ ), DESCA( MB_ ), MYPROW, DESCA( RSRC_ ),
     $       NPROW )
      MM = DESCA( M_ ) + DESCA( MB_ )*NPROW
      NN = 1
      LDIPIV = NROW + DESCA( MB_ )
      P0 = DESCA( RSRC_ )
      Q0 = MYPCOL
      IAONTEXT = DESCA( CTXT_ )
      INFO = 0
      CALL DESCINIT( DESCIPIV, MM, NN, DESCA( MB_ ), 1, P0, Q0,
     $               IAONTEXT, LDIPIV, INFO )
      CALL ASSERT( INFO.EQ.0, '** PFLUSOLVE: descinit returns info ',
     $             INFO )
      CALL PCLAPIV( 'Forward', 'Row', 'Col', M, NRHS, B, IB, JB, DESCB,
     $              IPIV, IA, 1, DESCIPIV, IDUM1 )
      DALPHA = ONE
      UPLO = 'L'
      TRANS = 'N'
      DIAG = 'U'
      SIDE = 'L'
      CALL PFCLATRSM( SIDE, UPLO, TRANS, DIAG, IODEV, M, NRHS, IA, JA,
     $                DALPHA, B, IB, JB, DESCB, WORK( IFREE ), NFREE,
     $                INFO )
*
*       Solve U*X = sub( B ), overwriting sub( B ) with X.
*
      DALPHA = ONE
      UPLO = 'U'
      TRANS = 'N'
      DIAG = 'N'
      SIDE = 'L'
      CALL PFCLATRSM( SIDE, UPLO, TRANS, DIAG, IODEV, M, NRHS, IA, JA,
     $                DALPHA, B, IB, JB, DESCB, WORK( IFREE ), NFREE,
     $                INFO )
      RETURN
      END
