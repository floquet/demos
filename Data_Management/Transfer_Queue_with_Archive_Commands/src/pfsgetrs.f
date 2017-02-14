      SUBROUTINE PFSGETRS( TRANS, N, NRHS, A, IA, JA, DESCA, IPIV, B,
     $                     IB, JB, DESCB, INFO )
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
*   PFGETRS solves a system of distributed linear equations
*
*                    op( sub( A ) ) * X = sub( B )
*
*   with a general N-by-N distributed matrix sub( A ) using the LU
*   factorization computed by PFGETRF.
*
*   sub( A ) denotes A(IA:IA+N-1,JA:JA+N-1), an out-of-core matrix
*   op( A ) = A, A**T or A**H
*   and sub( B ) denotes B(IB:IB+N-1,JB:JB+NRHS-1), an in-core matrix.
*
*
*
*
*  See ScaLAPACK PGETRS for more details.
*
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DT_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1 )
      INTEGER            CTXT_, M_, MB_, NB_
      PARAMETER          ( CTXT_ = 2, M_ = 3, MB_ = 5, NB_ = 6 )
      INTEGER            RSRC_, CSRC_
      PARAMETER          ( RSRC_ = 7, CSRC_ = 8 )
      INTEGER            IODEV_, SIZE_, DISK_BLOCK_CYCLIC_2D
      PARAMETER          ( IODEV_ = 10, SIZE_ = 11,
     $                   DISK_BLOCK_CYCLIC_2D = 601 )
      INTEGER            FLEN_
      PARAMETER          ( FLEN_ = 12 )
      INTEGER            LWORK
      PARAMETER          ( LWORK = 100*100 )
*     ..
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( FLEN_ ), DESCB( DLEN_ ), IPIV( * )
      REAL               A( * ), B( * )
*     ..
*     .. Local Scalars ..
      LOGICAL            ASIZEQUERY, HASWORK, ISVALID, NOTRAN
      INTEGER            ANEED, ANEED0, ANEED1, ASIZE, CSRC, ICONTXT,
     $                   ICTXT, IDUMMY, INEED, LDIPIV, LOCP, LOCQ, M,
     $                   MB, MIPIV, MM, MMB, MYPCOL, MYPROW, NB, NN,
     $                   NNB, NPCOL, NPROW, P0, Q0, RSRC, SAVEDT
      REAL               ONE
*     ..
*     .. Local Arrays ..
      INTEGER            DESCIP( DLEN_ ), IWORK( LWORK )
      REAL               WORK( LWORK )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            NUMROC
      EXTERNAL           LSAME, NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           ASSERT, BLACS_GRIDINFO, CHK1MAT, DESCINIT,
     $                   LAIO_INFO, PCHK1MAT, PFSTRSM, PSLAPIV, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, REAL
*     ..
*     .. Executable Statements ..
*
*  Perform error checking.
*
      M = N
      INFO = 0
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYPROW,
     $                     MYPCOL )
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 7*100+CTXT_ )
         RETURN
      ENDIF
      ICTXT = DESCB( CTXT_ )
      CALL BLACS_GRIDINFO( DESCB( CTXT_ ), NPROW, NPCOL, MYPROW,
     $                     MYPCOL )
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 11*100+CTXT_ )
         RETURN
      ENDIF
      ISVALID = ( LSAME( TRANS, 'T' ) .OR. LSAME( TRANS, 'C' ) .OR.
     $          LSAME( TRANS, 'N' ) )
      IF( .NOT.ISVALID ) THEN
         INFO = -1
         CALL PXERBLA( DESCA( CTXT_ ), 'PFSGETRS', 1 )
         RETURN
      ENDIF
      SAVEDT = DESCA( DT_ )
      DESCA( DT_ ) = BLOCK_CYCLIC_2D
      CALL CHK1MAT( N, 2, N, 2, IA, JA, DESCA, 7, INFO )
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PFSGETRS', 7 )
         RETURN
      ENDIF
      IDUMMY = 1
      CALL PCHK1MAT( N, 2, N, 2, IA, JA, DESCA, 7, IDUMMY, IDUMMY,
     $               IDUMMY, INFO )
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PFSGETRS', 7 )
         RETURN
      ENDIF
      DESCA( DT_ ) = SAVEDT
      SAVEDT = DESCB( DT_ )
      DESCB( DT_ ) = BLOCK_CYCLIC_2D
      CALL CHK1MAT( N, 2, NRHS, 3, IB, JB, DESCB, 11, INFO )
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCB( CTXT_ ), 'PFSGETRS', 11 )
         RETURN
      ENDIF
      IDUMMY = 1
      CALL PCHK1MAT( N, 2, NRHS, 3, IB, JB, DESCB, 11, IDUMMY, IDUMMY,
     $               IDUMMY, INFO )
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCB( CTXT_ ), 'PFSGETRS', 11 )
         RETURN
      ENDIF
      DESCB( DT_ ) = SAVEDT
      IF( DESCA( DT_ ).NE.DISK_BLOCK_CYCLIC_2D ) THEN
         INFO = -( 7*100+DT_ )
         CALL PXERBLA( DESCA( CTXT_ ), 'PFSGETRS', 7 )
         RETURN
      ENDIF
      INFO = 0
      CALL LAIO_INFO( DESCA( IODEV_ ), MM, NN, MMB, NNB, MB, NB, CSRC,
     $                RSRC, ICONTXT )
      IF( MB.NE.DESCA( MB_ ) ) THEN
         INFO = -( 7*100+MB_ )
      ENDIF
      IF( NB.NE.DESCA( NB_ ) ) THEN
         INFO = -( 7*100+NB_ )
      ENDIF
      IF( ICONTXT.NE.DESCA( CTXT_ ) ) THEN
         INFO = -( 7*100+CTXT_ )
      ENDIF
      IF( DESCA( RSRC_ ).NE.RSRC ) THEN
         INFO = -( 7*100+RSRC_ )
      ENDIF
      IF( DESCA( CSRC_ ).NE.CSRC ) THEN
         INFO = -( 7*100+CSRC_ )
      ENDIF
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PFSGETRS', 7 )
         RETURN
      ENDIF
      HASWORK = ( N.GE.1 ) .AND. ( NRHS.GE.1 )
      IF( .NOT.HASWORK ) THEN
         RETURN
      ENDIF
      ASIZE = DESCA( SIZE_ )
      ASIZEQUERY = ( ASIZE.EQ.-1 )
*
*     Note the use of mypcol instead of descA(CSRC_).
*     and extra descA(MB_) storage at end.
*
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYPROW,
     $                     MYPCOL )
      MIPIV = DESCA( M_ ) + DESCA( MB_ )*NPROW
      LOCP = NUMROC( DESCA( M_ ), DESCA( MB_ ), MYPROW, DESCA( RSRC_ ),
     $       NPROW )
      LDIPIV = LOCP + DESCA( MB_ )
      CALL DESCINIT( DESCIP, MIPIV, 1, DESCA( MB_ ), 1, DESCA( RSRC_ ),
     $               MYPCOL, DESCA( CTXT_ ), LDIPIV, INFO )
      CALL ASSERT( INFO.EQ.0,
     $             '** PFGETRS: descinit for ipiv returns info ', INFO )
*
*      Need at least one panel.
*
      P0 = MYPROW
      Q0 = MYPCOL
      MM = M
      NN = MIN( NNB, N )
      LOCP = NUMROC( MM, DESCA( MB_ ), MYPROW, P0, NPROW )
      LOCQ = NUMROC( NN, DESCA( NB_ ), MYPCOL, Q0, NPCOL )
      ANEED0 = MAX( 1, LOCP )*MAX( 1, LOCQ )
      MM = N
      NN = MIN( NNB, M )
      LOCP = NUMROC( MM, DESCA( MB_ ), MYPROW, P0, NPROW )
      LOCQ = NUMROC( NN, DESCA( NB_ ), MYPCOL, Q0, NPCOL )
      ANEED1 = MAX( 1, LOCP )*MAX( 1, LOCQ )
      ANEED = MAX( ANEED0, ANEED1 )
      INEED = 0
      IF( ANEED.GT.ASIZE ) THEN
         A( 1 ) = ANEED
         INFO = -( 100*7+SIZE_ )
         IF( .NOT.ASIZEQUERY ) THEN
            CALL PXERBLA( DESCA( CTXT_ ), 'PFSGETRS', 7 )
         ENDIF
         RETURN
      ENDIF
      ONE = REAL( 1 )
      NOTRAN = LSAME( TRANS, 'N' )
      IF( NOTRAN ) THEN
*
*
*         Solve sub( A ) * X = sub( B ).
*
*         Apply row interchanges to the right hand sides.
*
*
         CALL PSLAPIV( 'Forward', 'Row', 'Col', N, NRHS, B, IB, JB,
     $                 DESCB, IPIV, IA, 1, DESCIP, IWORK )
*
*
*         Solve L*X = sub( B ), overwriting sub( B ) with X.
*
*
         CALL PFSTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,
     $                 ONE, A, IA, JA, DESCA, B, IB, JB, DESCB, WORK,
     $                 LWORK, INFO )
         ISVALID = ( INFO.EQ.0 ) .OR. ( ASIZEQUERY .AND.
     $             ( INFO.EQ.-1111 ) )
         CALL ASSERT( ISVALID, '** PFGETRS: PFTRSM returns info ',
     $                INFO )
*
*
*         Solve U*X = sub( B ), overwriting sub( B ) with X.
*
*
         CALL PFSTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
     $                 NRHS, ONE, A, IA, JA, DESCA, B, IB, JB, DESCB,
     $                 WORK, LWORK, INFO )
         ISVALID = ( INFO.EQ.0 ) .OR. ( ASIZEQUERY .AND.
     $             ( INFO.EQ.-1111 ) )
         CALL ASSERT( ISVALID, '** PFGETRS: PFTRSM returns info ',
     $                INFO )
      ELSE
*
*
*         Solve sub( A )' * X = sub( B ).
*
*         Solve U'*X = sub( B ), overwriting sub( B ) with X.
*
*
         CALL PFSTRSM( 'Left', 'Upper', TRANS, 'Non-unit', N, NRHS, ONE,
     $                 A, IA, JA, DESCA, B, IB, JB, DESCB, WORK, LWORK,
     $                 INFO )
         ISVALID = ( INFO.EQ.0 ) .OR. ( ASIZEQUERY .AND.
     $             ( INFO.EQ.-1111 ) )
         CALL ASSERT( ISVALID, '** PFGETRS: PFTRSM returns info ',
     $                INFO )
*
*
*         Solve L'*X = sub( B ), overwriting sub( B ) with X.
*
*
         CALL PFSTRSM( 'Left', 'Lower', TRANS, 'Unit', N, NRHS, ONE, A,
     $                 IA, JA, DESCA, B, IB, JB, DESCB, WORK, LWORK,
     $                 INFO )
         ISVALID = ( INFO.EQ.0 ) .OR. ( ASIZEQUERY .AND.
     $             ( INFO.EQ.-1111 ) )
         CALL ASSERT( ISVALID, '** PFGETRS: PFTRSM returns info ',
     $                INFO )
*
*
*         Apply row interchanges to the solution vectors.
*
*
         CALL PSLAPIV( 'Backward', 'Row', 'Col', N, NRHS, B, IB, JB,
     $                 DESCB, IPIV, IA, 1, DESCIP, IWORK )
      ENDIF
      RETURN
      END
