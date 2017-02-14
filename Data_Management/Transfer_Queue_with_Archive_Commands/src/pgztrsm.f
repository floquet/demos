      SUBROUTINE PGZTRSM( SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A, IA,
     $                    JA, DESCA, B, IB, JB, DESCB )
*
*
*  -- ScaLAPACK auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     Oct 10, 1996
*
*
*
*  Purpose
*  =======
*
*  PGTRSM  solves one of the distributed matrix equations
*
*                 op( sub( A ) )*X = alpha*sub( B ),   or
*
*                 X*op( sub( A ) ) = alpha*sub( B ),
*
*  where sub( A ) denotes A(IA:IA+M-1,JA:JA+M-1)  if SIDE = 'L',
*        sub( A ) denotes A(IA:IA+N-1,JA:JA+N-1)  if SIDE = 'R',
*
*        sub( B ) denotes B(IB:IB+M-1,JB:JB+N-1),
*
*  alpha is a scalar, X and sub( B ) are an M-by-N distributed matrix,
*  sub( A ) is a unit, or non-unit, upper or lower triangular distribu-
*  ted matrix and op( A ) is one of
*
*     op( A ) = A   or   op( A ) = A'   or   op( A ) = conjg( A' ).
*
*  The distributed matrix X is overwritten on sub( B ).
*
*
*  This implementation attempts to overcome some alignment
*  constraints in PGRSM.
*
*
* Note:
* ====
*
* sub(C) is m x m and B is m x n is side = 'L'
* sub(C) is n x n and B is m x n is side = 'R'
*
* sub(C) is in-core matrix.
*
*
*
*
*
* Several cases to consider:
*
* forward or backward solve correspond to
* how the right-hand side B(*) is sweeped.
*
*
* L x = b        forward solve  axpyupdate   notransA
* L^t x = b      backward solve  dotupdate    transA
* U x = b        backward solve  axpyupdate   notransA
* U^t x = b      forward solve   dotupdate    transA
*
* x L = b        backward solve  dotupdate    notransA
* x L^t = b      forward solve   axpyupdate   transA
* x U = b        forward solve   dotupdate    notransA
* x U^t = b      backward solve  axpyupdate   transA
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DT_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1 )
      INTEGER            CTXT_, MB_, NB_
      PARAMETER          ( CTXT_ = 2, MB_ = 5, NB_ = 6 )
      INTEGER            RSRC_, CSRC_, LLD_
      PARAMETER          ( RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      INTEGER            LWORK
      PARAMETER          ( LWORK = 100*100 )
      LOGICAL            USE_PTRSM
      PARAMETER          ( USE_PTRSM = .false. )
*     ..
*     .. Scalar Arguments ..
      CHARACTER          DIAG, SIDE, TRANS, UPLO
      INTEGER            IA, IB, JA, JB, M, N
      COMPLEX*16         ALPHA
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( DLEN_ ), DESCB( DLEN_ )
      COMPLEX*16         A( * ), B( * )
*     ..
*     .. Local Scalars ..
      LOGICAL            HASWORK, ISAXPYUPDATE, ISBACKWARD, ISDOTUPDATE,
     $                   ISFORWARD, ISLEFT, ISLOWER, ISRIGHT, ISTRANS,
     $                   ISTRANSA, ISTRSM_ALIGNED, ISUPPER, ISVALID,
     $                   NOTRANS, NOTRANSA, USEATMP
      CHARACTER          TRANSA
      INTEGER            CDEST, CPROC, IADIAG, IAPOS, IBPOS, IBSIZE, IC,
     $                   ICOLA, ICOLB, ICONTXT, ICPOS, ICTXT, IDUMMY,
     $                   IIBPOS, IJOFF, INC, INEED, INFO, IP_APOS,
     $                   IP_ATMP, IROWA, IROWB, ISRC, ITMP, JADIAG,
     $                   JAPOS, JBPOS, JC, JCPOS, JEND, JENDLR, JINC,
     $                   JJBPOS, JSIZE, JSTART, KK, LCINDX, LCOFF, LDA,
     $                   LDATMP, LDB, LRINDX, LROFF, MM, MYID, MYPCOL,
     $                   MYPROW, NN, NPCOL, NPROC, NPROW, RDEST, RPROC,
     $                   SAVEDT
      COMPLEX*16         LALPHA, LBETA, ONE, ZERO
*     ..
*     .. Local Arrays ..
      INTEGER            DESCATMP( DLEN_ )
      COMPLEX*16         WORK( LWORK )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            INDXG2P
      EXTERNAL           LSAME, INDXG2P
*     ..
*     .. External Subroutines ..
      EXTERNAL           ASSERT, BLACS_GRIDINFO, BLACS_PINFO, CHK1MAT,
     $                   DESCINIT, INFOG2L, PCHK1MAT, PXERBLA, PZGEMM,
     $                   PZSCAL, PZTRSM, ZFILL, ZGERV2D, ZGESD2D
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, INT, MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
* ===========================================================
*
*  Check parameters.
*
      ONE = DCMPLX( DBLE( 1 ) )
      ZERO = DCMPLX( DBLE( 0 ) )
      IC = IA
      JC = JA
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
         CALL PXERBLA( ICONTXT, 'PGZTRSM', 1 )
         RETURN
      ENDIF
      ISVALID = ( LSAME( TRANS, 'T' ) .OR. LSAME( TRANS, 'C' ) .OR.
     $          LSAME( TRANS, 'N' ) )
      IF( .NOT.ISVALID ) THEN
         INFO = -3
         CALL PXERBLA( ICONTXT, 'PGZTRSM', 3 )
         RETURN
      ENDIF
      ISVALID = ( LSAME( UPLO, 'U' ) .OR. LSAME( UPLO, 'L' ) )
      IF( .NOT.ISVALID ) THEN
         INFO = -2
         CALL PXERBLA( ICONTXT, 'PGZTRSM', 2 )
         RETURN
      ENDIF
      ISVALID = ( LSAME( DIAG, 'U' ) .OR. LSAME( DIAG, 'N' ) )
      IF( .NOT.ISVALID ) THEN
         INFO = -4
         CALL PXERBLA( ICONTXT, 'PGZTRSM', 4 )
         RETURN
      ENDIF
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
            CALL PXERBLA( DESCA( CTXT_ ), 'PGZTRSM', 11 )
            RETURN
         ENDIF
         IDUMMY = 1
         CALL PCHK1MAT( M, 5, M, 5, IA, JA, DESCA, 11, IDUMMY, IDUMMY,
     $                  IDUMMY, INFO )
         IF( INFO.NE.0 ) THEN
            CALL PXERBLA( DESCA( CTXT_ ), 'PGZTRSM', 11 )
            RETURN
         ENDIF
         DESCA( DT_ ) = SAVEDT
      ELSE
         SAVEDT = DESCA( DT_ )
         DESCA( DT_ ) = BLOCK_CYCLIC_2D
         CALL CHK1MAT( N, 6, N, 6, IA, JA, DESCA, 11, INFO )
         IF( INFO.NE.0 ) THEN
            CALL PXERBLA( DESCA( CTXT_ ), 'PGZTRSM', 11 )
            RETURN
         ENDIF
         IDUMMY = 1
         CALL PCHK1MAT( N, 6, N, 6, IA, JA, DESCA, 11, IDUMMY, IDUMMY,
     $                  IDUMMY, INFO )
         IF( INFO.NE.0 ) THEN
            CALL PXERBLA( DESCA( CTXT_ ), 'PGZTRSM', 11 )
            RETURN
         ENDIF
         DESCA( DT_ ) = SAVEDT
      ENDIF
      SAVEDT = DESCB( DT_ )
      DESCB( DT_ ) = BLOCK_CYCLIC_2D
      CALL CHK1MAT( M, 5, N, 6, IB, JB, DESCB, 15, INFO )
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCB( CTXT_ ), 'PGZTRSM', 15 )
         RETURN
      ENDIF
      IDUMMY = 1
      CALL PCHK1MAT( M, 5, N, 6, IB, JB, DESCB, 15, IDUMMY, IDUMMY,
     $               IDUMMY, INFO )
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCB( CTXT_ ), 'PGZTRSM', 15 )
         RETURN
      ENDIF
      DESCB( DT_ ) = SAVEDT
*
* Check for quick return;
*
      HASWORK = ( M.GE.1 ) .AND. ( N.GE.1 )
      IF( .NOT.HASWORK ) THEN
         WORK( 1 ) = 0
         RETURN
      ENDIF
      CALL BLACS_PINFO( MYID, NPROC )
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYPROW,
     $                     MYPCOL )
      LDA = DESCA( LLD_ )
      LDB = DESCB( LLD_ )
      NOTRANS = LSAME( TRANS, 'N' )
      ISTRANS = ( .NOT.NOTRANS )
      ISLOWER = LSAME( UPLO, 'L' )
      ISUPPER = LSAME( UPLO, 'U' )
      ISLEFT = LSAME( SIDE, 'L' )
      ISRIGHT = LSAME( SIDE, 'R' )
      ISDOTUPDATE = ( ISLEFT .AND. ( ( ISLOWER .AND. ISTRANS ) .OR.
     $              ( ISUPPER .AND. ISTRANS ) ) ) .OR.
     $              ( ISRIGHT .AND. ( ( ISLOWER .AND. NOTRANS ) .OR.
     $              ( ISUPPER .AND. NOTRANS ) ) )
      ISTRANSA = ISTRANS
      ISFORWARD = ( ISLEFT .AND. ( ( ISLOWER .AND. NOTRANS ) .OR.
     $            ( ISUPPER .AND. ISTRANS ) ) ) .OR.
     $            ( ISRIGHT .AND. ( ( ISLOWER .AND. ISTRANS ) .OR.
     $            ( ISUPPER .AND. NOTRANS ) ) )
      ISAXPYUPDATE = ( .NOT.ISDOTUPDATE )
      ISBACKWARD = ( .NOT.ISFORWARD )
      NOTRANSA = ( .NOT.ISTRANSA )
      IF( NOTRANSA ) THEN
         TRANSA = 'N'
      ELSE
         TRANSA = 'C'
      ENDIF
      JINC = MAX( DESCA( MB_ ), DESCA( NB_ ) )
      CALL ASSERT( JINC.GE.1,
     $             '**  PGTRSM: need descA(N_) to be at least 1 ',
     $             JINC )
      IP_ATMP = 1
*
*  Calculate worst case requirement for Atmp.
*
      INEED = JINC*JINC
      CALL ASSERT( LWORK.GE.INEED, '**  PGTRSM: increase lwork to ',
     $             INEED )
*
* Check if PTRSM can handle it.
*
      IROWB = INDXG2P( IB, DESCB( MB_ ), MYPROW, DESCB( RSRC_ ), NPROW )
      ICOLB = INDXG2P( JB, DESCB( NB_ ), MYPCOL, DESCB( CSRC_ ), NPCOL )
      IROWA = INDXG2P( IA, DESCA( MB_ ), MYPROW, DESCA( RSRC_ ), NPROW )
      ICOLA = INDXG2P( JA, DESCA( NB_ ), MYPCOL, DESCA( CSRC_ ), NPCOL )
      ISTRSM_ALIGNED = ( ( LSAME( SIDE, 'L' ) ) .AND.
     $                 ( IROWA.EQ.IROWB ) ) .OR.
     $                 ( ( LSAME( SIDE, 'R' ) ) .AND.
     $                 ( ICOLA.EQ.ICOLB ) )
      IF( ( USE_PTRSM ) .AND. ( ISTRSM_ALIGNED ) ) THEN
         CALL PZTRSM( SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A, IA, JA,
     $                DESCA, B, IB, JB, DESCB )
         WORK( 1 ) = 0
         RETURN
*
*
* =====================================
* Basic algorithm:
*
* Use diagonal block to solve part of B.
* Update the rest of B.
* =====================================
*
*
      ENDIF
      IF( ISLEFT ) THEN
         JENDLR = M
         IJOFF = ( IC-1 )
      ELSE
         JENDLR = N
         IJOFF = ( JC-1 )
      ENDIF
      JSTART = 1
      JEND = JENDLR
   10 CONTINUE
      IF( .NOT.( ( ISFORWARD .AND. ( JSTART.GT.JENDLR ) ) .OR.
     $    ( ISBACKWARD .AND. ( JEND.LT.1 ) ) ) ) THEN
*
*
*  Operate on B( jstart:jend,1:n)       for side = 'L'
*  Operate on B( 1:m, jstart:jend)      for side = 'R'
*
*
         IF( ISFORWARD ) THEN
            JEND = MIN( JENDLR, JSTART+JINC-1 )
*
*               Attempt to be block aligned.
*
            IF( JEND.NE.JENDLR ) THEN
               ITMP = IJOFF + JEND
               ITMP = MAX( IJOFF+JSTART, INT( ITMP / JINC )*JINC )
               JEND = MAX( JSTART, ITMP-IJOFF )
            ENDIF
         ELSE
            JSTART = MAX( 1, JEND-JINC+1 )
*
*               Attempt to be block aligned.
*
            IF( JSTART.NE.1 ) THEN
               ITMP = IJOFF + JSTART
   20          CONTINUE
               IF( MOD( ITMP, JINC ).NE.1 ) THEN
                  ITMP = ITMP + 1
                  GOTO 20
               ENDIF
   30          CONTINUE
               JSTART = MAX( 1, MIN( JEND, ITMP-IJOFF ) )
            ENDIF
         ENDIF
         JSIZE = JEND - JSTART + 1
*
*       Position of diagonal.
*
         ICPOS = ( IC-1 ) + JSTART
         JCPOS = ( JC-1 ) + JSTART
         IADIAG = ICPOS
         JADIAG = JCPOS
         LROFF = MOD( IADIAG, JINC )
         IF( LROFF.EQ.0 ) THEN
            LROFF = JINC
         ENDIF
         LCOFF = MOD( JADIAG, JINC )
         IF( LCOFF.EQ.0 ) THEN
            LCOFF = JINC
         ENDIF
         IF( ISDOTUPDATE ) THEN
*
*            Perform dot product like computation
*
            IF( ISFORWARD ) THEN
*
*               Reference previously computed
*               solution  B(1:jstart-1)
*
               IBPOS = ( IB-1 ) + 1
               JBPOS = ( JB-1 ) + 1
               IBSIZE = JSTART - 1
            ELSE
               IF( ISBACKWARD ) THEN
*
*               Reference previously computed
*               solution  B( (jend+1):jendLR )
*
                  IF( ISLEFT ) THEN
                     IBPOS = ( IB-1 ) + ( JEND+1 )
                     JBPOS = ( JB-1 ) + 1
                  ELSE
                     IBPOS = ( IB-1 ) + 1
                     JBPOS = ( JB-1 ) + ( JEND+1 )
                  ENDIF
                  IBSIZE = JENDLR - ( JEND+1 ) + 1
               ENDIF
            ENDIF
            IF( ISLOWER ) THEN
               IAPOS = ( IC-1 ) + ( JEND+1 )
               JAPOS = ( JC-1 ) + JSTART
            ELSE
               IAPOS = ( IC-1 ) + 1
               JAPOS = ( JC-1 ) + JSTART
            ENDIF
*
*               Reminder:
*
*               PGEMM: C <- lalpha*op(A)*op(B) + lbeta*C
*
*               C is mm x nn
*
            LALPHA = -ONE
            LBETA = ONE
            IF( ISLEFT ) THEN
               ICPOS = ( IB-1 ) + JSTART
               JCPOS = ( JB-1 ) + 1
               MM = JSIZE
               NN = N
               KK = IBSIZE
            ELSE
               ICPOS = ( IB-1 ) + 1
               JCPOS = ( JB-1 ) + JSTART
               MM = M
               NN = JSIZE
               KK = IBSIZE
            ENDIF
            HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 ) .AND. ( KK.GE.1 )
            IF( HASWORK ) THEN
               IF( ISLEFT ) THEN
                  CALL PZGEMM( TRANSA, 'N', MM, NN, KK, LALPHA, A,
     $                         IAPOS, JAPOS, DESCA, B, IBPOS, JBPOS,
     $                         DESCB, LBETA, B, ICPOS, JCPOS, DESCB )
               ELSE
                  CALL PZGEMM( 'N', TRANSA, MM, NN, KK, LALPHA, B,
     $                         IBPOS, JBPOS, DESCB, A, IAPOS, JAPOS,
     $                         DESCA, LBETA, B, ICPOS, JCPOS, DESCB )
               ENDIF
            ENDIF
         ENDIF
*
*--------------------------------------------------------------
*    Apply inverse(A(jstart:jend,jstart:jend)) to
*               B(jstart:jend,*) for side= 'L'
*    Apply inverse(A(jstart:jend,jstart:jend)) to
*               B(*,jstart:jend) for side= 'R'
*--------------------------------------------------------------
*
         IF( ISLEFT ) THEN
            IBPOS = ( IB-1 ) + JSTART
            JBPOS = ( JB-1 ) + 1
            MM = JSIZE
            NN = N
         ELSE
            IBPOS = ( IB-1 ) + 1
            JBPOS = ( JB-1 ) + JSTART
            MM = M
            NN = JSIZE
         ENDIF
         LALPHA = ONE
         HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
         CALL ASSERT( HASWORK,
     $                '** PGTRSM: internal error with PTRSM, mm ', MM )
         IROWB = INDXG2P( IBPOS, DESCB( MB_ ), MYPROW, DESCB( RSRC_ ),
     $           NPROW )
         ICOLB = INDXG2P( JBPOS, DESCB( NB_ ), MYPCOL, DESCB( CSRC_ ),
     $           NPCOL )
         IROWA = INDXG2P( IADIAG, DESCA( MB_ ), MYPROW, DESCA( RSRC_ ),
     $           NPROW )
         ICOLA = INDXG2P( JADIAG, DESCA( NB_ ), MYPCOL, DESCA( CSRC_ ),
     $           NPCOL )
         ISTRSM_ALIGNED = ( ( LSAME( SIDE, 'L' ) ) .AND.
     $                    ( IROWA.EQ.IROWB ) ) .OR.
     $                    ( ( LSAME( SIDE, 'R' ) ) .AND.
     $                    ( ICOLA.EQ.ICOLB ) )
         USEATMP = .NOT.( USE_PTRSM .AND. ISTRSM_ALIGNED )
         IF( USEATMP ) THEN
*
*           Copy part of matrix into  a temporary array
*           to get around alignment limitations in PTSRM.
*
            LDATMP = JINC
            CALL ZFILL( JINC*JINC, ZERO, WORK( IP_ATMP ), 1 )
*
*               Create a fake descriptor for this scalapack array.
*               note Atmp(1,1) is owned by the same processor
*               holding B(ibpos,jbpos).
*
            CALL DESCINIT( DESCATMP, JSIZE, JSIZE, DESCA( MB_ ),
     $                     DESCA( NB_ ), IROWB, ICOLB, DESCB( CTXT_ ),
     $                     MAX( 1, LDATMP ), INFO )
            CALL ASSERT( INFO.EQ.0, '**  PGTRSM: descinit return info ',
     $                   INFO )
            IF( ( MYPROW.EQ.IROWA ) .AND. ( MYPCOL.EQ.ICOLA ) ) THEN
* Package and send
               CALL INFOG2L( IADIAG, JADIAG, DESCA, NPROW, NPCOL,
     $                       MYPROW, MYPCOL, LRINDX, LCINDX, RPROC,
     $                       CPROC )
               CALL ASSERT( ( RPROC.EQ.MYPROW ) .AND.
     $                      ( CPROC.EQ.MYPCOL ),
     $             '**  PGTRSM: internal error in calling infog2l,myid '
     $                      , MYID )
               ISRC = LRINDX + ( LCINDX-1 )*DESCA( LLD_ )
               ISVALID = ( 1.LE.ISRC )
               CALL ASSERT( ISVALID,
     $                  '** PGTRSM: internal error invalid isrc, isrc= '
     $                      , ISRC )
               RDEST = IROWB
               CDEST = ICOLB
               CALL ZGESD2D( DESCA( CTXT_ ), JSIZE, JSIZE, A( ISRC ),
     $                       DESCA( LLD_ ), RDEST, CDEST )
            ENDIF
            IF( ( MYPROW.EQ.IROWB ) .AND. ( MYPCOL.EQ.ICOLB ) ) THEN
               IP_APOS = IP_ATMP + ( LROFF-1 ) + ( LCOFF-1 )*LDATMP
               CALL ZGERV2D( DESCA( CTXT_ ), JSIZE, JSIZE,
     $                       WORK( IP_APOS ), LDATMP, IROWA, ICOLA )
            ENDIF
            CALL ASSERT( MOD( LROFF, JINC ).EQ.MOD( IBPOS, JINC ),
     $                   '** PGTRSM:  invalid lroff = ', LROFF )
            CALL ASSERT( MOD( LCOFF, JINC ).EQ.MOD( JBPOS, JINC ),
     $                   '** PGTRSM: invalid lcoff = ', LCOFF )
            CALL PZTRSM( SIDE, UPLO, TRANS, DIAG, MM, NN, LALPHA,
     $                   WORK( IP_ATMP ), LROFF, LCOFF, DESCATMP, B,
     $                   IBPOS, JBPOS, DESCB )
         ELSE
*
*               Matrices are TRSM aligned.
*               No need to perform copy.
*
            CALL PZTRSM( SIDE, UPLO, TRANS, DIAG, MM, NN, LALPHA, A,
     $                   IADIAG, JADIAG, DESCA, B, IBPOS, JBPOS, DESCB )
         ENDIF
* ------------------------------------------------------
* use the newly computed X(:) to update the rest of the
* entries in B
* ------------------------------------------------------
         IF( ISAXPYUPDATE ) THEN
*
*           Column update in B
*
            IF( ISFORWARD ) THEN
*
*               The rest of the vector (jend+1):jendLR
*
               IF( ISLEFT ) THEN
                  IBPOS = ( IB-1 ) + ( JEND+1 )
                  JBPOS = ( JB-1 ) + 1
               ELSE
                  IBPOS = ( IB-1 ) + 1
                  JBPOS = ( JB-1 ) + ( JEND+1 )
               ENDIF
               IBSIZE = JENDLR - ( JEND+1 ) + 1
            ELSE
*
*               The rest of the vector 1:(jstart-1)
*
               IBPOS = ( IB-1 ) + 1
               JBPOS = ( JB-1 ) + 1
               IBSIZE = JSTART - 1
            ENDIF
            IF( ISLOWER ) THEN
               IAPOS = ( IC-1 ) + JEND + 1
               JAPOS = ( JC-1 ) + JSTART
            ELSE
               IF( ISUPPER ) THEN
                  IAPOS = ( IC-1 ) + 1
                  JAPOS = ( JC-1 ) + JSTART
               ENDIF
            ENDIF
            LALPHA = -ONE
            LBETA = ONE
            IF( ISLEFT ) THEN
               MM = IBSIZE
               NN = N
               KK = JSIZE
               IIBPOS = ( IB-1 ) + JSTART
               JJBPOS = ( JB-1 ) + 1
            ELSE
               MM = M
               NN = IBSIZE
               KK = JSIZE
               IIBPOS = ( IB-1 ) + 1
               JJBPOS = ( JB-1 ) + JSTART
            ENDIF
            HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 ) .AND. ( KK.GE.1 )
            IF( HASWORK ) THEN
               IF( ISLEFT ) THEN
                  CALL PZGEMM( TRANSA, 'N', MM, NN, KK, LALPHA, A,
     $                         IAPOS, JAPOS, DESCA, B, IIBPOS, JJBPOS,
     $                         DESCB, LBETA, B, IBPOS, JBPOS, DESCB )
               ELSE
                  CALL PZGEMM( 'N', TRANSA, MM, NN, KK, LALPHA, B,
     $                         IIBPOS, JJBPOS, DESCB, A, IAPOS, JAPOS,
     $                         DESCA, LBETA, B, IBPOS, JBPOS, DESCB )
               ENDIF
            ENDIF
         ENDIF
* end if (isaxpyupdate)
* prepare for next iteration
         IF( ISFORWARD ) THEN
            JSTART = JEND + 1
         ELSE
            JEND = JSTART - 1
         ENDIF
         GOTO 10
      ENDIF
   40 CONTINUE
* end while
* -----------------------------------
*  need to rescale the solution by alpha
* -----------------------------------
      IF( ALPHA.NE.ONE ) THEN
         DO 50 JSTART = 1, N
            IBPOS = ( IB-1 ) + 1
            JBPOS = ( JB-1 ) + JSTART
            INC = 1
            CALL PZSCAL( M, ALPHA, B, IBPOS, JBPOS, DESCB, INC )
   50    CONTINUE
   60    CONTINUE
      ENDIF
      WORK( 1 ) = INEED
      RETURN
      END
