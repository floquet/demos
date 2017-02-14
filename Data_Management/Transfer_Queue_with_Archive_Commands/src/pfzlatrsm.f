      SUBROUTINE PFZLATRSM( SIDE, UPLO, TRANS, DIAG, IODEV, M, N, IC,
     $                      JC, ALPHA, B, IB, JB, DESCB, WORK, LWORK,
     $                      INFO )
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
* Triangular solve
*
* op( sub(C) ) X = alpha sub(B)
*
* Note:
*
* sub(C) is m x m and B is m x n is side = 'L'
* sub(C) is n x n and B is m x n is side = 'R'
*
* sub(C) is out-of-core matrix
*
*
*
* Several cases to consider:
*
* forward or backward solve correspond to
* how the right-hand side B(*) is sweeped.
*
* L x = b        forward solve   axpyupdate   notransA
* L^t x = b      backward solve  dotupdate    transA
* U x = b        backward solve  axpyupdate   notransA
* U^t x = b      forward solve   dotupdate    transA
*
* x L = b        backward solve  dotupdate    notransA
* x L^t = b      forward solve   axpyupdate   transA
* x U = b        forward solve   dotupdate    notransA
* x U^t = b      backward solve  axpyupdate   transA
*
*
*
*     .. Parameters ..
      INTEGER            DLEN_
      PARAMETER          ( DLEN_ = 9 )
      INTEGER            CTXT_, MB_, NB_
      PARAMETER          ( CTXT_ = 2, MB_ = 5, NB_ = 6 )
      INTEGER            RSRC_, CSRC_
      PARAMETER          ( RSRC_ = 7, CSRC_ = 8 )
*     ..
*     .. Scalar Arguments ..
      CHARACTER          DIAG, SIDE, TRANS, UPLO
      INTEGER            IB, IC, INFO, IODEV, JB, JC, LWORK, M, N
      COMPLEX*16         ALPHA
*     ..
*     .. Array Arguments ..
      INTEGER            DESCB( DLEN_ )
      COMPLEX*16         B( * ), WORK( * )
*     ..
*     .. Local Scalars ..
      LOGICAL            HASWORK, ISAXPYUPDATE, ISBACKWARD, ISDOTUPDATE,
     $                   ISFORWARD, ISLEFT, ISLOWER, ISRIGHT, ISTRANS,
     $                   ISTRANSA, ISUPPER, ISVALID, LWORKQUERY,
     $                   NOTRANS, NOTRANSA
      CHARACTER          TRANSA
      INTEGER            CSRC, I, IA, IADIAG, IBPOS, IBSIZE, ICOLB,
     $                   ICOLC, ICONTXT, ICPOS, IIBPOS, IIC, INC, INEED,
     $                   IP_AMAT, IROWB, IROWC, ISIZE, JA, JADIAG,
     $                   JBPOS, JCPOS, JEND, JENDLR, JINC, JJBPOS, JJC,
     $                   JSIZE, JSTART, KK, LDA, LOCP, LOCQ, MB, MM,
     $                   MMB, MYID, MYPCOL, MYPROW, M_C, NB, NN, NNB,
     $                   NPCOL, NPROC, NPROW, N_C, P0, Q0, RSRC
      COMPLEX*16         LALPHA, LBETA, ONE, ZERO
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA( DLEN_ )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            INDXG2P, NUMROC
      EXTERNAL           LSAME, INDXG2P, NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           ASSERT, BLACS_GRIDINFO, BLACS_PINFO, DESCINIT,
     $                   LAIO_INFO, PFMAXSIZE, PXERBLA, PZGEMM, PZSCAL,
     $                   PZTRSM, ZLAREAD
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, ICHAR, MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
* ===========================================================
      ONE = DCMPLX( DBLE( 1 ) )
      ZERO = DCMPLX( DBLE( 0 ) )
      CALL BLACS_PINFO( MYID, NPROC )
      CALL LAIO_INFO( IODEV, MM, NN, MMB, NNB, MB, NB, CSRC, RSRC,
     $                ICONTXT )
      CALL BLACS_GRIDINFO( ICONTXT, NPROW, NPCOL, MYPROW, MYPCOL )
      IP_AMAT = 1
      NOTRANS = LSAME( TRANS, 'N' )
      ISTRANS = ( .NOT.NOTRANS )
      ISVALID = ( ISTRANS .OR. NOTRANS )
      IF( .NOT.ISVALID ) THEN
         WRITE( *, FMT = * )'** PFLATRSM: invalid value for trans ',
     $      TRANS
         STOP '** error in PFLATRSM ** '
      ENDIF
      ISLOWER = LSAME( UPLO, 'L' )
      ISUPPER = LSAME( UPLO, 'U' )
      ISVALID = ( ISLOWER .OR. ISUPPER )
      CALL ASSERT( ISVALID, '** PFLATRSM: invalid uplo, ichar(uplo) =  '
     $             , ICHAR( UPLO ) )
      ISLEFT = LSAME( SIDE, 'L' )
      ISRIGHT = LSAME( SIDE, 'R' )
      ISVALID = ISLEFT .OR. ISRIGHT
      CALL ASSERT( ISVALID, '** PFLATRSM: invalid side, ichar(side) =  '
     $             , ICHAR( SIDE ) )
      ISDOTUPDATE = ( ISLEFT .AND. ( ( ISLOWER .AND. ISTRANS ) .OR.
     $              ( ISUPPER .AND. ISTRANS ) ) ) .OR.
     $              ( ISRIGHT .AND. ( ( ISLOWER .AND. NOTRANS ) .OR.
     $              ( ISUPPER .AND. NOTRANS ) ) )
      ISTRANSA = ( ISLEFT .AND. ( ( ISLOWER .AND. ISTRANS ) .OR.
     $           ( ISUPPER .AND. ISTRANS ) ) ) .OR.
     $           ( ISRIGHT .AND. ( ( ISLOWER .AND. ISTRANS ) .OR.
     $           ( ISUPPER .AND. ISTRANS ) ) )
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
      IF( ISLEFT ) THEN
*
*       sub(C) is m x m
*
         M_C = M
         N_C = M
         JENDLR = M
         JINC = MIN( M, NNB )
      ELSE
*
*       sub(C) is n x n
*
         M_C = N
         N_C = N
         JENDLR = N
         JINC = MIN( N, NNB )
      ENDIF
*
*       Check whether there is sufficient work space.
*       Note the use of p0,q0 to overestimate storage.
*
      P0 = MYPROW
      Q0 = MYPCOL
      LOCP = NUMROC( M_C, MB, MYPROW, P0, NPROW )
      LOCQ = NUMROC( JINC, NB, MYPCOL, Q0, NPCOL )
      INEED = LOCP*LOCQ
      LWORKQUERY = ( LWORK.EQ.-1 )
      IF( LWORK.LT.INEED ) THEN
         WORK( 1 ) = INEED
         INFO = -16
         IF( .NOT.LWORKQUERY ) THEN
            CALL PXERBLA( DESCB( CTXT_ ), 'PFZLATRSM', 16 )
         ENDIF
         RETURN
      ENDIF
      CALL PFMAXSIZE( 'C', LWORK, M_C, JINC, MB, NB, P0, Q0, ICONTXT,
     $                INFO )
      CALL ASSERT( JINC.GE.MIN( N_C, NNB ),
     $             '** PFLATRSM: pfmaxsize returns info = ', INFO )
      JINC = MIN( JINC, N_C )
      IF( JINC.NE.N_C ) THEN
         IF( JINC.GE.NNB ) THEN
            JINC = JINC - MOD( JINC, NNB )
         ENDIF
      ENDIF
*
*
* =====================================
* Basic algorithm:
*
* Bring in a part of out of core matrix.
* Use diagonal block to solve part of B.
* Update the rest of B.
* =====================================
*
*
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
            IF( JEND.NE.JENDLR ) THEN
*
*               Try to be block aligned.
*
               IF( ISUPPER ) THEN
                  JJC = ( JC-1 ) + JEND
   20             CONTINUE
                  IF( MOD( JJC, NB ).NE.0 ) THEN
                     JJC = JJC - 1
                     GOTO 20
                  ENDIF
   30             CONTINUE
                  JEND = JJC - ( JC-1 )
               ELSE
                  IIC = ( IC-1 ) + JEND
   40             CONTINUE
                  IF( MOD( IIC, MB ).NE.0 ) THEN
                     IIC = IIC - 1
                     GOTO 40
                  ENDIF
   50             CONTINUE
                  JEND = IIC - ( IC-1 )
               ENDIF
               JEND = MIN( JENDLR, MAX( JEND, JSTART ) )
            ENDIF
         ELSE
            JSTART = MAX( 1, JEND-JINC+1 )
            IF( JSTART.NE.1 ) THEN
*
*               Try to be block aligned.
*
               IF( ISUPPER ) THEN
                  JJC = ( JC-1 ) + JSTART
   60             CONTINUE
                  IF( MOD( JJC, NB ).NE.1 ) THEN
                     JJC = JJC + 1
                     GOTO 60
                  ENDIF
   70             CONTINUE
                  JSTART = JJC - ( JC-1 )
               ELSE
                  IIC = ( IC-1 ) + JSTART
   80             CONTINUE
                  IF( MOD( IIC, MB ).NE.1 ) THEN
                     IIC = IIC + 1
                     GOTO 80
                  ENDIF
   90             CONTINUE
                  JSTART = IIC - ( IC-1 )
               ENDIF
               JSTART = MAX( 1, MIN( JSTART, JEND ) )
            ENDIF
         ENDIF
         JSIZE = JEND - JSTART + 1
*      ---------------------------------------------------
*      bring in a piece from out-of-core matrix
*      ---------------------------------------------------
         IF( ISLOWER ) THEN
            ICPOS = ( IC-1 ) + JSTART
            JCPOS = ( JC-1 ) + JSTART
            ISIZE = JENDLR - JSTART + 1
* jstart:m
            IA = 1
            JA = 1
            IADIAG = 1
            JADIAG = 1
         ELSE
            IF( ISUPPER ) THEN
               ICPOS = ( IC-1 ) + 1
               JCPOS = ( JC-1 ) + JSTART
               ISIZE = JEND
* 1:jend
               IA = 1
               JA = 1
               IADIAG = JSTART
               JADIAG = 1
            ENDIF
         ENDIF
         IROWC = INDXG2P( ICPOS, MB, MYPROW, RSRC, NPROW )
         ICOLC = INDXG2P( JCPOS, NB, MYPCOL, CSRC, NPCOL )
         P0 = IROWC
         Q0 = ICOLC
         IF( ISLEFT ) THEN
            IBPOS = ( IB-1 ) + JSTART
            JBPOS = ( JB-1 ) + 1
         ELSE
            IBPOS = ( IB-1 ) + 1
            JBPOS = ( JB-1 ) + JSTART
         ENDIF
         IROWB = INDXG2P( IBPOS, DESCB( MB_ ), MYPROW, DESCB( RSRC_ ),
     $           NPROW )
         ICOLB = INDXG2P( JBPOS, DESCB( NB_ ), MYPROW, DESCB( CSRC_ ),
     $           NPCOL )
         IF( ISLEFT ) THEN
            P0 = IROWB
         ELSE
            Q0 = ICOLB
         ENDIF
         LOCP = NUMROC( ISIZE, MB, MYPROW, P0, NPROW )
         LOCQ = NUMROC( JSIZE, NB, MYPCOL, Q0, NPCOL )
         IF( .NOT.LWORKQUERY ) THEN
            CALL ASSERT( LOCP*LOCQ.LE.LWORK,
     $                   '** PFLATRSM: insufficient lwork ', LWORK )
         ENDIF
         LDA = MAX( 1, LOCP )
         CALL DESCINIT( DESCA, ISIZE, JSIZE, MB, NB, P0, Q0, ICONTXT,
     $                  LDA, INFO )
         INFO = 0
         MM = ISIZE
         NN = JSIZE
         CALL ZLAREAD( IODEV, ISIZE, JSIZE, ICPOS, JCPOS,
     $                 WORK( IP_AMAT ), IA, JA, DESCA, INFO )
         CALL ASSERT( INFO.EQ.0, '** PFLATRSM: LAREAD returns info ',
     $                INFO )
         IF( ISDOTUPDATE ) THEN
* perform dot product like computation
            IF( ISFORWARD ) THEN
* reference previously computed solution
               IBPOS = ( IB-1 ) + 1
               JBPOS = ( JB-1 ) + 1
               IBSIZE = JSTART - 1
            ELSE
               IF( ISBACKWARD ) THEN
* reference previously computed solution
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
               IA = JSIZE + 1
               JA = 1
            ELSE
               IA = 1
               JA = 1
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
                  CALL PZGEMM( TRANSA, 'N', MM, NN, KK, LALPHA,
     $                         WORK( IP_AMAT ), IA, JA, DESCA, B, IBPOS,
     $                         JBPOS, DESCB, LBETA, B, ICPOS, JCPOS,
     $                         DESCB )
               ELSE
                  CALL PZGEMM( 'N', TRANSA, MM, NN, KK, LALPHA, B,
     $                         IBPOS, JBPOS, DESCB, WORK( IP_AMAT ), IA,
     $                         JA, DESCA, LBETA, B, ICPOS, JCPOS,
     $                         DESCB )
               ENDIF
            ENDIF
         ENDIF
*
* --------------------------------------------------------------
* Apply inverse( C(jstart:jend, jstart:jend) ) to B(jstart:jend)
* --------------------------------------------------------------
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
     $                '** PFLATRSM: internal error with PTRSM, mm ',
     $                MM )
         CALL PZTRSM( SIDE, UPLO, TRANS, DIAG, MM, NN, LALPHA,
     $                WORK( IP_AMAT ), IADIAG, JADIAG, DESCA, B, IBPOS,
     $                JBPOS, DESCB )
* ------------------------------------------------------
* use the newly computed X(:) to update the rest of the
* entries in B
* ------------------------------------------------------
         IF( ISAXPYUPDATE ) THEN
* column update in B
            IF( ISFORWARD ) THEN
* the rest of the vector
               IF( ISLEFT ) THEN
                  IBPOS = ( IB-1 ) + ( JEND+1 )
                  JBPOS = ( JB-1 ) + 1
               ELSE
                  IBPOS = ( IB-1 ) + 1
                  JBPOS = ( JB-1 ) + ( JEND+1 )
               ENDIF
               IBSIZE = JENDLR - ( JEND+1 ) + 1
            ELSE
* the rest of the vector
               IBPOS = ( IB-1 ) + 1
               JBPOS = ( JB-1 ) + 1
               IBSIZE = JSTART - 1
            ENDIF
            IF( ISLOWER ) THEN
               IA = JSIZE + 1
               JA = 1
            ELSE
               IF( ISUPPER ) THEN
                  IA = 1
                  JA = 1
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
                  CALL PZGEMM( TRANSA, 'N', MM, NN, KK, LALPHA,
     $                         WORK( IP_AMAT ), IA, JA, DESCA, B,
     $                         IIBPOS, JJBPOS, DESCB, LBETA, B, IBPOS,
     $                         JBPOS, DESCB )
               ELSE
                  CALL PZGEMM( 'N', TRANSA, MM, NN, KK, LALPHA, B,
     $                         IIBPOS, JJBPOS, DESCB, WORK( IP_AMAT ),
     $                         IA, JA, DESCA, LBETA, B, IBPOS, JBPOS,
     $                         DESCB )
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
  100 CONTINUE
* end while
* -----------------------------------
*  need to rescale the solution by alpha
* -----------------------------------
      IF( ALPHA.NE.ONE ) THEN
         DO 110 I = 1, N
            IBPOS = ( IB-1 ) + 1
            JBPOS = ( JB-1 ) + I
            INC = 1
            CALL PZSCAL( M, ALPHA, B, IBPOS, JBPOS, DESCB, INC )
  110    CONTINUE
  120    CONTINUE
      ENDIF
      WORK( 1 ) = INEED
      RETURN
      END
