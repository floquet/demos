      SUBROUTINE PFCGETF2( M, N, CWORK, IC, JC, DESCC, IPIV, INFO )
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
*  PFGETF2 computes an LU factorization of a general M-by-N
*  out-of-core matrix sub( C ) = C(IC:IC+M-1,JC:JC+N-1) using
*  partial pivoting with row interchanges.
*
*  The factorization has the form sub( A ) = P * L * U, where P is a
*  permutation matrix, L is lower triangular with unit diagonal
*  elements (lower trapezoidal if m > n), and U is upper triangular
*  (upper trapezoidal if m < n).
*
*  This is the left-looking Parallel version of the algorithm.
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
      LOGICAL            FULLPIVOTB
      PARAMETER          ( FULLPIVOTB = .false. )
      LOGICAL            USE_UNROLL
      PARAMETER          ( USE_UNROLL = .true. )
*     ..
*     .. Scalar Arguments ..
      INTEGER            IC, INFO, JC, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCC( FLEN_ ), IPIV( * )
      COMPLEX            CWORK( * )
*     ..
*     .. Local Scalars ..
      LOGICAL            HASWORK, ISALIGNED, ISVALID, LWORKQUERY
      INTEGER            AEND, ANEED, BNEED, BPIVNEED, CSRC, CSRCIO,
     $                   FIRSTSTART, I, IA, IAEND, IASIZE, IASTART, IB,
     $                   ICOLB, ICOLC, ICONTXT, ICTXT, IDUMMY, IENDA,
     $                   IENDB, IFREE, IIC, IINFO, IIP, INC, IODEV,
     $                   IPROC, IP_AMAT, IP_BMAT, IP_BPIV, IP_BPIVEND,
     $                   IROWB, IROWC, ISIZEA, ISIZEB, ISTARTA, ISTARTB,
     $                   JA, JB, JENDA, JENDB, JJC, JJEND, JJSTART,
     $                   JSIZEA, JSIZEB, JSTARTA, JSTARTB, KK,
     $                   LASTINDEX, LDA, LDB, LDPIV, LOCPA, LOCPB,
     $                   LOCQA, LOCQB, LWORK, MB, MINMN, MM, MMB,
     $                   MYPCOL, MYPROW, NB, NCOLA, NCOLB, NFREE, NN,
     $                   NNB, NPANELS, NPCOL, NPROW, NROWA, NROWB, P0,
     $                   Q0, RSRC, RSRCIO, SAVEDT
      COMPLEX            ALPHA, BETA, ONE
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA( DLEN_ ), DESCB( DLEN_ ),
     $                   DESCIPIV( DLEN_ )
*     ..
*     .. External Functions ..
      INTEGER            ICEIL, INDXG2P, NUMROC
      EXTERNAL           ICEIL, INDXG2P, NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           ASSERT, BLACS_GRIDINFO, CHK1MAT, CLAREAD,
     $                   CLAWRITE, DESCINIT, INFOG1L, LAIO_INFO, PCGEMM,
     $                   PCGETRF, PCHK1MAT, PCLAPIV, PCTRSM, PFMAXSIZE,
     $                   PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, INT, MAX, MIN, MOD, REAL
*     ..
*     .. Executable Statements ..
      ONE = CMPLX( REAL( 1 ) )
      MINMN = MIN( M, N )
      LWORK = DESCC( SIZE_ )
      LWORKQUERY = ( LWORK.EQ.-1 )
      IODEV = DESCC( IODEV_ )
      ICTXT = DESCC( CTXT_ )
      CALL BLACS_GRIDINFO( DESCC( CTXT_ ), NPROW, NPCOL, MYPROW,
     $                     MYPCOL )
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 6*100+CTXT_ )
         RETURN
      ENDIF
      SAVEDT = DESCC( DT_ )
      DESCC( DT_ ) = BLOCK_CYCLIC_2D
      CALL CHK1MAT( M, 1, N, 2, IC, JC, DESCC, 6, INFO )
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCC( CTXT_ ), 'PFCGETF2', 6 )
         RETURN
      ENDIF
      IDUMMY = 1
      CALL PCHK1MAT( M, 1, N, 2, IC, JC, DESCC, 6, IDUMMY, IDUMMY,
     $               IDUMMY, INFO )
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCC( CTXT_ ), 'PFCGETF2', 6 )
         RETURN
      ENDIF
      DESCC( DT_ ) = SAVEDT
      IF( DESCC( DT_ ).NE.DISK_BLOCK_CYCLIC_2D ) THEN
         INFO = -( 6*100+DT_ )
         CALL PXERBLA( DESCC( CTXT_ ), 'PFCGETF2', 6 )
         RETURN
      ENDIF
      INFO = 0
      CALL LAIO_INFO( DESCC( IODEV_ ), MM, NN, MMB, NNB, MB, NB, CSRC,
     $                RSRC, ICONTXT )
      IF( MB.NE.DESCC( MB_ ) ) THEN
         INFO = -( 6*100+MB_ )
      ENDIF
      IF( NB.NE.DESCC( NB_ ) ) THEN
         INFO = -( 6*100+NB_ )
      ENDIF
      IF( ICONTXT.NE.DESCC( CTXT_ ) ) THEN
         INFO = -( 6*100+CTXT_ )
      ENDIF
      IF( DESCC( RSRC_ ).NE.RSRC ) THEN
         INFO = -( 6*100+RSRC_ )
      ENDIF
      IF( DESCC( CSRC_ ).NE.CSRC ) THEN
         INFO = -( 6*100+CSRC_ )
      ENDIF
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCC( CTXT_ ), 'PFCGETF2', 6 )
         RETURN
      ENDIF
      CALL LAIO_INFO( IODEV, MM, NN, MMB, NNB, MB, NB, CSRCIO, RSRCIO,
     $                ICONTXT )
      CALL BLACS_GRIDINFO( ICONTXT, NPROW, NPCOL, MYPROW, MYPCOL )
*
* Going to reuse storage in piv.
* Double check there is sufficient space.
*
      IF( IC.NE.1 ) THEN
         CALL INFOG1L( IC, MB, NPROW, MYPROW, RSRCIO, IIP, IPROC )
         IP_BPIV = IIP
         CALL INFOG1L( ( IC-1 )+M, MB, NPROW, MYPROW, RSRCIO, IIP,
     $                 IPROC )
         IF( IPROC.NE.MYPROW ) THEN
            IIP = IIP - 1
         ENDIF
         IP_BPIVEND = IIP
         IROWC = INDXG2P( IC, MB, MYPROW, RSRCIO, NPROW )
         BPIVNEED = NUMROC( M, MB, MYPROW, IROWC, NPROW )
         CALL ASSERT( ( IP_BPIVEND-IP_BPIV+1 ).GE.BPIVNEED,
     $                '**  PFxGETF2: insufficient storage in ipiv(*) ',
     $                ( IP_BPIVEND-IP_BPIV+1 ) )
      ENDIF
*
*  First panel is special case.
*  Try to use all available memory.
*
      ANEED = 0
      INFO = 0
*
*   Align B to perform local I/O operations.
*
      IROWC = INDXG2P( IC, MB, MYPROW, RSRCIO, NPROW )
      ICOLC = INDXG2P( JC, NB, MYPCOL, CSRCIO, NPCOL )
      P0 = IROWC
      Q0 = ICOLC
      NROWB = M
      NCOLB = -1
      CALL PFMAXSIZE( 'C', LWORK, NROWB, NCOLB, MB, NB, P0, Q0, ICONTXT,
     $                INFO )
      ISVALID = ( INFO.EQ.0 ) .OR. ( LWORKQUERY .AND. ( INFO.EQ.-2 ) )
      CALL ASSERT( ISVALID, '**  PFxGETF2**: pfmaxsize returns ', INFO )
      JSTARTB = 1
      JENDB = MIN( N, JSTARTB+NCOLB-1 )
      IF( JENDB.NE.N ) THEN
         CALL ASSERT( NCOLB.GE.NNB, '** PFGETF2: ncolB < nnb, ncolB =  '
     $                , NCOLB )
*
*       Try to align to block boundary.
*
         JJC = ( JC-1 ) + JENDB
         JJC = MAX( ( JC-1 )+1, INT( JJC / NNB )*NNB )
         JENDB = JJC - ( JC-1 )
      ENDIF
      JSIZEB = JENDB - JSTARTB + 1
      NCOLB = JSIZEB
      LOCPB = NUMROC( NROWB, MB, MYPROW, P0, NPROW )
      LOCQB = NUMROC( NCOLB, NB, MYPCOL, Q0, NPCOL )
      LDB = MAX( 1, LOCPB )
      IFREE = 1
      BNEED = LOCPB*LOCQB
      IP_BMAT = IFREE
      IFREE = IFREE + BNEED
      CALL ASSERT( LWORK.GE.BNEED,
     $             '**  PFxGETF2: Bneed > lwork, Bneed = ', BNEED )
      CALL DESCINIT( DESCB, NROWB, NCOLB, MB, NB, P0, Q0, ICONTXT, LDB,
     $               INFO )
      CALL ASSERT( INFO.EQ.0, '**  PFxGETF2: descinit returns info = ',
     $             INFO )
*
*   In place use of storage for pivot vector
*
      IF( IC.EQ.1 ) THEN
         IP_BPIV = 1
      ELSE
         CALL INFOG1L( IC, MB, NPROW, MYPROW, RSRCIO, IIP, IPROC )
         IP_BPIV = IIP
      ENDIF
*
*   Read into panel.
*
      INFO = 0
      CALL CLAREAD( IODEV, M, NCOLB, IC, JC, CWORK( IP_BMAT ), 1, 1,
     $              DESCB, INFO )
      CALL ASSERT( INFO.EQ.0, '**  PFxGETF2: LAREAD returns info = ',
     $             INFO )
*
*   Perform LU factorization.
*
      IINFO = 0
      CALL PCGETRF( M, NCOLB, CWORK( IP_BMAT ), 1, 1, DESCB,
     $              IPIV( IP_BPIV ), IINFO )
      IF( IINFO.GT.0 ) THEN
*
*        Matrix may be singular
*
         INFO = IINFO
         LASTINDEX = ( JSTARTB-1 ) + IINFO
         GOTO 90
      ENDIF
      CALL ASSERT( INFO.EQ.0, '** PGETF2: PGETRF returns info ', INFO )
*
*   Write out factors.
*
      INFO = 0
      CALL CLAWRITE( IODEV, M, NCOLB, IC, JC, CWORK( IP_BMAT ), 1, 1,
     $               DESCB, INFO )
      CALL ASSERT( INFO.EQ.0, '**  PFxGETF2: LAWRITE returns info = ',
     $             INFO )
      IF( JENDB.EQ.N ) THEN
*
*         whole matrix fit in core.
*         all done.
*
         INFO = 0
         LASTINDEX = M
         GOTO 90
      ENDIF
*
*  Factor the rest of matrix.
*
*
*        Need storage for one (m x nnb) panel A,
*        and storage for another (m x nnb) panel B
*
*
*        Note the use of p0,q0 to overestimage
*        storage requirement
*
      P0 = MYPROW
      Q0 = MYPCOL
      NROWA = M
      NCOLA = MIN( N, NNB )
      LOCPA = NUMROC( NROWA, MB, MYPROW, P0, NPROW )
      LOCQA = NUMROC( NCOLA, NB, MYPROW, Q0, NPCOL )
      ANEED = LOCPA*LOCQA
*
*        Minimum requirement for B
*        is another m x nnb pannel
*
      P0 = MYPROW
      Q0 = MYPCOL
      NROWB = M
      NCOLB = MIN( ( N-( JENDB+1 )+1 ), NNB )
      LOCPB = NUMROC( NROWB, MB, MYPROW, P0, NPROW )
      LOCQB = NUMROC( NCOLB, NB, MYPCOL, Q0, NPCOL )
      BNEED = LOCPB*LOCQB
      IF( ANEED+BNEED.GT.LWORK ) THEN
         INFO = -( 6*100+SIZE_ )
         IF( .NOT.LWORKQUERY ) THEN
            CALL PXERBLA( ICONTXT, 'PFCGETF2', 6 )
         ENDIF
         LASTINDEX = ( JSTARTB-1 )
         GOTO 90
      ENDIF
*
*        The rest of storage for B.
*
      NROWB = M
      NCOLB = -1
      NFREE = LWORK - ANEED
      CALL PFMAXSIZE( 'C', NFREE, NROWB, NCOLB, MB, NB, P0, Q0, ICONTXT,
     $                INFO )
      CALL ASSERT( INFO.EQ.0, '**  PFxGETF2: pfmaxsize returns info = ',
     $             INFO )
      NPANELS = 1
      IF( ( N-( JENDB+1 )+1 ).GE.NNB ) THEN
*
*               Need at least one (m x nnb) panel.
*
         NPANELS = INT( NCOLB / NNB )
         CALL ASSERT( NPANELS.GE.1,
     $                '**  PFxGETF2: npanels <= 0, npanels ', NPANELS )
         NCOLB = NPANELS*NNB
      ENDIF
      IFREE = 1
      IP_AMAT = IFREE
      IFREE = IFREE + ANEED
      IP_BMAT = IFREE
      JSTARTB = JENDB + 1
      FIRSTSTART = JSTARTB
      CALL ASSERT( MOD( FIRSTSTART-1, NNB ).EQ.0,
     $             '** PFGETF2: invalid firststart ', FIRSTSTART )
*
*               Align B for local I/O
*
      IIC = ( IC-1 ) + IC
      JJC = ( JC-1 ) + JSTARTB
      IROWB = INDXG2P( IIC, MB, MYPROW, RSRCIO, NPROW )
      ICOLB = INDXG2P( JJC, NB, MYPCOL, CSRCIO, NPCOL )
      P0 = IROWB
      Q0 = ICOLB
      LOCPB = NUMROC( M, MB, MYPROW, P0, NPROW )
      LOCQB = NUMROC( NCOLB, NB, MYPCOL, Q0, NPCOL )
      LDB = MAX( 1, LOCPB )
      BNEED = LOCPB*LOCQB
      CALL ASSERT( BNEED.LE.( LWORK-ANEED ),
     $             '**  PFxGETF2: invalid Bneed = ', BNEED )
      CALL DESCINIT( DESCB, M, NCOLB, MB, NB, P0, Q0, ICONTXT, LDB,
     $               INFO )
      CALL ASSERT( INFO.EQ.0, '**  PFxGETF2: descinit returns info = ',
     $             INFO )
*
*               Align panel A to B
*               I/O should also be local.
*
      AEND = MIN( MINMN, ( JSTARTB-1 ) )
      NROWA = M
      NCOLA = MIN( NNB, AEND )
      LOCPA = NUMROC( NROWA, MB, MYPROW, P0, NPROW )
      LOCQA = NUMROC( NCOLA, NB, MYPCOL, Q0, NPCOL )
      LDA = MAX( 1, LOCPA )
      CALL ASSERT( LOCPA*LOCQA.LE.ANEED,
     $             '**  PFxGETF2: internal error, LocpA*LocqA > Aneed ',
     $             LOCPA*LOCQA )
      CALL DESCINIT( DESCA, NROWA, NCOLA, MB, NB, P0, Q0, ICONTXT, LDA,
     $               INFO )
      CALL ASSERT( INFO.EQ.0, '**  PFxGETF2: descinit returns info = ',
     $             INFO )
*
*           Note the use of mypcol and not descB(CSRC_)
*           in setting up descipiv(*)
*
      LOCPB = NUMROC( DESCB( M_ ), DESCB( MB_ ), MYPROW, DESCB( RSRC_ ),
     $        NPROW )
      LDPIV = MAX( 1, LOCPB+MB )
      CALL DESCINIT( DESCIPIV, DESCB( M_ )+NPROW*DESCB( MB_ ), 1,
     $               DESCB( MB_ ), 1, DESCB( RSRC_ ), MYPCOL,
     $               DESCB( CTXT_ ), LDPIV, INFO )
*
*        Note jstartB <= n, should work even for
*        rectangular matrix.
*
   10 CONTINUE
      IF( JSTARTB.LE.N ) THEN
         JENDB = MIN( N, JSTARTB+DESCB( N_ )-1 )
         IF( JENDB.NE.N ) THEN
            JJC = ( JC-1 ) + JENDB
            JJC = MAX( ( JC-1 )+JSTARTB, INT( JJC / NNB )*NNB )
            JENDB = JJC - ( JC-1 )
         ENDIF
         JSIZEB = JENDB - JSTARTB + 1
*
*               Read in panel B.
*
         IIC = ( IC-1 ) + IC
         JJC = ( JC-1 ) + JSTARTB
         CALL CLAREAD( IODEV, M, JSIZEB, IIC, JJC, CWORK( IP_BMAT ), 1,
     $                 1, DESCB, INFO )
         CALL ASSERT( INFO.EQ.0, '**  PFxGETF2: LAREAD returns info = ',
     $                INFO )
*
*        Perform permutation.
*
         IF( FULLPIVOTB ) THEN
*
*             Perform all permutations.
*
            MM = MIN( M, ( JSTARTB-1 ) )
            NN = JSIZEB
         ELSE
*
*           Perform permutations for special case
*           in first 2 panels.
*
            MM = MIN( M, FIRSTSTART-1 )
            NN = JSIZEB
         ENDIF
         HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
         IF( HASWORK ) THEN
            CALL PCLAPIV( 'Forward', 'Rowise', 'Columwise', MM, NN,
     $                    CWORK( IP_BMAT ), 1, 1, DESCB,
     $                    IPIV( IP_BPIV ), 1, 1, DESCIPIV, IDUMMY )
         ENDIF
         AEND = MIN( MINMN, ( JSTARTB-1 ) )
         JSTARTA = 1
   20    CONTINUE
         IF( JSTARTA.LE.AEND ) THEN
            ISTARTA = JSTARTA
            IENDA = M
            ISIZEA = IENDA - ISTARTA + 1
            IF( ISIZEA.LE.0 ) THEN
*
*                       Rectangular matrix.
*                       No more updates.
*
               GOTO 50
            ENDIF
            NCOLA = DESCA( N_ )
            JENDA = MIN( AEND, JSTARTA+NCOLA-1 )
            IF( JENDA.NE.( JSTARTB-1 ) ) THEN
               JJC = ( JC-1 ) + JENDA
               JJC = MAX( ( JC-1 )+JSTARTA, INT( JJC / NNB )*NNB )
               JENDA = JJC - ( JC-1 )
            ENDIF
            JENDA = MAX( JSTARTA, JENDA )
            JSIZEA = JENDA - JSTARTA + 1
*
*                 Read in panel A
*
            IIC = ( IC-1 ) + ISTARTA
            JJC = ( JC-1 ) + JSTARTA
            IA = ISTARTA
            JA = 1
            CALL CLAREAD( IODEV, ISIZEA, JSIZEA, IIC, JJC,
     $                    CWORK( IP_AMAT ), IA, JA, DESCA, INFO )
            CALL ASSERT( INFO.EQ.0, '**  PFxGETF2: LAREAD returns info '
     $                   , INFO )
*
*        Permute lower factors
*
            IF( JSTARTA.LE.( FIRSTSTART-1 ) ) THEN
*
*                        Special case.
*
               JJSTART = FIRSTSTART
            ELSE
               JJSTART = ( FIRSTSTART-1 ) +
     $                   ICEIL( ( JSTARTA-( FIRSTSTART-1 ) ), NCOLB )*
     $                   NCOLB + 1
            ENDIF
            IF( FULLPIVOTB ) THEN
*
*                        Permute panel A to match panel B.
*
               MM = MIN( M, ( JSTARTB-1 ) ) - JJSTART + 1
               NN = JSIZEA
               HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
               IF( HASWORK ) THEN
                  CALL PCLAPIV( 'Forward', 'Rowwise', 'Columnwise', MM,
     $                          NN, CWORK( IP_AMAT ), JJSTART, 1, DESCB,
     $                          IPIV( IP_BPIV ), JJSTART, 1, DESCIPIV,
     $                          IDUMMY )
               ENDIF
            ELSE
*
*                         Permute panel B to match panel A.
*
               ISALIGNED = ( 1.EQ.MOD( ( JSTARTA-( FIRSTSTART-1 ) ),
     $                     NPANELS*NNB ) )
               IF( ( JSTARTA.GE.FIRSTSTART ) .AND. ( ISALIGNED ) ) THEN
*
*                        Past the first 2 panels.
*
                  JJEND = MAX( JJSTART-1, JENDA )
                  MM = JJEND - JSTARTA + 1
                  NN = JSIZEB
                  HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
                  IF( HASWORK ) THEN
                     CALL PCLAPIV( 'Forward', 'Rowwise', 'Columnwise',
     $                             MM, NN, CWORK( IP_BMAT ), JSTARTA, 1,
     $                             DESCB, IPIV( IP_BPIV ), JSTARTA, 1,
     $                             DESCIPIV, IDUMMY )
                  ENDIF
               ENDIF
            ENDIF
*
*  Update panel B with panel A.
*
            IF( USE_UNROLL ) THEN
*
*        Unroll computation to treat one mb by mb block at a time.
*
               MB = DESCA( MB_ )
               INC = MB
               DO 30 IASTART = JSTARTA, JENDA, INC
                  IAEND = MIN( JENDA, IASTART+INC-1 )
                  IF( IAEND.NE.JENDA ) THEN
*
*               Align to block boundary.
*
                     IAEND = MAX( IASTART, INT( IAEND / MB )*MB )
                  ENDIF
                  IASIZE = IAEND - IASTART + 1
                  ALPHA = ONE
                  MM = IASIZE
                  NN = JSIZEB
                  IA = IASTART
                  JA = IASTART - JSTARTA + 1
                  IB = IASTART
                  JB = 1
                  HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
                  IF( HASWORK ) THEN
                     CALL PCTRSM( 'Leftside', 'LowerTriangular',
     $                            'NoTranspose', 'UnitDiagonal', MM, NN,
     $                            ALPHA, CWORK( IP_AMAT ), IA, JA,
     $                            DESCA, CWORK( IP_BMAT ), IB, JB,
     $                            DESCB )
                  ENDIF
                  MM = M - ( IAEND+1 ) + 1
                  NN = JSIZEB
                  KK = IASIZE
                  IA = IAEND + 1
                  JA = IASTART - JSTARTA + 1
                  IB = IASTART
                  JB = 1
                  HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 ) .AND.
     $                      ( KK.GE.1 )
                  IF( HASWORK ) THEN
                     ALPHA = -ONE
                     BETA = ONE
                     CALL PCGEMM( 'Notrans', 'NoTrans', MM, NN, KK,
     $                            ALPHA, CWORK( IP_AMAT ), IA, JA,
     $                            DESCA, CWORK( IP_BMAT ), IB, JB,
     $                            DESCB, BETA, CWORK( IP_BMAT ), IA, JB,
     $                            DESCB )
                  ENDIF
   30          CONTINUE
   40          CONTINUE
            ELSE
*
*     Update upper part.
*
*     B(jstartA:jendA, 1:jsizeB) <-
*           A(jstartA:jendA, 1:jsizeA) \ B( jstartA:jendA, 1:jsizeB )
*
               IA = JSTARTA
               JA = 1
               IB = JSTARTA
               JB = 1
               MM = JSIZEA
               NN = JSIZEB
               HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
               IF( HASWORK ) THEN
                  ALPHA = ONE
                  CALL PCTRSM( 'Leftside', 'LowerTriangular',
     $                         'NoTranspose', 'UnitDiagonal', MM, NN,
     $                         ALPHA, CWORK( IP_AMAT ), IA, JA, DESCA,
     $                         CWORK( IP_BMAT ), IB, JB, DESCB )
               ENDIF
*
*       Update the rest of the column.
*
*       B( jendA+1:m, 1:jsizeB) <-
*               B( jendA+1:m, 1:jsizeB) -
*                       A(jendA+1:m, 1:jsizeA) *
*                               B(jstartA:jendA, 1:jsizeB)
*
               MM = M - ( JENDA+1 ) + 1
               NN = JSIZEB
               KK = JSIZEA
               HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 ) .AND. ( KK.GE.1 )
               IF( HASWORK ) THEN
                  ALPHA = -ONE
                  BETA = ONE
                  CALL PCGEMM( 'Notrans', 'NoTrans', MM, NN, KK, ALPHA,
     $                         CWORK( IP_AMAT ), JENDA+1, 1, DESCA,
     $                         CWORK( IP_BMAT ), JSTARTA, 1, DESCB,
     $                         BETA, CWORK( IP_BMAT ), JENDA+1, 1,
     $                         DESCB )
               ENDIF
            ENDIF
            JSTARTA = JENDA + 1
            GOTO 20
         ENDIF
   50    CONTINUE
* end while
*
*        All updates are done.
*        Factor panel B.
*
         MM = M - JSTARTB + 1
         NN = JSIZEB
         HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
         IF( HASWORK ) THEN
            IINFO = 0
            IB = JSTARTB
            JB = 1
            CALL PCGETRF( MM, NN, CWORK( IP_BMAT ), IB, JB, DESCB,
     $                    IPIV( IP_BPIV ), IINFO )
            IF( IINFO.GT.0 ) THEN
*
*                Check for near singularity.
*
               INFO = ( JSTARTB-1 ) + IINFO
               LASTINDEX = ( JSTARTB-1 ) + IINFO
               GOTO 90
            ENDIF
            CALL ASSERT( INFO.EQ.0, '** PGETF2: PGETRF returns info ',
     $                   INFO )
         ENDIF
*
*        Write out panel B.
*
         IIC = ( IC-1 ) + IC
         JJC = ( JC-1 ) + JSTARTB
         CALL CLAWRITE( IODEV, M, JSIZEB, IIC, JJC, CWORK( IP_BMAT ), 1,
     $                  1, DESCB, INFO )
         CALL ASSERT( INFO.EQ.0, '**  PFxGETF2: LAWRITE returns info = '
     $                , INFO )
         JSTARTB = JENDB + 1
         GOTO 10
      ENDIF
   60 CONTINUE
*
*   Final pass to reorder factors.
*   reuse storage for Bpiv(*)
*
      JSTARTB = FIRSTSTART
      MINMN = MIN( M, N )
*
*   Note jstartB <= min(m,n), not n.
*   Since only the Lower factors need to be reordered.
*
   70 CONTINUE
      IF( JSTARTB.LE.MINMN ) THEN
         JENDB = MIN( MINMN, JSTARTB+DESCB( N_ )-1 )
         IF( JENDB.NE.MINMN ) THEN
            JJC = ( JC-1 ) + JENDB
            JJC = MAX( ( JC-1 )+JSTARTB, INT( JJC / NNB )*NNB )
            JENDB = JJC - ( JC-1 )
         ENDIF
         JSIZEB = JENDB - JSTARTB + 1
         ISTARTB = JENDB + 1
         IENDB = M
         ISIZEB = IENDB - ISTARTB + 1
         MM = MIN( M, N ) - ISTARTB + 1
         NN = JSIZEB
         HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
         IF( HASWORK ) THEN
            IIC = ( IC-1 ) + ISTARTB
            JJC = ( JC-1 ) + JSTARTB
            IB = ISTARTB
            JB = 1
            MM = M - ISTARTB + 1
            NN = JSIZEB
            CALL CLAREAD( IODEV, MM, NN, IIC, JJC, CWORK( IP_BMAT ),
     $                    ISTARTB, JB, DESCB, INFO )
            CALL ASSERT( INFO.EQ.0,
     $                   '**  PFxGETF2: LAREAD returns info = ', INFO )
            MM = MIN( M, N ) - ISTARTB + 1
            NN = JSIZEB
            CALL PCLAPIV( 'Forward', 'Rowise', 'Columwise', MM, NN,
     $                    CWORK( IP_BMAT ), IB, 1, DESCB,
     $                    IPIV( IP_BPIV ), IB, 1, DESCIPIV, IDUMMY )
            MM = M - ISTARTB + 1
            NN = JSIZEB
            CALL CLAWRITE( IODEV, MM, NN, IIC, JJC, CWORK( IP_BMAT ),
     $                     IB, JB, DESCB, INFO )
            CALL ASSERT( INFO.EQ.0,
     $                   '**  PFxGETF2: LAWRITE returns info = ', INFO )
         ENDIF
         JSTARTB = JENDB + 1
         GOTO 70
      ENDIF
   80 CONTINUE
*
*  Special case for first 2 panels.
*
      IP_BMAT = 1
      IROWC = INDXG2P( IC, MB, MYPROW, RSRCIO, NPROW )
      ICOLC = INDXG2P( JC, NB, MYPCOL, CSRCIO, NPCOL )
      P0 = IROWC
      Q0 = ICOLC
      NCOLB = FIRSTSTART - 1
      LOCPB = NUMROC( M, MB, MYPROW, P0, NPROW )
      LOCQB = NUMROC( NCOLB, NB, MYPCOL, Q0, NPCOL )
      LDB = MAX( 1, LOCPB )
      CALL DESCINIT( DESCB, M, NCOLB, MB, NB, P0, Q0, ICONTXT, LDB,
     $               INFO )
*
*            Note the use of mypcol and not descB(CSRC_)
*            in setting up descipiv(*)
*
      LOCPB = NUMROC( DESCB( M_ ), DESCB( MB_ ), MYPROW, DESCB( RSRC_ ),
     $        NPROW )
      LDPIV = MAX( 1, LOCPB+MB )
      CALL DESCINIT( DESCIPIV, DESCB( M_ )+NPROW*DESCB( MB_ ), 1,
     $               DESCB( MB_ ), 1, DESCB( RSRC_ ), MYPCOL,
     $               DESCB( CTXT_ ), LDPIV, INFO )
      JSTARTB = 1
      JENDB = FIRSTSTART - 1
      JSIZEB = JENDB - JSTARTB + 1
      ISTARTB = JENDB + 1
      IENDB = M
      ISIZEB = IENDB - ISTARTB + 1
      MM = MIN( M, N ) - ISTARTB + 1
      NN = JSIZEB
      HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
      IF( HASWORK ) THEN
         IIC = ( IC-1 ) + ISTARTB
         JJC = ( JC-1 ) + JSTARTB
         IB = ISTARTB
         JB = 1
         MM = M - ISTARTB + 1
         CALL CLAREAD( IODEV, MM, NN, IIC, JJC, CWORK( IP_BMAT ), IB,
     $                 JB, DESCB, INFO )
         CALL ASSERT( INFO.EQ.0, '**  PFxGETF2: LAREAD returns info = ',
     $                INFO )
         MM = MIN( M, N ) - ISTARTB + 1
         CALL PCLAPIV( 'Forward', 'Rowise', 'Columwise', MM, NN,
     $                 CWORK( IP_BMAT ), IB, 1, DESCB, IPIV( IP_BPIV ),
     $                 IB, 1, DESCIPIV, IDUMMY )
         MM = M - ISTARTB + 1
         CALL CLAWRITE( IODEV, MM, NN, IIC, JJC, CWORK( IP_BMAT ), IB,
     $                  JB, DESCB, INFO )
         CALL ASSERT( INFO.EQ.0, '**  PFxGETF2: LAWRITE returns info = '
     $                , INFO )
      ENDIF
*
*   All done.
*
      INFO = 0
      LASTINDEX = M
      GOTO 90
   90 CONTINUE
      CALL INFOG1L( ( IC-1 )+LASTINDEX, MB, NPROW, MYPROW, RSRCIO, IIP,
     $              IPROC )
      IF( IPROC.NE.MYPROW ) THEN
         IIP = IIP - 1
      ENDIF
      IP_BPIVEND = IIP
      IF( IC.NE.1 ) THEN
*
*        retore offset in pivot vector.
*
         DO 100 I = IP_BPIV, IP_BPIVEND
            IPIV( I ) = IPIV( I ) + ( IC-1 )
  100    CONTINUE
  110    CONTINUE
      ENDIF
      CWORK( 1 ) = ANEED + BNEED
      RETURN
      END
