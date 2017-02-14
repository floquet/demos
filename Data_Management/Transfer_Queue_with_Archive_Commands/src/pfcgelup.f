      SUBROUTINE PFCGELUP( M, N, CWORK, IC, JC, DESCC, IPIV, JSIZEB,
     $                     BMAT, IB, JB, DESCB, FIRSTSTART, NCOLB )
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
*  Perform update by out-of-core lower triangular factors.
*  The lower factors are stored on disk in a partially pivoted
*  manner. Pivoting is applied to the target array to match the
*  pivot order of the lower triangular panels.
*
*  submatrix  B( ib:ib+m-1, jb:jb+jsizeB-1) is modified.
*
*  The routine is used in PFGELUP and PFGETRS.
*
*  Note that firststart, ncolB are parameters associated with
*  the PFGETF2 factorization.
*
*
*     .. Parameters ..
      INTEGER            DLEN_
      PARAMETER          ( DLEN_ = 9 )
      INTEGER            CTXT_, M_, N_, MB_
      PARAMETER          ( CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5 )
      INTEGER            RSRC_
      PARAMETER          ( RSRC_ = 7 )
      INTEGER            IODEV_, SIZE_
      PARAMETER          ( IODEV_ = 10, SIZE_ = 11 )
      INTEGER            FLEN_
      PARAMETER          ( FLEN_ = 12 )
*     ..
*     .. Scalar Arguments ..
      INTEGER            FIRSTSTART, IB, IC, JB, JC, JSIZEB, M, N, NCOLB
*     ..
*     .. Array Arguments ..
      INTEGER            DESCB( DLEN_ ), DESCC( FLEN_ ), IPIV( * )
      COMPLEX            BMAT( * ), CWORK( * )
*     ..
*     .. Local Scalars ..
      LOGICAL            HASWORK, ISALIGNED
      INTEGER            AEND, CSRCIO, IA, IAEND, IASIZE, IASTART,
     $                   ICOLA, ICOLC, ICONTXT, IDUMMY, IENDA, IIB, IIC,
     $                   IIP, INC, INFO, IODEV, IPROC, IP_AMAT, IP_BPIV,
     $                   IROWA, IROWC, ISIZEA, ISTARTA, JA, JENDA, JJB,
     $                   JJC, JJEND, JJSTART, JSIZEA, JSTARTA, KK, LDA,
     $                   LDIPIV, LOCP, LWORK, MB, MINMN, MIPIV, MM, MMB,
     $                   MYPCOL, MYPROW, NB, NCOLA, NN, NNB, NPANELS,
     $                   NPCOL, NPROW, NROWA, P0, Q0, RSRCIO
      COMPLEX            ALPHA, BETA, ONE, ZERO
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA( DLEN_ ), DESCIPIV( DLEN_ )
*     ..
*     .. External Functions ..
      INTEGER            ICEIL, INDXG2P, NUMROC
      EXTERNAL           ICEIL, INDXG2P, NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           ASSERT, BLACS_GRIDINFO, CLAREAD, DESCINIT,
     $                   INFOG1L, LAIO_INFO, PCGEMM, PCLAPIV, PCTRSM,
     $                   PFMAXSIZE
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, INT, MAX, MIN, MOD, REAL
*     ..
*     .. Executable Statements ..
      ONE = CMPLX( REAL( 1 ) )
      ZERO = CMPLX( REAL( 0 ) )
      IODEV = DESCC( IODEV_ )
      LWORK = DESCC( SIZE_ )
      MINMN = MIN( M, N )
      IP_AMAT = 1
*
*     Note the use of mypcol instead of descC(CSRC_).
*     and extra descC(MB_) storage at end.
*
      CALL LAIO_INFO( IODEV, MM, NN, MMB, NNB, MB, NB, CSRCIO, RSRCIO,
     $                ICONTXT )
      CALL BLACS_GRIDINFO( DESCC( CTXT_ ), NPROW, NPCOL, MYPROW,
     $                     MYPCOL )
      NPANELS = MAX( 1, INT( NCOLB / NNB ) )
      CALL ASSERT( MOD( NCOLB, NNB ).EQ.0,
     $             '**  PFxGELUP **: mod(ncolB,nnb) != 0, ncolB= ',
     $             NCOLB )
*
*       Initialize descriptor for panel A.
*
      IROWC = INDXG2P( IC, MB, MYPROW, RSRCIO, NPROW )
      ICOLC = INDXG2P( JC, NB, MYPCOL, CSRCIO, NPCOL )
      IROWA = IROWC
      ICOLA = ICOLC
      P0 = IROWC
      Q0 = ICOLC
      NROWA = M
      NCOLA = -1
      INFO = 0
      CALL PFMAXSIZE( 'C', LWORK, NROWA, NCOLA, MB, NB, P0, Q0, ICONTXT,
     $                INFO )
      CALL ASSERT( INFO.EQ.0, '**  PFxGELUP **: pfmaxsize returns ',
     $             INFO )
      CALL ASSERT( NCOLA.GE.NNB,
     $             '**  PFxGELUP **: ncolA < nnb, ncolA = ', NCOLA )
      NCOLA = MAX( 1, MIN( MINMN, NNB ) )
      LDA = MAX( 1, NUMROC( NROWA, MB, MYPROW, P0, NPROW ) )
      CALL DESCINIT( DESCA, NROWA, NCOLA, MB, NB, P0, Q0, ICONTXT, LDA,
     $               INFO )
*
*       Initialize descriptor for pivot vector.
*
      MIPIV = DESCC( M_ ) + DESCC( MB_ )*NPROW
      LOCP = NUMROC( DESCC( M_ ), DESCC( MB_ ), MYPROW, DESCC( RSRC_ ),
     $       NPROW )
      LDIPIV = LOCP + DESCC( MB_ )
      CALL DESCINIT( DESCIPIV, MIPIV, 1, DESCC( MB_ ), 1,
     $               DESCC( RSRC_ ), MYPCOL, DESCC( CTXT_ ), LDIPIV,
     $               INFO )
      CALL ASSERT( INFO.EQ.0,
     $             '** PFGETRS: descinit for ipiv returns info ', INFO )
*
*  Perform initial permutation for
*  special case with initial double panels.
*
      IF( IC.EQ.1 ) THEN
         IP_BPIV = 1
      ELSE
         CALL INFOG1L( IC, MB, NPROW, MYPROW, RSRCIO, IIP, IPROC )
         IP_BPIV = IIP
      ENDIF
      MM = MIN( M, FIRSTSTART-1 )
      NN = JSIZEB
      CALL PCLAPIV( 'Forward', 'Rowise', 'Columwise', MM, NN, BMAT, IB,
     $              JB, DESCB, IPIV( IP_BPIV ), 1, 1, DESCIPIV, IDUMMY )
*
*   use panel A to hold factors
*
      AEND = MINMN
      JSTARTA = 1
   10 CONTINUE
      IF( JSTARTA.LE.AEND ) THEN
         ISTARTA = JSTARTA
         IENDA = M
         ISIZEA = IENDA - ISTARTA + 1
         IF( ISIZEA.LE.0 ) THEN
*
*                        Rectangular matrix.
*                        No more updates.
*
            GOTO 50
         ENDIF
         NCOLA = DESCA( N_ )
         JENDA = MIN( AEND, JSTARTA+NCOLA-1 )
         IF( JENDA.NE.AEND ) THEN
            JJC = ( JC-1 ) + JENDA
            JJC = MAX( ( JC-1 )+JSTARTA, INT( JJC / NNB )*NNB )
            JENDA = JJC - ( JC-1 )
         ENDIF
         JENDA = MAX( JSTARTA, JENDA )
         JSIZEA = JENDA - JSTARTA + 1
*
*                  Read in panel A
*
         IIC = ( IC-1 ) + ISTARTA
         JJC = ( JC-1 ) + JSTARTA
         IA = ISTARTA
         JA = 1
         CALL CLAREAD( IODEV, ISIZEA, JSIZEA, IIC, JJC,
     $                 CWORK( IP_AMAT ), IA, JA, DESCA, INFO )
         CALL ASSERT( INFO.EQ.0, '**  PFxGELUP : LAREAD returns info ',
     $                INFO )
*
*        Permute lower factors
*
         IF( JSTARTA.LE.( FIRSTSTART-1 ) ) THEN
*
*                         Special case.
*
            JJSTART = FIRSTSTART
         ELSE
            JJSTART = ( FIRSTSTART-1 ) +
     $                ICEIL( ( JSTARTA-( FIRSTSTART-1 ) ), NCOLB )*
     $                NCOLB + 1
         ENDIF
*
*                          Permute panel B to match panel A.
*
         ISALIGNED = ( 1.EQ.MOD( ( JSTARTA-( FIRSTSTART-1 ) ),
     $               NPANELS*NNB ) )
         IF( ( JSTARTA.GE.FIRSTSTART ) .AND. ISALIGNED ) THEN
*
*                         Past the first 2 panels.
*
            JJEND = MAX( JJSTART-1, JENDA )
            MM = JJEND - JSTARTA + 1
            NN = JSIZEB
            HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
            IF( HASWORK ) THEN
               IIB = ( IB-1 ) + JSTARTA
               JJB = ( JB-1 ) + 1
               CALL PCLAPIV( 'Forward', 'Rowwise', 'Columnwise', MM, NN,
     $                       BMAT, IIB, JJB, DESCB, IPIV( IP_BPIV ),
     $                       JSTARTA, 1, DESCIPIV, IDUMMY )
            ENDIF
         ENDIF
*
*  Update panel B with panel A.
*
*
*         Unroll computation to treat one mb by mb block at a time.
*
         MB = DESCA( MB_ )
         INC = MB
         IASTART = JSTARTA
   20    CONTINUE
         IF( .NOT.( IASTART.LE.JENDA ) )
     $      GOTO 40
         IAEND = MIN( JENDA, IASTART+INC-1 )
         IF( IAEND.NE.JENDA ) THEN
*
*                Align to block boundary.
*
            IAEND = MAX( IASTART, INT( IAEND / MB )*MB )
         ENDIF
         IASIZE = IAEND - IASTART + 1
         CALL ASSERT( IASIZE.GE.1, '**  PFxGELUP : iasize <= 0, iasize '
     $                , IASIZE )
         ALPHA = ONE
         MM = IASIZE
         NN = JSIZEB
         IA = IASTART
         JA = IASTART - JSTARTA + 1
         IIB = ( IB-1 ) + IASTART
         JJB = ( JB-1 ) + 1
         HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
         IF( HASWORK ) THEN
            CALL PCTRSM( 'Leftside', 'LowerTriangular', 'NoTranspose',
     $                   'UnitDiagonal', MM, NN, ALPHA,
     $                   CWORK( IP_AMAT ), IA, JA, DESCA, BMAT, IIB,
     $                   JJB, DESCB )
         ENDIF
         MM = M - ( IAEND+1 ) + 1
         NN = JSIZEB
         KK = IASIZE
         IA = IAEND + 1
         JA = IASTART - JSTARTA + 1
         IIB = ( IB-1 ) + IASTART
         JJB = ( JB-1 ) + 1
         HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 ) .AND. ( KK.GE.1 )
         IF( HASWORK ) THEN
            ALPHA = -ONE
            BETA = ONE
            CALL PCGEMM( 'Notrans', 'NoTrans', MM, NN, KK, ALPHA,
     $                   CWORK( IP_AMAT ), IA, JA, DESCA, BMAT, IIB,
     $                   JJB, DESCB, BETA, BMAT, ( IB-1 )+IA, JJB,
     $                   DESCB )
         ENDIF
   30    CONTINUE
         IASTART = IAEND + 1
         GOTO 20
   40    CONTINUE
         JSTARTA = JENDA + 1
         GOTO 10
      ENDIF
   50 CONTINUE
* end while
      RETURN
      END
