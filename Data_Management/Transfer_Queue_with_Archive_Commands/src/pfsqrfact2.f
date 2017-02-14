      SUBROUTINE PFSQRFACT2( IODEV, M, N, CWORK, IC, JC, DESCC, TAUC,
     $                       WORK, LWORK )
*
*
*  -- ScaLAPACK auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     Oct 10, 1996
*
*
* Purpose:
* ========
*
* Perform out of core QR factorization.
*
* A, B are temporary scalapack arrays used in the factorization.
*
* Matrix A is used in the application of Householder transformations.
*
* Matrix B is used to accumulated updates and perform
* in-core QR factorization.
*
*
*     .. Parameters ..
      INTEGER            DLEN_
      PARAMETER          ( DLEN_ = 9 )
      INTEGER            CTXT_
      PARAMETER          ( CTXT_ = 2 )
      INTEGER            RSRC_, CSRC_
      PARAMETER          ( RSRC_ = 7, CSRC_ = 8 )
      INTEGER            IODEV_, SIZE_
      PARAMETER          ( IODEV_ = 10, SIZE_ = 11 )
      INTEGER            FLEN_
      PARAMETER          ( FLEN_ = 12 )
*     ..
*     .. Scalar Arguments ..
      INTEGER            IC, IODEV, JC, LWORK, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCC( FLEN_ )
      REAL               CWORK( * ), TAUC( * ), WORK( LWORK )
*     ..
*     .. Local Scalars ..
      LOGICAL            ALL_INVOLVED, HASWORK, ISVALID
      INTEGER            AEND, ANEED, BNEED, CSIZE, CSRCIO, CWORKNEED,
     $                   FIRSTSTART, IAPOS, IBPOS, ICOLA, ICOLC,
     $                   ICONTXT, ICPOS, IENDA, IENDB, IFREE, INFO,
     $                   IP_AMAT, IP_BMAT, IP_TAU, IP_TAUA, IP_TAUB,
     $                   IROWA, IROWB, IROWC, ISIZEA, ISIZEB, ISTARTA,
     $                   ISTARTB, JAPOS, JBPOS, JCPOS, JENDA, JENDB,
     $                   JSIZEA, JSIZEB, JSTARTA, JSTARTB, KK, LDA, LDB,
     $                   LINFO, LOCP, LOCQ, MB, MM, MMB, MYPCOL, MYPROW,
     $                   NB, NCOLA, NCOLB, NFREE, NN, NNB, NPCOL, NPROW,
     $                   P0, Q0, QR_NEED, RSRCIO, TAUA_DIM, TAUB_DIM,
     $                   TAU_NEED, WORKNEED
      REAL               DZERO
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA( DLEN_ ), DESCB( DLEN_ )
*     ..
*     .. External Functions ..
      INTEGER            INDXG2P, NUMROC
      EXTERNAL           INDXG2P, NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           ASSERT, BLACS_GRIDINFO, DESCINIT, LAIO_INFO,
     $                   PFMAXSIZE, PFSCOPYTAU, PSGEQRF, PSORMQR, SFILL,
     $                   SGEBR2D, SGEBS2D, SLAREAD, SLAWRITE
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT, MAX, MIN, REAL
*     ..
*     .. Executable Statements ..
      DZERO = REAL( 0 )
      IODEV = DESCC( IODEV_ )
      CALL LAIO_INFO( IODEV, MM, NN, MMB, NNB, MB, NB, CSRCIO, RSRCIO,
     $                ICONTXT )
      CALL BLACS_GRIDINFO( ICONTXT, NPROW, NPCOL, MYPROW, MYPCOL )
      IFREE = 1
      IP_BMAT = 1
*
*         First pass, special case,
*         use all memory for panel.
*
      CSIZE = DESCC( SIZE_ )
      IROWC = INDXG2P( IC, MB, MYPROW, RSRCIO, NPROW )
      ICOLC = INDXG2P( JC, NB, MYPCOL, CSRCIO, NPCOL )
      P0 = IROWC
      Q0 = ICOLC
      NCOLB = -1
      CALL PFMAXSIZE( 'Column', CSIZE, M, NCOLB, MB, NB, P0, Q0,
     $                ICONTXT, INFO )
      CALL ASSERT( INFO.EQ.0, 'PFxQRFACT2 : pfmaxsize returns info ',
     $             INFO )
*
*        Align to block boundary.
*
      JSTARTB = 1
      JENDB = MIN( N, NCOLB )
      IF( JENDB.NE.N ) THEN
         JCPOS = ( JC-1 ) + JENDB
         JCPOS = MAX( JC, INT( JCPOS / NNB )*NNB )
         JENDB = JCPOS - ( JC-1 )
         JENDB = MAX( JSTARTB, JENDB )
      ENDIF
      JSIZEB = JENDB - JSTARTB + 1
      P0 = IROWC
      Q0 = ICOLC
      LOCP = NUMROC( M, MB, MYPROW, P0, NPROW )
      LOCQ = NUMROC( JSIZEB, NB, MYPCOL, Q0, NPCOL )
      LDB = MAX( 1, LOCP )
      CWORKNEED = LOCP*LOCQ
      ISVALID = ( LDB*LOCQ.LE.CSIZE )
      CALL ASSERT( ISVALID, 'PFxQRFACT2 : ldB*Locq > Csize, ldB*Locq = '
     $             , LDB*LOCQ )
      TAUB_DIM = LOCQ
      IF( TAUB_DIM.GT.LWORK ) THEN
         INFO = -10
         WORK( 1 ) = TAUB_DIM
         CWORK( 1 ) = CWORKNEED
         RETURN
      ENDIF
      IP_TAUB = IFREE
      IFREE = IFREE + LOCQ
      NFREE = LWORK - IFREE + 1
      IP_TAUA = IP_TAUB
      TAU_NEED = LOCQ
      CALL DESCINIT( DESCB, M, JSIZEB, MB, NB, P0, Q0, ICONTXT, LDB,
     $               INFO )
      CALL ASSERT( INFO.EQ.0, 'PFxQRFACT2 : descinit returns info ',
     $             INFO )
      ICPOS = ( IC-1 ) + 1
      JCPOS = ( JC-1 ) + 1
      CALL SLAREAD( IODEV, M, JSIZEB, ICPOS, JCPOS, CWORK( IP_BMAT ), 1,
     $              1, DESCB, INFO )
      CALL ASSERT( INFO.EQ.0, 'PFxQRFACT2 : LAREAD returns info ',
     $             INFO )
      LINFO = 0
      CALL PSGEQRF( M, JSIZEB, CWORK( IP_BMAT ), 1, 1, DESCB,
     $              WORK( IP_TAUB ), WORK( IFREE ), NFREE, LINFO )
*
*        If only a subset of processors are involved,
*        need to perform broadcast of tau vector.
*
      ALL_INVOLVED = ( M.GE.NPROW*MB )
      IF( .NOT.ALL_INVOLVED ) THEN
         LOCQ = NUMROC( JSIZEB, NB, MYPCOL, DESCB( CSRC_ ), NPCOL )
         HASWORK = ( LOCQ.GE.1 )
         IF( HASWORK ) THEN
            IF( MYPROW.EQ.DESCB( RSRC_ ) ) THEN
               CALL SGEBS2D( DESCB( CTXT_ ), 'Column', 'I', LOCQ, 1,
     $                       WORK( IP_TAUB ), LOCQ )
            ELSE
               CALL SGEBR2D( DESCB( CTXT_ ), 'Column', 'I', LOCQ, 1,
     $                       WORK( IP_TAUB ), LOCQ, DESCB( RSRC_ ),
     $                       MYPCOL )
            ENDIF
         ENDIF
      ENDIF
      JCPOS = ( JC-1 ) + JSTARTB
      CALL PFSCOPYTAU( JSIZEB, 1, WORK( IP_TAUB ), DESCB, JCPOS, TAUC,
     $                 DESCC )
      QR_NEED = INT( WORK( IFREE ) )
      WORKNEED = QR_NEED + TAU_NEED
      IF( LINFO.GE.1 ) THEN
         INFO = LINFO
         WORK( 1 ) = WORKNEED
         CWORK( 1 ) = CWORKNEED
         RETURN
      ENDIF
      ICPOS = ( IC-1 ) + 1
      JCPOS = ( JC-1 ) + 1
      CALL SLAWRITE( IODEV, M, JSIZEB, ICPOS, JCPOS, CWORK( IP_BMAT ),
     $               1, 1, DESCB, INFO )
      CALL ASSERT( INFO.EQ.0, 'PFxQRFACT2 : LAWRITE returns info ',
     $             INFO )
      FIRSTSTART = JENDB + 1
      IF( FIRSTSTART.GE.N ) THEN
*
*          The matrix fit in core. All done.
*
         INFO = 0
         WORK( 1 ) = WORKNEED
         CWORK( 1 ) = CWORKNEED
         RETURN
      ENDIF
*
*         Panel A need one m * nnb pannel,
*         the rest is storage  for panel B.
*
      IROWC = INDXG2P( IC, MB, MYPROW, RSRCIO, NPROW )
      ICOLC = INDXG2P( FIRSTSTART, NB, MYPCOL, CSRCIO, NPCOL )
      P0 = IROWC
      Q0 = ICOLC
      NCOLA = MIN( NNB, N )
      LOCP = NUMROC( M, MB, MYPROW, P0, NPROW )
      LOCQ = NUMROC( NCOLA, NB, MYPCOL, Q0, NPCOL )
      LDA = MAX( LOCP, 1 )
      TAUA_DIM = LOCQ
      ANEED = LOCP*LOCQ
      IP_AMAT = 1
      IP_BMAT = IP_AMAT + ANEED
      CALL PFMAXSIZE( 'Column', CSIZE-ANEED, M, NCOLB, MB, NB, P0, Q0,
     $                ICONTXT, INFO )
      CALL ASSERT( INFO.EQ.0, 'PFxQRFACT2 : pfmaxsize returns info ',
     $             INFO )
      LOCP = NUMROC( M, MB, MYPROW, P0, NPROW )
      LOCQ = NUMROC( MIN( N, NNB ), NB, MYPCOL, Q0, NPCOL )
      BNEED = LOCP*LOCQ
      TAUB_DIM = LOCQ
      LOCQ = NUMROC( NCOLB, NB, MYPCOL, Q0, NPCOL )
      LDB = MAX( 1, LOCQ )
      BNEED = MAX( BNEED, LOCP*LOCQ )
      TAUB_DIM = MAX( TAUB_DIM, LOCQ )
      TAU_NEED = MAX( TAUA_DIM, TAUB_DIM )
      IP_TAU = IFREE
      IFREE = IFREE + TAU_NEED
      NFREE = LWORK - IFREE + 1
      ISVALID = ( ANEED+BNEED.LE.CSIZE )
      IF( .NOT.ISVALID ) THEN
         INFO = -( 7*100+SIZE_ )
         CWORK( 1 ) = MAX( CWORKNEED, ANEED+BNEED )
         WORK( 1 ) = MAX( WORKNEED, TAU_NEED )
         RETURN
      ENDIF
      ISVALID = ( TAU_NEED.LT.LWORK )
      IF( .NOT.ISVALID ) THEN
         INFO = -10
         CWORK( 1 ) = MAX( CWORKNEED, ANEED+BNEED )
         WORK( 1 ) = MAX( WORKNEED, TAU_NEED )
         RETURN
      ENDIF
*
*        Reuse storage.
*
      IP_TAUA = IP_TAU
      IP_TAUB = IP_TAU
      JSTARTB = FIRSTSTART
   10 CONTINUE
      IF( JSTARTB.LE.N ) THEN
         JENDB = MIN( N, JSTARTB+NCOLB-1 )
         IF( JENDB.NE.N ) THEN
*
*                Align to block boundary.
*
            JCPOS = ( JC-1 ) + JENDB
            JCPOS = MAX( ( JC-1 )+JSTARTB, INT( JCPOS / NNB )*NNB )
            JENDB = JCPOS - ( JC-1 )
            JENDB = MAX( JSTARTB, JENDB )
         ENDIF
         JSIZEB = JENDB - JSTARTB + 1
         ICPOS = ( IC-1 ) + 1
         JCPOS = ( JC-1 ) + JSTARTB
         IROWC = INDXG2P( ICPOS, MB, MYPROW, RSRCIO, NPROW )
         ICOLC = INDXG2P( JCPOS, NB, MYPCOL, CSRCIO, NPCOL )
         P0 = IROWC
         Q0 = ICOLC
         LOCP = NUMROC( M, MB, MYPROW, P0, NPROW )
         LDB = MAX( 1, LOCP )
         LOCQ = NUMROC( JSIZEB, NB, MYPCOL, Q0, NPCOL )
         ISVALID = ( LOCP*LOCQ.LE.BNEED )
         CALL ASSERT( ISVALID,
     $          'PFxQRFACT2 : internal error, Locp*Locq > Bneed, Bneed '
     $                , BNEED )
         CALL DESCINIT( DESCB, M, JSIZEB, MB, NB, P0, Q0, ICONTXT, LDB,
     $                  INFO )
         LOCP = NUMROC( M, MB, MYPROW, P0, NPROW )
         LDA = MAX( 1, LOCP )
         CALL DESCINIT( DESCA, M, NCOLA, MB, NB, P0, Q0, ICONTXT, LDA,
     $                  INFO )
*
*           Read in panel B.
*
         ICPOS = ( IC-1 ) + 1
         JCPOS = ( JC-1 ) + JSTARTB
         INFO = 0
         CALL SLAREAD( IODEV, M, JSIZEB, ICPOS, JCPOS, CWORK( IP_BMAT ),
     $                 1, 1, DESCB, INFO )
         CALL ASSERT( INFO.EQ.0,
     $                'PFxQRFACT2 : laread for B returns info ', INFO )
         JSTARTA = 1
         AEND = MIN( MIN( M, N ), JSTARTB-1 )
   20    CONTINUE
         IF( JSTARTA.LE.AEND ) THEN
            JENDA = MIN( JSTARTA+NCOLA-1, AEND )
            IF( JENDA.NE.AEND ) THEN
               JCPOS = ( JC-1 ) + JENDA
               JCPOS = MAX( ( JC-1 )+JSTARTA, INT( JCPOS / NNB )*NNB )
               JENDA = JCPOS - ( JC-1 )
               JENDA = MAX( JSTARTA, JENDA )
            ENDIF
            JSIZEA = JENDA - JSTARTA + 1
            ICPOS = ( IC-1 ) + JSTARTA
            JCPOS = ( JC-1 ) + JSTARTA
            IAPOS = JSTARTA
            JAPOS = 1
*
*                Double check alignment, should
*                perform local I/O.
*
            IROWC = INDXG2P( ICPOS, MB, MYPROW, RSRCIO, NPROW )
            ICOLC = INDXG2P( JCPOS, NB, MYPCOL, CSRCIO, NPCOL )
            IROWA = INDXG2P( IAPOS, MB, MYPROW, DESCA( RSRC_ ), NPROW )
            ICOLA = INDXG2P( JAPOS, NB, MYPCOL, DESCA( CSRC_ ), NPCOL )
            ISVALID = ( IROWC.EQ.IROWA ) .AND. ( ICOLC.EQ.ICOLA )
            CALL ASSERT( ISVALID,
     $                   'PFxQRFACT2 : misaligned A, jstartA = ',
     $                   JSTARTA )
            ISTARTA = JSTARTA
            IENDA = M
            ISIZEA = IENDA - ISTARTA + 1
            INFO = 0
            CALL SLAREAD( IODEV, ISIZEA, JSIZEA, ICPOS, JCPOS,
     $                    CWORK( IP_AMAT ), IAPOS, JAPOS, DESCA, INFO )
            CALL ASSERT( INFO.EQ.0,
     $                   'PFxQRFACT2 : laread for A returns info ',
     $                   INFO )
            CALL SFILL( TAUA_DIM, DZERO, WORK( IP_TAUA ), 1 )
            CALL PFSCOPYTAU( JSIZEA, JCPOS, TAUC, DESCC, JAPOS,
     $                       WORK( IP_TAUA ), DESCA )
*
*                Apply Householder transformation to
*                the lower part of panel B.
*
            HASWORK = ( ISIZEA.GE.1 ) .AND. ( JSIZEA.GE.1 ) .AND.
     $                ( JSIZEB.GE.1 )
            IF( HASWORK ) THEN
               MM = ISIZEA
               NN = JSIZEB
               KK = JSIZEA
               IBPOS = IAPOS
               JBPOS = 1
               LINFO = 0
               CALL PSORMQR( 'Leftside', 'T', MM, NN, KK,
     $                       CWORK( IP_AMAT ), IAPOS, JAPOS, DESCA,
     $                       WORK( IP_TAUA ), CWORK( IP_BMAT ), IBPOS,
     $                       JBPOS, DESCB, WORK( IFREE ), NFREE, LINFO )
               QR_NEED = MAX( QR_NEED, INT( WORK( IFREE ) ) )
               IF( LINFO.NE.0 ) THEN
                  IF( LINFO.EQ.16 ) THEN
                     INFO = -10
                     CWORK( 1 ) = MAX( CWORKNEED, ANEED+BNEED )
                     WORK( 1 ) = MAX( QR_NEED+TAU_NEED, WORKNEED )
                     RETURN
                  ENDIF
               ENDIF
               CALL ASSERT( LINFO.EQ.0,
     $                      'PFxQRFACT2 : PORMQR returns info = ',
     $                      LINFO )
            ENDIF
            JSTARTA = JENDA + 1
            GOTO 20
         ENDIF
   30    CONTINUE
*
*               All updates are done.
*               Perform QR factorization.
*
*
*         Matrix may be rectangular. n > m
*
         ISTARTB = JSTARTB
         IENDB = M
         ISIZEB = IENDB - ISTARTB + 1
         MM = ISIZEB
         NN = JSIZEB
         HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
         IF( HASWORK ) THEN
            IBPOS = ISTARTB
            JBPOS = 1
            LINFO = 0
            CALL PSGEQRF( MM, NN, CWORK( IP_BMAT ), IBPOS, JBPOS, DESCB,
     $                    WORK( IP_TAUB ), WORK( IFREE ), NFREE, LINFO )
            QR_NEED = MAX( QR_NEED, INT( WORK( IFREE ) ) )
            IF( LINFO.NE.0 ) THEN
               IF( LINFO.GE.1 ) THEN
                  INFO = ( JSTARTB-1 ) + LINFO
               ELSE
                  INFO = -10
               ENDIF
               WORK( 1 ) = MAX( WORKNEED, TAU_NEED+QR_NEED )
               CWORK( 1 ) = MAX( CWORKNEED, ANEED+BNEED )
               RETURN
            ENDIF
*
*         If only a subset of processors are involved,
*         need to perform broadcast of tau vector.
*
            ALL_INVOLVED = ( MM.GE.NPROW*MB )
            IF( .NOT.ALL_INVOLVED ) THEN
               LOCQ = NUMROC( NN, NB, MYPCOL, DESCB( CSRC_ ), NPCOL )
               HASWORK = ( LOCQ.GE.1 )
               IF( HASWORK ) THEN
                  IROWB = INDXG2P( IBPOS, MB, MYPROW, DESCB( RSRC_ ),
     $                    NPROW )
                  IF( MYPROW.EQ.IROWB ) THEN
                     CALL SGEBS2D( DESCB( CTXT_ ), 'Column', 'I', LOCQ,
     $                             1, WORK( IP_TAUB ), LOCQ )
                  ELSE
                     CALL SGEBR2D( DESCB( CTXT_ ), 'Column', 'I', LOCQ,
     $                             1, WORK( IP_TAUB ), LOCQ, IROWB,
     $                             MYPCOL )
                  ENDIF
               ENDIF
            ENDIF
            JCPOS = ( JC-1 ) + JSTARTB
            CALL PFSCOPYTAU( JSIZEB, 1, WORK( IP_TAUB ), DESCB, JCPOS,
     $                       TAUC, DESCC )
         ENDIF
*
*            Write out panel B.
*
         ICPOS = ( IC-1 ) + 1
         JCPOS = ( JC-1 ) + JSTARTB
         INFO = 0
         CALL SLAWRITE( IODEV, M, JSIZEB, ICPOS, JCPOS,
     $                  CWORK( IP_BMAT ), 1, 1, DESCB, INFO )
         CALL ASSERT( INFO.EQ.0,
     $                'PFxQRFACT2 : lawrite for B returns info ', INFO )
         JSTARTB = JENDB + 1
         GOTO 10
      ENDIF
   40 CONTINUE
      CWORK( 1 ) = MAX( CWORKNEED, ANEED+BNEED )
      WORK( 1 ) = MAX( WORKNEED, QR_NEED+TAU_NEED )
      INFO = 0
      RETURN
      END
