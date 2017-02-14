      SUBROUTINE PFCOCLU( IODEV, NNB, M, N, IC, JC, DESCC, IPIV, LDPIV,
     $                    WORK, LWORK, INFO )
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
* Internal auxiliary routine to
* perform out-of-core LU factorization.
*
*
*     .. Parameters ..
      INTEGER            DLEN_
      PARAMETER          ( DLEN_ = 9 )
      INTEGER            CTXT_, N_, MB_, NB_
      PARAMETER          ( CTXT_ = 2, N_ = 4, MB_ = 5, NB_ = 6 )
      INTEGER            RSRC_, CSRC_
      PARAMETER          ( RSRC_ = 7, CSRC_ = 8 )
*     ..
*     .. Scalar Arguments ..
      INTEGER            IC, INFO, IODEV, JC, LDPIV, LWORK, M, N, NNB
*     ..
*     .. Array Arguments ..
      INTEGER            DESCC( DLEN_ ), IPIV( * )
      COMPLEX            WORK( * )
*     ..
*     .. Local Scalars ..
      LOGICAL            ALLDONE, ATEND, HASWORK, SPACEOK
      INTEGER            FIRSTWIDTH, GNCOLA, GNCOLB, IA, IB, ICPOS,
     $                   IDUM, IFREE, IIC, IP_AMAT, IP_BMAT, ISPACE, JA,
     $                   JB, JCPOS, JENDA, JENDB, JJC, JSIZEA, JSIZEB,
     $                   JSTARTA, JSTARTB, KK, LDA, LDB, MINMN, MM,
     $                   MYID, MYPCOL, MYPROW, NCOLA, NCOLB, NFREE, NN,
     $                   NPCOL, NPROW, NROW, SBWID
      COMPLEX            ALPHA, BETA, ONE
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA( DLEN_ ), DESCB( DLEN_ ),
     $                   DESCIPIV( DLEN_ )
*     ..
*     .. External Functions ..
      INTEGER            BLACS_PNUM, NUMROC
      EXTERNAL           BLACS_PNUM, NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           ASSERT, BLACS_GRIDINFO, CALCOLSIZE, CLAREAD,
     $                   CLAWRITE, DESCINIT, PCGEMM, PCGETRF, PCLAPIV,
     $                   PCTRSM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, MAX, MIN, MOD, REAL
*     ..
*     .. Executable Statements ..
*
*       Get information on processor mesh
*
      CALL BLACS_GRIDINFO( DESCC( CTXT_ ), NPROW, NPCOL, MYPROW,
     $                     MYPCOL )
      MYID = BLACS_PNUM( DESCC( CTXT_ ), MYPROW, MYPCOL )
      ONE = CMPLX( REAL( 1 ) )
      MINMN = MIN( M, N )
*
* Initialize descriptors
*
      NROW = NUMROC( MINMN, DESCC( MB_ ), MYPROW, DESCC( RSRC_ ),
     $       NPROW )
      ISPACE = NROW + DESCC( MB_ )*NPROW
      SPACEOK = ( ISPACE.LE.LDPIV )
      CALL ASSERT( SPACEOK, '** PFOCLU: increase ldpiv to ', ISPACE )
      SPACEOK = ( 1.LE.LDPIV ) .AND. ( DESCC( MB_ )+NROW.LE.LDPIV )
      CALL ASSERT( SPACEOK, 'increase ldpiv to ', DESCC( MB_ )+NROW )
*
* Note the use of mypcol as source processor, not descC( CSRC_)
*
      CALL DESCINIT( DESCIPIV, M+NPROW*DESCC( MB_ ), 1, DESCC( MB_ ), 1,
     $               DESCC( RSRC_ ), MYPCOL, DESCC( CTXT_ ), LDPIV,
     $               INFO )
      CALL ASSERT( INFO.EQ.0, '** PFOCLU: descinit returns info ',
     $             INFO )
*
* Descriptor for panel A
*
      GNCOLA = N
      NROW = NUMROC( M, DESCC( MB_ ), MYPROW, DESCC( RSRC_ ), NPROW )
      LDA = MAX( 1, NROW )
      CALL DESCINIT( DESCA, M, GNCOLA, DESCC( MB_ ), DESCC( NB_ ),
     $               DESCC( RSRC_ ), DESCC( CSRC_ ), DESCC( CTXT_ ),
     $               LDA, INFO )
      CALL ASSERT( INFO.EQ.0, '** PFOCLU: descinit returns info ',
     $             INFO )
*
* Descriptor for panel B
*
      GNCOLB = N
      LDB = MAX( 1, NROW )
      CALL DESCINIT( DESCB, M, GNCOLB, DESCC( MB_ ), DESCC( NB_ ),
     $               DESCC( RSRC_ ), DESCC( CSRC_ ), DESCC( CTXT_ ),
     $               LDB, INFO )
      CALL ASSERT( INFO.EQ.0, '** PFOCLU: descinit returns info ',
     $             INFO )
      IFREE = 1
      NFREE = LWORK
      CALL CALCOLSIZE( NFREE, M, DESCB, GNCOLB )
      GNCOLB = MAX( 1, MIN( N, GNCOLB ) )
      ATEND = ( GNCOLB.EQ.N )
      IF( .NOT.ATEND ) THEN
*
*       Align to block boundary
*
         IF( GNCOLB.GE.NNB ) THEN
            GNCOLB = GNCOLB - MOD( GNCOLB, NNB )
         ENDIF
         IF( GNCOLB.GE.DESCC( NB_ ) ) THEN
            GNCOLB = GNCOLB - MOD( GNCOLB, DESCC( NB_ ) )
         ENDIF
      ENDIF
      DESCB( N_ ) = GNCOLB
      NCOLB = NUMROC( GNCOLB, DESCC( NB_ ), MYPCOL, DESCC( CSRC_ ),
     $        NPCOL )
      IP_BMAT = IFREE
      IFREE = IFREE + NCOLB*NROW
      NFREE = LWORK - IFREE + 1
      SPACEOK = ( NFREE.GE.0 )
      CALL ASSERT( SPACEOK, '** PFOCLU: insufficient work space, need ',
     $             ( NCOLB*NROW ) )
*
* First panel is special case.
* Use all of available storage.
*
      CALL CLAREAD( IODEV, M, GNCOLB, IC, JC, WORK( IP_BMAT ), 1, 1,
     $              DESCB, INFO )
      CALL ASSERT( INFO.EQ.0, '** PFOCLU: dlread returns info ', INFO )
      CALL PCGETRF( M, GNCOLB, WORK( IP_BMAT ), 1, 1, DESCB, IPIV,
     $              INFO )
      CALL ASSERT( INFO.EQ.0, '** PFOCLU: PGETRF returns info ', INFO )
      ALLDONE = ( GNCOLB.EQ.N )
      IF( .NOT.ALLDONE ) THEN
*
*       Need to undo permutation
*       before writing out to disk.
*
         MM = MIN( GNCOLB, MINMN )
         NN = GNCOLB
         HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
         IF( HASWORK ) THEN
            CALL PCLAPIV( 'B', 'R', 'C', MM, NN, WORK( IP_BMAT ), 1, 1,
     $                    DESCB, IPIV, 1, 1, DESCIPIV, IDUM )
         ENDIF
      ENDIF
      CALL CLAWRITE( IODEV, M, GNCOLB, IC, JC, WORK( IP_BMAT ), 1, 1,
     $               DESCB, INFO )
      CALL ASSERT( INFO.EQ.0, '** PFOCLU: dlwrite returns info ', INFO )
      IF( ALLDONE ) THEN
         INFO = 0
         RETURN
      ENDIF
      FIRSTWIDTH = GNCOLB
*
*       Allocate space for panels
*
      IFREE = 1
      NFREE = LWORK
      GNCOLA = NNB
      DESCA( N_ ) = GNCOLA
      NCOLA = NUMROC( GNCOLA, DESCC( NB_ ), MYPCOL, DESCC( CSRC_ ),
     $        NPCOL )
*
*       Allocate storage for panel A
*
      IP_AMAT = IFREE
      IFREE = IFREE + ( NROW*NCOLA )
      NFREE = LWORK - IFREE + 1
      CALL ASSERT( IFREE.GE.0,
     $             '** PFOCLU: insufficient storage, require extra ',
     $             NFREE )
*
*       Allocate the rest as storage for panel B
*
      CALL CALCOLSIZE( NFREE, M, DESCB, GNCOLB )
      IF( GNCOLB.GE.NNB ) THEN
         GNCOLB = GNCOLB - MOD( GNCOLB, NNB )
      ENDIF
      IF( GNCOLB.GE.DESCC( NB_ ) ) THEN
         GNCOLB = GNCOLB - MOD( GNCOLB, DESCC( NB_ ) )
      ENDIF
      DESCB( N_ ) = GNCOLB
      CALL ASSERT( GNCOLB.GE.DESCB( NB_ ),
     $             '** PFOCLU: insufficient storage, gncolB ', GNCOLB )
      NCOLB = NUMROC( GNCOLB, DESCC( NB_ ), MYPCOL, DESCC( CSRC_ ),
     $        NPCOL )
      IP_BMAT = IFREE
      IFREE = IFREE + NROW*NCOLB
      NFREE = LWORK - IFREE + 1
      CALL ASSERT( IFREE.GE.0,
     $             '** PFOCLU: insufficient storage, require extra ',
     $             NFREE )
*
*       Loop over column superblocks
*
      SBWID = GNCOLB
      DO 30 JSTARTB = FIRSTWIDTH + 1, N, SBWID
         JENDB = MIN( N, JSTARTB+SBWID-1 )
         JSIZEB = JENDB - JSTARTB + 1
*
*  Read in current block
*
         ICPOS = ( IC-1 ) + 1
         JCPOS = ( JC-1 ) + JSTARTB
         CALL CLAREAD( IODEV, M, JSIZEB, ICPOS, JCPOS, WORK( IP_BMAT ),
     $                 1, 1, DESCB, INFO )
         CALL ASSERT( INFO.EQ.0, '** PFOCLU: laread returns info ',
     $                INFO )
*
*  Apply previous interchanges to current block
*
         MM = MIN( MINMN, ( JSTARTB-1 ) )
         NN = JSIZEB
         HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
         IF( HASWORK ) THEN
            CALL PCLAPIV( 'F', 'R', 'C', MM, NN, WORK( IP_BMAT ), 1, 1,
     $                    DESCB, IPIV, 1, 1, DESCIPIV, IDUM )
         ENDIF
*
*  Loop over superblocks to left of current superblock/panel
*
         DO 10 JSTARTA = 1, MIN( MINMN, JSTARTB-1 ), DESCA( N_ )
            JENDA = MIN( MIN( MINMN, JSTARTB-1 ),
     $              JSTARTA+DESCA( N_ )-1 )
            JSIZEA = JENDA - JSTARTA + 1
*
* Read in temporary superblock (panel A)
*
            ICPOS = ( IC-1 ) + 1
            JCPOS = ( JC-1 ) + JSTARTA
            MM = M
            NN = JSIZEA
            HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
            IF( HASWORK ) THEN
               CALL CLAREAD( IODEV, M, JSIZEA, ICPOS, JCPOS,
     $                       WORK( IP_AMAT ), 1, 1, DESCA, INFO )
               CALL ASSERT( INFO.EQ.0, '** PFOCLU: laread returns info '
     $                      , INFO )
            ENDIF
*
*  Apply interchanges to temporary superblock (panel A)
*
            MM = MIN( MINMN, ( JSTARTB-1 ) )
            NN = JSIZEA
            HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
            IF( HASWORK ) THEN
               CALL PCLAPIV( 'F', 'R', 'C', MM, NN, WORK( IP_AMAT ), 1,
     $                       1, DESCA, IPIV, 1, 1, DESCIPIV, IDUM )
            ENDIF
*
* Determine sub-block K of current superblock
*
*  L11*U12 = A12,   L11*U13 = A13
*
* [ A22 ] <- [ A22 ] - [ L21 ]*U12
* [ A23 ] <- [ A23 ]   [ L23 ]
*
            MM = JSIZEA
            NN = JSIZEB
            ALPHA = ONE
            IA = JSTARTA
            JA = 1
            IB = JSTARTA
            JB = 1
*
*  Solve L11*U12 = A12
*
            HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
            IF( HASWORK ) THEN
               CALL PCTRSM( 'Left', 'Lower', 'No transpose', 'Unit', MM,
     $                      NN, ALPHA, WORK( IP_AMAT ), IA, JA, DESCA,
     $                      WORK( IP_BMAT ), IB, JB, DESCB )
            ENDIF
*
*  Update sub-blocks below sub-block K of
*  current superblock (U12)
*
            MM = M - ( JENDA+1 ) + 1
            NN = JSIZEB
            KK = JSIZEA
            ALPHA = -ONE
            BETA = ONE
            IA = JENDA + 1
            JA = 1
            IB = JSTARTA
            JB = 1
            IIC = JENDA + 1
            JJC = 1
            HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 ) .AND. ( KK.GE.1 )
            IF( HASWORK ) THEN
               CALL PCGEMM( 'No transpose', 'No transpose', MM, NN, KK,
     $                      ALPHA, WORK( IP_AMAT ), IA, JA, DESCA,
     $                      WORK( IP_BMAT ), IB, JB, DESCB, BETA,
     $                      WORK( IP_BMAT ), IIC, JJC, DESCB )
            ENDIF
   10    CONTINUE
   20    CONTINUE
* end do jstartA
* ======================================
* all previous updates are done
* factor diagonal and subdiagonal blocks and
* test for exact singularity
* ========================================
         MM = M - JSTARTB + 1
         NN = JSIZEB
         IB = JSTARTB
         JB = 1
         HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
         IF( HASWORK ) THEN
            CALL PCGETRF( MM, NN, WORK( IP_BMAT ), IB, JB, DESCB, IPIV,
     $                    INFO )
            CALL ASSERT( INFO.EQ.0, '** PFOCLU: PGETRF returns info ',
     $                   INFO )
         ENDIF
* ------------------------------------
* undo interchanges for current block
* ------------------------------------
         MM = MIN( MINMN, JENDB )
         NN = JSIZEB
         HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
         IF( HASWORK ) THEN
            CALL PCLAPIV( 'B', 'R', 'C', MM, NN, WORK( IP_BMAT ), 1, 1,
     $                    DESCB, IPIV, 1, 1, DESCIPIV, IDUM )
         ENDIF
* ---------------------------
* write current block to file
* ---------------------------
         MM = M
         NN = JSIZEB
         ICPOS = ( IC-1 ) + 1
         JCPOS = ( JC-1 ) + JSTARTB
         HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
         IF( HASWORK ) THEN
            CALL CLAWRITE( IODEV, MM, NN, ICPOS, JCPOS, WORK( IP_BMAT ),
     $                     1, 1, DESCB, INFO )
            CALL ASSERT( INFO.EQ.0, '** PFOCLU: lawrite returns info ',
     $                   INFO )
         ENDIF
   30 CONTINUE
   40 CONTINUE
* end do jstartB
* =======================================
* Factorization is complete but still need to
* apply interchanges to matrix columns 1:min(m,n)
* =======================================
*
*    Repartition and use all memory
*    also try to align to i/o boundary
*
      IFREE = 1
      NFREE = LWORK
      CALL CALCOLSIZE( NFREE, M, DESCB, GNCOLB )
      GNCOLB = MAX( 1, MIN( N, GNCOLB ) )
      IF( GNCOLB.GE.NNB ) THEN
         GNCOLB = GNCOLB - MOD( GNCOLB, NNB )
      ENDIF
      IF( GNCOLB.GE.DESCB( NB_ ) ) THEN
         GNCOLB = GNCOLB - MOD( GNCOLB, DESCB( NB_ ) )
      ENDIF
      NROW = NUMROC( M, DESCB( MB_ ), MYPROW, DESCB( RSRC_ ), NPROW )
      NCOLB = NUMROC( GNCOLB, DESCB( NB_ ), MYPCOL, DESCB( CSRC_ ),
     $        NPCOL )
      IP_BMAT = IFREE
      IFREE = IFREE + NROW*NCOLB
      NFREE = LWORK - IFREE + 1
      CALL ASSERT( NFREE.GE.0, '** PFOCLU: internal error, nfree ',
     $             NFREE )
      DESCB( N_ ) = GNCOLB
      SBWID = GNCOLB
      DO 50 JSTARTB = 1, N, SBWID
         JENDB = MIN( N, JSTARTB+SBWID-1 )
         JSIZEB = JENDB - JSTARTB + 1
*
*  Read in panel
*
         MM = M
         NN = JSIZEB
         ICPOS = ( IC-1 ) + 1
         JCPOS = ( JC-1 ) + JSTARTB
         HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
         IF( HASWORK ) THEN
            CALL CLAREAD( IODEV, MM, NN, ICPOS, JCPOS, WORK( IP_BMAT ),
     $                    1, 1, DESCB, INFO )
            CALL ASSERT( INFO.EQ.0, '** PFOCLU: laread returns info ',
     $                   INFO )
         ENDIF
*
* Perform permutation
*
         MM = MIN( M, N )
         NN = JSIZEB
         HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
         IF( HASWORK ) THEN
            CALL PCLAPIV( 'F', 'R', 'C', MM, NN, WORK( IP_BMAT ), 1, 1,
     $                    DESCB, IPIV, 1, 1, DESCIPIV, IDUM )
         ENDIF
*
*  Write out panel
*
         MM = M
         NN = JSIZEB
         ICPOS = ( IC-1 ) + 1
         JCPOS = ( JC-1 ) + JSTARTB
         HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
         IF( HASWORK ) THEN
            CALL CLAWRITE( IODEV, MM, NN, ICPOS, JCPOS, WORK( IP_BMAT ),
     $                     1, 1, DESCB, INFO )
            CALL ASSERT( INFO.EQ.0, '** PFOCLU: lawrite returns info ',
     $                   INFO )
         ENDIF
   50 CONTINUE
   60 CONTINUE
*
* All done
*
      INFO = 0
      RETURN
      END
