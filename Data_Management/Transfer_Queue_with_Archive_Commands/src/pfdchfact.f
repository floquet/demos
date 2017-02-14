      SUBROUTINE PFDCHFACT( IODEV, NNB, M, IC, JC, DESCC, WORK, LWORK,
     $                      INFO )
*
*
*  -- ScaLAPACK auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     Oct 10, 1996
*
*
* Purpose
* =======
*
* Auxiliary routine to perform out-of-core Cholesky factorization.
* This routine works on the lower half of matrix but
* may require disk storage for full matrix.
*
* A 'Left-look' variant of the Cholesky factorization
* is implemented.
*
* The routine attempts to use a variable width column panel to
* fully utilize all of temporary workspace in performing
* the factorization.
*
*
*
*     .. Parameters ..
      INTEGER            DLEN_
      PARAMETER          ( DLEN_ = 9 )
      INTEGER            CTXT_, M_, N_, MB_, NB_
      PARAMETER          ( CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6 )
      INTEGER            RSRC_, CSRC_
      PARAMETER          ( RSRC_ = 7, CSRC_ = 8 )
      LOGICAL            USE_PTRSM
      PARAMETER          ( USE_PTRSM = .false. )
*     ..
*     .. Scalar Arguments ..
      INTEGER            IC, INFO, IODEV, JC, LWORK, M, NNB
*     ..
*     .. Array Arguments ..
      INTEGER            DESCC( DLEN_ )
      DOUBLE PRECISION   WORK( * )
*     ..
*     .. Local Scalars ..
      LOGICAL            HASWORK, ISVALID
      CHARACTER          DIAG, SIDE, TRANS, UPLO
      INTEGER            ASIZE, BSIZE, GNCOLA, GNCOLB, IA, IA2, IB,
     $                   ICPOS, IENDA, IENDB, IFREE, IP_AMAT, IP_BMAT,
     $                   IP_SAVE, ISIZEA, ISIZEB, ISTARTA, ISTARTB, JA,
     $                   JA2, JB, JCPOS, JENDA, JENDB, JJEND, JJINC,
     $                   JJSIZE, JJSTART, JSIZEA, JSIZEB, JSTARTA,
     $                   JSTARTB, KK, LDA, LDB, LINFO, MM, MYID, MYPCOL,
     $                   MYPROW, NCOLA, NCOLB, NDIMA, NDIMB, NFREE, NN,
     $                   NPCOL, NPROC, NPROW, NROWA, NROWB
      DOUBLE PRECISION   ALPHA, BETA, ONE, ZERO
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA( DLEN_ ), DESCB( DLEN_ )
*     ..
*     .. External Functions ..
      INTEGER            NUMROC
      EXTERNAL           NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           ASSERT, BLACS_GRIDINFO, BLACS_PINFO, CHCALSIZE,
     $                   DESCINIT, DLAREAD, DLAWRITE, PDGEMM, PDPOTRF,
     $                   PDSYRK, PDTRSM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN
*     ..
*     .. Executable Statements ..
      ONE = DBLE( 1 )
      ZERO = DBLE( 0 )
      CALL BLACS_PINFO( MYID, NPROC )
      CALL BLACS_GRIDINFO( DESCC( CTXT_ ), NPROW, NPCOL, MYPROW,
     $                     MYPCOL )
      IFREE = 1
      NFREE = LWORK - IFREE + 1
*
*       Check parameters
*
      ISVALID = ( 1.LE.IC ) .AND. ( 1.LE.JC ) .AND. ( M.GE.1 ) .AND.
     $          ( IC+M-1.LE.DESCC( M_ ) ) .AND.
     $          ( JC+M-1.LE.DESCC( N_ ) )
      IF( .NOT.ISVALID ) THEN
         WRITE( *, FMT = 9999 )MYID, M, IC, JC, DESCC( M_ ),
     $      DESCC( N_ )
 9999    FORMAT( ' myid ', I4, ' m,ic,jc ', 3( 1X, I6 ), ' M,N ',
     $         2( 1X, I6 ) )
      ENDIF
      CALL ASSERT( ISVALID, '** PFCHFACT: error in parameters, myid ',
     $             MYID )
      IP_SAVE = IFREE
      NFREE = LWORK - IFREE + 1
      JSTARTB = 1
   10 CONTINUE
      IF( .NOT.( JSTARTB.LE.M ) )
     $   GOTO 90
      ISTARTB = JSTARTB
      IENDB = M
      ISIZEB = IENDB - ISTARTB + 1
*
*        Allocate storage for variable width
*        column panels A(*) and B(*)
*
      ISTARTA = ISTARTB
      IENDA = M
      ISIZEA = IENDA - ISTARTA + 1
      CALL CHCALSIZE( DESCC, M, NNB, JSTARTB, NFREE, GNCOLA, GNCOLB )
      NROWB = NUMROC( ISIZEB, DESCC( MB_ ), MYPROW, DESCC( RSRC_ ),
     $        NPROW )
      NROWA = NUMROC( ISIZEA, DESCC( MB_ ), MYPROW, DESCC( RSRC_ ),
     $        NPROW )
      NCOLA = NUMROC( GNCOLA, DESCC( NB_ ), MYPCOL, DESCC( CSRC_ ),
     $        NPCOL )
      NDIMA = NCOLA
      NCOLB = NUMROC( GNCOLB, DESCC( NB_ ), MYPCOL, DESCC( CSRC_ ),
     $        NPCOL )
      NDIMB = NCOLB
      LDA = NROWA
      LDB = NROWB
      ASIZE = LDA*NDIMA
      BSIZE = LDB*NDIMB
      IP_AMAT = IFREE
      IFREE = IFREE + ASIZE
      NFREE = LWORK - IFREE + 1
      CALL ASSERT( NFREE.GE.0, '** PFCHFACT: increase lwork by ',
     $             ABS( NFREE )+2 )
      IP_BMAT = IFREE
      IFREE = IFREE + BSIZE
      NFREE = LWORK - IFREE + 1
      CALL ASSERT( NFREE.GE.0, '** PFCHFACT: increase lwork by ',
     $             ABS( NFREE )+2 )
      INFO = 0
      CALL DESCINIT( DESCA, ISIZEB, GNCOLA, DESCC( MB_ ), DESCC( NB_ ),
     $               DESCC( RSRC_ ), DESCC( CSRC_ ), DESCC( CTXT_ ),
     $               MAX( 1, LDA ), INFO )
      CALL ASSERT( INFO.EQ.0,
     $             '** PFCHFACT: descinit for A returns info ', INFO )
      INFO = 0
      CALL DESCINIT( DESCB, ISIZEB, GNCOLB, DESCC( MB_ ), DESCC( NB_ ),
     $               DESCC( RSRC_ ), DESCC( CSRC_ ), DESCC( CTXT_ ),
     $               MAX( 1, LDB ), INFO )
      CALL ASSERT( INFO.EQ.0,
     $             '** PFCHFACT: descinit for B returns info ', INFO )
*
*         --------------------------
*         Ready to start computation.
*         --------------------------
*
      JENDB = MIN( M, JSTARTB+DESCB( N_ )-1 )
      JSIZEB = JENDB - JSTARTB + 1
*
*       Read in part of matrix into B
*
      ICPOS = ( IC-1 ) + ISTARTB
      JCPOS = ( JC-1 ) + JSTARTB
      INFO = 0
      CALL DLAREAD( IODEV, ISIZEB, JSIZEB, ICPOS, JCPOS,
     $              WORK( IP_BMAT ), 1, 1, DESCB, INFO )
      CALL ASSERT( INFO.EQ.0, '** PFCHFACT: laread of B returns info ',
     $             INFO )
*
*         Left-look algorithm with variable panel width.
*         Bring in previously factored entries and perform update.
*
      JSTARTA = 1
   20 CONTINUE
      IF( .NOT.( JSTARTA.LE.JSTARTB-1 ) )
     $   GOTO 40
      JENDA = MIN( JSTARTB-1, JSTARTA+DESCA( N_ )-1 )
      JSIZEA = JENDA - JSTARTA + 1
* read in factor
      IA = 1
      JA = 1
      ICPOS = ( IC-1 ) + ISTARTA
      JCPOS = ( JC-1 ) + JSTARTA
      INFO = 0
      CALL DLAREAD( IODEV, ISIZEA, JSIZEA, ICPOS, JCPOS,
     $              WORK( IP_AMAT ), IA, JA, DESCA, INFO )
      CALL ASSERT( INFO.EQ.0, '** PFCHFACT: laread of A returns info ',
     $             INFO )
*
*               Perform update.
*               B <- B - A*sub(A)^t
*
*
*               Use 2 calls to update separately the
*               triangular part and rectangular part of B.
*
      UPLO = 'L'
      TRANS = 'N'
      NN = JSIZEB
      KK = JSIZEA
      IA = 1
      JA = 1
      IB = 1
      JB = 1
      ALPHA = ( -ONE )
      BETA = ONE
*
*               Reminder:
*
*               C <- alpha*A*A^t + beta*C
*               where
*               C is nn by nn, A is nn by kk
*
      CALL PDSYRK( UPLO, TRANS, NN, KK, ALPHA, WORK( IP_AMAT ), IA, JA,
     $             DESCA, BETA, WORK( IP_BMAT ), IB, JB, DESCB )
      MM = ISIZEB - JSIZEB
      NN = JSIZEB
      KK = JSIZEA
      HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 ) .AND. ( KK.GE.1 )
      IF( HASWORK ) THEN
*
*                       C <- alpha op(A)*op(B) + beta*C,
*                       C is mm by nn
*
         ALPHA = ( -ONE )
         BETA = ONE
         IA = 1 + JSIZEB
         JA = 1
         IA2 = 1
         JA2 = JA
         IB = 1 + JSIZEB
         JB = 1
         CALL PDGEMM( 'N', 'T', MM, NN, KK, ALPHA, WORK( IP_AMAT ), IA,
     $                JA, DESCA, WORK( IP_AMAT ), IA2, JA2, DESCA, BETA,
     $                WORK( IP_BMAT ), IB, JB, DESCB )
      ENDIF
* end if (haswork)
   30 CONTINUE
      JSTARTA = JSTARTA + DESCA( N_ )
      GOTO 20
   40 CONTINUE
* end for jstartA
*
*       All updates are completed
*
      IF( USE_PTRSM ) THEN
*
*               Factor diagonal block
*
         UPLO = 'L'
         IB = 1
         JB = 1
         LINFO = 0
         CALL PDPOTRF( UPLO, JSIZEB, WORK( IP_BMAT ), IB, JB, DESCB,
     $                 LINFO )
         IF( INFO.GE.1 ) THEN
*
*               Part of matrix is not positive definite.
*
            INFO = ( JSTARTB-1 ) + LINFO
            RETURN
         ENDIF
         CALL ASSERT( INFO.EQ.0, '** PFCHFACT: PPOTRF returns info ',
     $                INFO )
*
*               Modify the rest of column.
*
         ISTARTB = JSTARTB + JSIZEB
         IENDB = M
         ISIZEB = IENDB - ISTARTB + 1
         HASWORK = ( ISIZEB.GE.1 )
         IF( HASWORK ) THEN
            SIDE = 'R'
            UPLO = 'L'
            TRANS = 'T'
            DIAG = 'N'
            ALPHA = ONE
            IA = 1
            JA = 1
            IB = JSIZEB + 1
            JB = 1
            CALL PDTRSM( SIDE, UPLO, TRANS, DIAG, ISIZEB, JSIZEB, ALPHA,
     $                   WORK( IP_BMAT ), IA, JA, DESCB,
     $                   WORK( IP_BMAT ), IB, JB, DESCB )
         ENDIF
* end if haswork
      ELSE
*
*               Perform cholesky factorization and
*               update together by unrolling the computation.
*               Seems to give better performance than
*               one call to PTRSM.
*
*               This may be considered a 'Right-looking' algorithm
*               to factor the in-core panel.
*
*
         ISTARTB = JSTARTB
         IENDB = M
         ISIZEB = IENDB - ISTARTB + 1
         JJINC = DESCB( NB_ )
         JJSTART = 1
   50    CONTINUE
         IF( .NOT.( JJSTART.LE.JSIZEB ) )
     $      GOTO 70
         JJEND = MIN( JSIZEB, JJSTART+JJINC-1 )
         JJSIZE = JJEND - JJSTART + 1
*
*               Cholesky factorization of diagonal block
*
         UPLO = 'L'
         IB = JJSTART
         JB = JJSTART
         NN = JJSIZE
         LINFO = 0
         CALL PDPOTRF( UPLO, NN, WORK( IP_BMAT ), IB, JB, DESCB, LINFO )
         IF( INFO.GE.1 ) THEN
*
*                       Part of matrix is not positive definite.
*
            INFO = ( JSTARTB-1 ) + ( JJSTART-1 ) + LINFO
            RETURN
         ENDIF
         CALL ASSERT( INFO.EQ.0, '** PFCHFACT: PPOTRF return info ',
     $                INFO )
*
*               Update current column
*
         SIDE = 'R'
         UPLO = 'L'
         TRANS = 'T'
         DIAG = 'N'
         ALPHA = ONE
         IA = JJSTART
         JA = JJSTART
         IB = JJEND + 1
         JB = JJSTART
         NN = ISIZEB - ( JJEND+1 ) + 1
         KK = JJSIZE
         HASWORK = ( NN.GE.1 ) .AND. ( KK.GE.1 )
         IF( HASWORK ) THEN
            CALL PDTRSM( SIDE, UPLO, TRANS, DIAG, NN, KK, ALPHA,
     $                   WORK( IP_BMAT ), IA, JA, DESCB,
     $                   WORK( IP_BMAT ), IB, JB, DESCB )
         ENDIF
*
*               Update  diagonal block to right
*
         UPLO = 'L'
         TRANS = 'N'
         ALPHA = ( -ONE )
         BETA = ONE
         IA = JJEND + 1
         JA = JJSTART
         IB = JJEND + 1
         JB = JJEND + 1
         NN = JSIZEB - ( JJEND+1 ) + 1
         KK = JJSIZE
         HASWORK = ( NN.GE.1 ) .AND. ( KK.GE.1 )
         IF( HASWORK ) THEN
            CALL PDSYRK( UPLO, TRANS, NN, KK, ALPHA, WORK( IP_BMAT ),
     $                   IA, JA, DESCB, BETA, WORK( IP_BMAT ), IB, JB,
     $                   DESCB )
         ENDIF
*
*               Update the rest of the right panel
*
         MM = ISIZEB - JSIZEB
         NN = JSIZEB - ( JJEND+1 ) + 1
         KK = JJSIZE
         ALPHA = ( -ONE )
         BETA = ONE
         IA = ISIZEB - MM + 1
         JA = JJSTART
         IA2 = JJEND + 1
         JA2 = JA
         IB = IA
         JB = ( JJEND+1 )
         HASWORK = ( ( MM.GE.1 ) .AND. ( NN.GE.1 ) .AND. ( KK.GE.1 ) )
         IF( HASWORK ) THEN
            CALL PDGEMM( 'N', 'T', MM, NN, KK, ALPHA, WORK( IP_BMAT ),
     $                   IA, JA, DESCB, WORK( IP_BMAT ), IA2, JA2,
     $                   DESCB, BETA, WORK( IP_BMAT ), IB, JB, DESCB )
         ENDIF
   60    CONTINUE
         JJSTART = JJEND + 1
         GOTO 50
   70    CONTINUE
* end for jjstart
      ENDIF
* end if (usePTRSM)
*
*       Write out the recently computed factors
*       back out to disk
*
      ISTARTB = JSTARTB
      IENDB = M
      ISIZEB = IENDB - ISTARTB + 1
      ICPOS = ( IC-1 ) + ISTARTB
      JCPOS = ( JC-1 ) + JSTARTB
      INFO = 0
      CALL DLAWRITE( IODEV, ISIZEB, JSIZEB, ICPOS, JCPOS,
     $               WORK( IP_BMAT ), 1, 1, DESCB, INFO )
      CALL ASSERT( INFO.EQ.0, '** PFCHFACT: lawrite of B returns info ',
     $             INFO )
*
*       Deallocate storage
*
      IFREE = IP_SAVE
      NFREE = LWORK - IFREE + 1
   80 CONTINUE
      JSTARTB = JENDB + 1
      GOTO 10
   90 CONTINUE
* end for
      RETURN
      END
