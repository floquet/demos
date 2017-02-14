      SUBROUTINE PFDUCHFACT( IODEV, M, IC, JC, WORK, LWORK, INFO )
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
* This routine uses the upper triangular part of matrix.
*
*
*
*     .. Parameters ..
      INTEGER            DLEN_
      PARAMETER          ( DLEN_ = 9 )
      LOGICAL            USE_PTRSM
      PARAMETER          ( USE_PTRSM = .false. )
*     ..
*     .. Scalar Arguments ..
      INTEGER            IC, INFO, IODEV, JC, LWORK, M
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   WORK( * )
*     ..
*     .. Local Scalars ..
      LOGICAL            HASWORK, LWORKQUERY
      INTEGER            ANEED, BNEED, CSRCIO, IB, ICOLA, ICOLB,
     $                   ICONTXT, IEND, IENDA, IENDB, IFREE, IIC, INCA,
     $                   INCB, INEED, IP_AMAT, IP_BMAT, IROWA, IROWB,
     $                   ISIZE, ISIZEA, ISIZEB, ISTART, ISTARTA,
     $                   ISTARTB, JB, JENDA, JENDB, JJC, JSIZEA, JSIZEB,
     $                   JSTARTA, JSTARTB, KK, LDA, LDB, LINFO, LOCP,
     $                   LOCPA, LOCPB, LOCQA, LOCQB, LWORK0, LWORK1, MB,
     $                   MM, MMB, MYPCOL, MYPROW, N, NB, NCOLA, NCOLB,
     $                   NFREE, NN, NNB, NPCOL, NPROW, NROWA, NROWB, P0,
     $                   Q0, RSRCIO
      DOUBLE PRECISION   ALPHA, BETA, ONE, ZERO
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA( DLEN_ ), DESCB( DLEN_ )
*     ..
*     .. External Functions ..
      INTEGER            INDXG2P, NUMROC
      EXTERNAL           INDXG2P, NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           ASSERT, BLACS_GRIDINFO, DESCINIT, DLAREAD,
     $                   DLAWRITE, LAIO_INFO, PDGEMM, PDPOTRF, PDSYRK,
     $                   PDTRSM, PFMAXSIZE, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
      N = M
      ONE = DBLE( 1 )
      ZERO = DBLE( 0 )
      CALL LAIO_INFO( IODEV, MM, NN, MMB, NNB, MB, NB, CSRCIO, RSRCIO,
     $                ICONTXT )
      CALL BLACS_GRIDINFO( ICONTXT, NPROW, NPCOL, MYPROW, MYPCOL )
*
* Determine storage requirement.
*
*
*   Storage for two mmb x n panels.
*
      P0 = MYPROW
      Q0 = MYPCOL
      LOCPB = NUMROC( MMB, MB, MYPROW, P0, NPROW )
      LOCQB = NUMROC( N, NB, MYPCOL, Q0, NPCOL )
      LWORK0 = MAX( LOCPB, 1 )*MAX( LOCQB, 1 )*2
*
*      Storage to hold entire matrix.
*
      LOCP = NUMROC( M, MB, MYPROW, P0, NPROW )
      LWORK1 = MAX( LOCP, 1 )*MAX( LOCQB, 1 )
      INEED = MIN( LWORK1, LWORK0 )
      LWORKQUERY = ( LWORK.EQ.-1 )
      IF( LWORK.LT.INEED ) THEN
         INFO = -6
         WORK( 1 ) = INEED
         IF( .NOT.LWORKQUERY ) THEN
            CALL PXERBLA( ICONTXT, 'PFDUCHFACT', -INFO )
         ENDIF
         RETURN
      ENDIF
*
*    Align B to perform local I/O operations.
*
      IROWB = INDXG2P( IC, MB, MYPROW, RSRCIO, NPROW )
      ICOLB = INDXG2P( JC, NB, MYPCOL, CSRCIO, NPCOL )
      P0 = IROWB
      Q0 = ICOLB
      NCOLB = N
      NROWB = -1
      CALL PFMAXSIZE( 'R', LWORK, NROWB, NCOLB, MB, NB, P0, Q0, ICONTXT,
     $                INFO )
      CALL ASSERT( INFO.EQ.0, '** PFxUCHFACT:  pfmaxsize returns ',
     $             INFO )
      NROWB = MIN( NROWB, M )
      IF( NROWB.NE.M ) THEN
         IF( NROWB.GE.MMB ) THEN
            NROWB = NROWB - MOD( NROWB, MMB )
         ENDIF
         IF( NROWB.GE.MB ) THEN
            NROWB = NROWB - MOD( NROWB, MB )
         ENDIF
      ENDIF
*
*       Allocate storage.
*
      IFREE = 1
      IP_BMAT = IFREE
      LOCPB = NUMROC( NROWB, MB, MYPROW, P0, NPROW )
      LDB = MAX( 1, LOCPB )
      CALL DESCINIT( DESCB, NROWB, NCOLB, MB, NB, P0, Q0, ICONTXT, LDB,
     $               INFO )
      CALL ASSERT( INFO.EQ.0, '** PFxUCHFACT: descinit returns info ',
     $             INFO )
      ISTARTB = 1
      INCB = NROWB
      IENDB = MIN( M, ISTARTB+INCB-1 )
      ISIZEB = IENDB - ISTARTB + 1
      JSTARTB = ISTARTB
      JENDB = M
      JSIZEB = JENDB - JSTARTB + 1
*
*       Read in
*
*       B <- C(istartB:iendB, jstartB:jendB )
*
      INFO = 0
      IIC = ( IC-1 ) + ISTARTB
      JJC = ( JC-1 ) + JSTARTB
      CALL DLAREAD( IODEV, ISIZEB, JSIZEB, IIC, JJC, WORK( IP_BMAT ), 1,
     $              1, DESCB, INFO )
      CALL ASSERT( INFO.EQ.0, '** PFxUCHFACT: LAREAD returns info = ',
     $             INFO )
      IF( USE_PTRSM ) THEN
*
*       Factor diagonal block
*
*        B(1:isizeB,1:isizeB) <- chol(B(1:isizeB,1:isizeB))
*
         LINFO = 0
         CALL PDPOTRF( 'Upper', ISIZEB, WORK( IP_BMAT ), 1, 1, DESCB,
     $                 LINFO )
         IF( LINFO.GE.1 ) THEN
*
*         Part of matrix is not positive definite.
*
            WORK( 1 ) = INEED
            INFO = ( ISTARTB-1 ) + LINFO
            RETURN
         ENDIF
*
*       Update the rest of B.
*
*       B(1:isizeB, (isizeB+1):jsizeB ) <-
*           B(1:isizeB,1:isizeB)'\B(1:isizeB, (isizeB+1):jsizeB )
*
         MM = ISIZEB
         NN = JSIZEB - ( ISIZEB+1 ) + 1
         ALPHA = ONE
         IB = 1
         JB = ( ISIZEB+1 )
         HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
         IF( HASWORK ) THEN
            CALL PDTRSM( 'Leftside', 'UpperTriangular', 'T',
     $                   'NonunitDiagonal', MM, NN, ALPHA,
     $                   WORK( IP_BMAT ), 1, 1, DESCB, WORK( IP_BMAT ),
     $                   IB, JB, DESCB )
         ENDIF
      ELSE
*
*        Perform cholesky factorization and
*        update together by unrolling the computation.
*        Seems to give better performance than
*        one call to PTRSM.
*
*        This may be considered a 'Right-looking' algorithm
*        to factor the in-core panel.
*
*
*
*        factor B(1:isizeB, 1:jsizeB ).
*
         DO 10 ISTART = 1, ISIZEB, MB
            IEND = MIN( ISIZEB, ISTART+MB-1 )
            ISIZE = IEND - ISTART + 1
*
*          Factor diagonal block
*
*          B(istart:iend,istart:iend).
*
            LINFO = 0
            CALL PDPOTRF( 'Upper', ISIZE, WORK( IP_BMAT ), ISTART,
     $                    ISTART, DESCB, LINFO )
            IF( LINFO.GE.1 ) THEN
*
*                Matrix may not be positive definite.
*
               INFO = ( ISTARTB-1 ) + ISTART
               WORK( 1 ) = INEED
               RETURN
            ENDIF
*
*       Update rest of row/column.
*
*       B(istart:iend, (iend+1):jsizeB ) =
*               B(istart:iend,istart:iend)' \
*                       B(istart:iend, (iend+1):jsizeB )
*
            MM = ISIZE
            NN = JSIZEB - ( IEND+1 ) + 1
            ALPHA = ONE
            HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
            IF( HASWORK ) THEN
               CALL PDTRSM( 'Leftside', 'UpperTriangular', 'T',
     $                      'NonunitDiagonal', MM, NN, ALPHA,
     $                      WORK( IP_BMAT ), ISTART, ISTART, DESCB,
     $                      WORK( IP_BMAT ), ISTART, ( IEND+1 ), DESCB )
            ENDIF
*
*  Symmetric update of diagonal block.
*
*  B( (iend+1):isizeB, (iend+1):isizeB ) =
*       B( (iend+1):isizeB, (iend+1):isizeB ) -
*          B( istart:iend, (iend+1):isizeB )' *
*               B(istart:iend, (iend+1):isizeB)
*
            NN = ISIZEB - ( IEND+1 ) + 1
            KK = ISIZE
            ALPHA = -ONE
            BETA = ONE
            HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
            IF( HASWORK ) THEN
               CALL PDSYRK( 'Upper', 'T', NN, KK, ALPHA,
     $                      WORK( IP_BMAT ), ISTART, ( IEND+1 ), DESCB,
     $                      BETA, WORK( IP_BMAT ), ( IEND+1 ),
     $                      ( IEND+1 ), DESCB )
            ENDIF
*
* Update the rest
*
* B( (iend+1):isizeB, (isizeB+1):jsizeB ) <-
*       B( istart:iend, (iend+1):isizeB )' *
*                       B(istart:iend,  (isizeB+1):jsizeB )
*
            MM = ISIZEB - ( IEND+1 ) + 1
            NN = JSIZEB - ( ISIZEB+1 ) + 1
            KK = ISIZE
            ALPHA = -ONE
            BETA = ONE
            HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 ) .AND. ( KK.GE.1 )
            IF( HASWORK ) THEN
               CALL PDGEMM( 'T', 'NoTranspose', MM, NN, KK, ALPHA,
     $                      WORK( IP_BMAT ), ISTART, ( IEND+1 ), DESCB,
     $                      WORK( IP_BMAT ), ISTART, ( ISIZEB+1 ),
     $                      DESCB, BETA, WORK( IP_BMAT ), ( IEND+1 ),
     $                      ( ISIZEB+1 ), DESCB )
            ENDIF
   10    CONTINUE
   20    CONTINUE
* end do istart
      ENDIF
      IIC = ( IC-1 ) + ISTARTB
      JJC = ( JC-1 ) + JSTARTB
      CALL DLAWRITE( IODEV, ISIZEB, JSIZEB, IIC, JJC, WORK( IP_BMAT ),
     $               1, 1, DESCB, INFO )
      CALL ASSERT( INFO.EQ.0, '** PFxUCHFACT: LAWRITE returns info = ',
     $             INFO )
*
*  Factor the rest of matrix.
*
      ISTARTB = NROWB + 1
   30 CONTINUE
      IF( ISTARTB.LE.M ) THEN
         IFREE = 1
         JSTARTB = ISTARTB
         JENDB = N
         JSIZEB = JENDB - JSTARTB + 1
         IIC = ( IC-1 ) + ISTARTB
         JJC = ( JC-1 ) + JSTARTB
         IROWB = INDXG2P( IIC, MB, MYPROW, RSRCIO, NPROW )
         ICOLB = INDXG2P( JJC, NB, MYPCOL, CSRCIO, NPCOL )
*
*               Partition storage for variable width panel B.
*               panelA is mmb x jendB
*
*
*           Note the use of p0 = myprow, q0 = mypcol
*           to over estimate storage requirements.
*
         P0 = MYPROW
         Q0 = MYPCOL
         NROWA = MIN( MMB, ISTARTB-1 )
         NCOLA = MIN( N, MAX( 1, JSIZEB ) )
         LOCPA = NUMROC( NROWA, MB, MYPROW, P0, NPROW )
         LOCQA = NUMROC( NCOLA, NB, MYPCOL, Q0, NPCOL )
         ANEED = LOCPA*LOCQA
         IP_AMAT = IFREE
         IFREE = IFREE + ANEED
         NFREE = LWORK - ANEED
         CALL ASSERT( NFREE.GE.1, '** PFxUCHFACT: nfree <= 0 ', NFREE )
*
*           The rest of storage goes to B.
*
         P0 = IROWB
         Q0 = ICOLB
         NCOLB = JSIZEB
         NROWB = -1
         CALL PFMAXSIZE( 'Row', NFREE, NROWB, NCOLB, MB, NB, P0, Q0,
     $                   ICONTXT, INFO )
         NROWB = MIN( NROWB, ( M-ISTARTB+1 ) )
         IF( NROWB.NE.( M-ISTARTB+1 ) ) THEN
            IF( NROWB.GE.MMB ) THEN
               NROWB = NROWB - MOD( NROWB, MMB )
            ENDIF
            IF( NROWB.GE.MB ) THEN
               NROWB = NROWB - MOD( NROWB, MB )
            ENDIF
         ENDIF
*
*            Double check storage.
*
         LOCPB = NUMROC( NROWB, MB, MYPROW, P0, NPROW )
         LOCQB = NUMROC( NCOLB, NB, MYPCOL, Q0, NPCOL )
         BNEED = LOCPB*LOCQB
         CALL ASSERT( BNEED.LE.NFREE,
     $                '** PFxUCHFACT: need more storage, Bneed = ',
     $                BNEED )
         IP_BMAT = IFREE
         IFREE = IFREE + BNEED
         LDB = MAX( 1, LOCPB )
         CALL DESCINIT( DESCB, NROWB, NCOLB, MB, NB, P0, Q0, ICONTXT,
     $                  LDB, INFO )
         CALL ASSERT( INFO.EQ.0,
     $                '** PFxUCHFACT: descinit returns info =', INFO )
*
*         Read in
*         B <- C(istartB:iendB, jstartB:jendB )
*
         INCB = NROWB
         IENDB = MIN( M, ISTARTB+INCB-1 )
         ISIZEB = IENDB - ISTARTB + 1
         INFO = 0
         IIC = ( IC-1 ) + ISTARTB
         JJC = ( JC-1 ) + JSTARTB
         CALL DLAREAD( IODEV, ISIZEB, JSIZEB, IIC, JJC, WORK( IP_BMAT ),
     $                 1, 1, DESCB, INFO )
         CALL ASSERT( INFO.EQ.0, '** PFxUCHFACT: LAREAD returns info = '
     $                , INFO )
         ISTARTA = 1
         INCA = MMB
   40    CONTINUE
         IF( ISTARTA.LE.( ISTARTB-1 ) ) THEN
            IENDA = MIN( ( ISTARTB-1 ), ISTARTA+INCA-1 )
            ISIZEA = IENDA - ISTARTA + 1
            JSTARTA = JSTARTB
            JENDA = JENDB
            JSIZEA = JENDA - JSTARTA + 1
            IIC = ( IC-1 ) + ISTARTA
            JJC = ( JC-1 ) + JSTARTB
*
*               Align panelA to perform local I/O operations.
*
            IROWA = INDXG2P( IIC, MB, MYPROW, RSRCIO, NPROW )
            ICOLA = INDXG2P( JJC, NB, MYPCOL, CSRCIO, NPCOL )
            P0 = IROWA
            Q0 = ICOLA
            NROWA = ISIZEA
            NCOLA = JSIZEB
            INFO = 0
            LOCPA = NUMROC( NROWA, MB, MYPROW, P0, NPROW )
            LDA = MAX( 1, LOCPA )
            CALL DESCINIT( DESCA, NROWA, NCOLA, MB, NB, P0, Q0, ICONTXT,
     $                     LDA, INFO )
            CALL ASSERT( INFO.EQ.0,
     $                   '** PFxUCHFACT: descinit returns info =',
     $                   INFO )
*
*                 Read in
*                 A <- C(istartA:iendA,  jstartB:jendB );
*
*
            INFO = 0
            IIC = ( IC-1 ) + ISTARTA
            JJC = ( JC-1 ) + JSTARTB
            CALL DLAREAD( IODEV, ISIZEA, JSIZEB, IIC, JJC,
     $                    WORK( IP_AMAT ), 1, 1, DESCA, INFO )
            CALL ASSERT( INFO.EQ.0,
     $                   '** PFxUCHFACT: LAREAD returns info ', INFO )
*
*               Symmetric update to diagonal block.
*
*               B(1:isizeB,1:isizeB) <- B(1:isizeB,1:isizeB) -
*                       A(1:isizeA,1:isizeB)'*
*                               A(1:isizeA,1:isizeB)
*
            ALPHA = -ONE
            BETA = ONE
            NN = ISIZEB
            KK = ISIZEA
            HASWORK = ( NN.GE.1 ) .AND. ( KK.GE.1 )
            IF( HASWORK ) THEN
               CALL PDSYRK( 'Upper', 'T', NN, KK, ALPHA,
     $                      WORK( IP_AMAT ), 1, 1, DESCA, BETA,
     $                      WORK( IP_BMAT ), 1, 1, DESCB )
            ENDIF
*
*               Update remaining part of B.
*
*               B(1:isizeB, (isizeB+1):jsizeB ) =
*                       B(1:isizeB, (isizeB+1):jsizeB ) -
*                       A(1:isizeA,1:isizeB)'*
*                               A(1:isizeA, (isizeB+1):jsizeB)
*
*
            MM = ISIZEB
            NN = JSIZEB - ( ISIZEB+1 ) + 1
            KK = ISIZEA
            ALPHA = -ONE
            BETA = ONE
            HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 ) .AND. ( KK.GE.1 )
            IF( HASWORK ) THEN
               CALL PDGEMM( 'T', 'NoTranspose', MM, NN, KK, ALPHA,
     $                      WORK( IP_AMAT ), 1, 1, DESCA,
     $                      WORK( IP_AMAT ), 1, ( ISIZEB+1 ), DESCA,
     $                      BETA, WORK( IP_BMAT ), 1, ( ISIZEB+1 ),
     $                      DESCB )
            ENDIF
            ISTARTA = IENDA + 1
            GOTO 40
         ENDIF
   50    CONTINUE
* end while (istartA <= (istartB-1) )
*
*        All updates are done.
*
         IF( USE_PTRSM ) THEN
*
*        Factor diagonal block.
*        B(1:isizeB,1:isizeB) <- chol( B(1:isizeB,1:isizeB) );
*
            LINFO = 0
            CALL PDPOTRF( 'Upper', ISIZEB, WORK( IP_BMAT ), 1, 1, DESCB,
     $                    LINFO )
            IF( LINFO.GE.1 ) THEN
*
*               Part of matrix is not positive definite.
*
               WORK( 1 ) = INEED
               INFO = ( ISTARTB-1 ) + LINFO
               RETURN
            ENDIF
*
*         Modify the rest of B.
*
*         B(1:isizeB, (isizeB+1):jsizeB ) =
*               B(1:isizeB,1:isizeB)' \
*                       B(1:isizeB, (isizeB+1):jsizeB )
*
*
            MM = ISIZEB
            NN = JSIZEB - ( ISIZEB+1 ) + 1
            ALPHA = ONE
            HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
            IF( HASWORK ) THEN
               CALL PDTRSM( 'Leftside', 'UpperPart', 'T',
     $                      'NonunitDiagonal', MM, NN, ALPHA,
     $                      WORK( IP_BMAT ), 1, 1, DESCB,
     $                      WORK( IP_BMAT ), 1, ( ISIZEB+1 ), DESCB )
            ENDIF
         ELSE
*
*        Perform cholesky factorization and
*        update together by unrolling the computation.
*        Seems to give better performance than
*        one call to PTRSM.
*
*        This may be considered a 'Right-looking' algorithm
*        to factor the in-core panel.
*
*
*
*        factor B(1:isizeB, 1:jsizeB ).
*
            DO 60 ISTART = 1, ISIZEB, MB
               IEND = MIN( ISIZEB, ISTART+MB-1 )
               ISIZE = IEND - ISTART + 1
*
*          Factor diagonal block B(istart:iend,istart:iend).
*
               LINFO = 0
               CALL PDPOTRF( 'Upper', ISIZE, WORK( IP_BMAT ), ISTART,
     $                       ISTART, DESCB, LINFO )
               IF( LINFO.GE.1 ) THEN
*
*                Matrix may not be positive definite.
*
                  INFO = ( ISTARTB-1 ) + ISTART
                  WORK( 1 ) = INEED
                  RETURN
               ENDIF
*
*       Update rest of column.
*       B(istart:iend, (iend+1):jsizeB ) =
*               B(istart:iend,istart:iend)' \
*                       B(istart:iend, (iend+1):jsizeB )
*
               MM = ISIZE
               NN = JSIZEB - ( IEND+1 ) + 1
               ALPHA = ONE
               HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
               IF( HASWORK ) THEN
                  CALL PDTRSM( 'Leftside', 'UpperTriangular', 'T',
     $                         'NonunitDiagonal', MM, NN, ALPHA,
     $                         WORK( IP_BMAT ), ISTART, ISTART, DESCB,
     $                         WORK( IP_BMAT ), ISTART, ( IEND+1 ),
     $                         DESCB )
               ENDIF
*
*  Symmetric update of diagonal block.
*
*  B( (iend+1):isizeB, (iend+1):isizeB ) =
*       B( (iend+1):isizeB, (iend+1):isizeB ) -
*          B( istart:iend, (iend+1):isizeB )' *
*               B(istart:iend, (iend+1):isizeB)
*
               NN = ISIZEB - ( IEND+1 ) + 1
               KK = ISIZE
               ALPHA = -ONE
               BETA = ONE
               HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
               IF( HASWORK ) THEN
                  CALL PDSYRK( 'Upper', 'T', NN, KK, ALPHA,
     $                         WORK( IP_BMAT ), ISTART, ( IEND+1 ),
     $                         DESCB, BETA, WORK( IP_BMAT ), ( IEND+1 ),
     $                         ( IEND+1 ), DESCB )
               ENDIF
*
* Update the rest
*
* B( (iend+1):isizeB, (isizeB+1):jsizeB ) <-
*       B( istart:iend, (iend+1):isizeB )' *
*                       B(istart:iend,  (isizeB+1):jsizeB )
*
               MM = ISIZEB - ( IEND+1 ) + 1
               NN = JSIZEB - ( ISIZEB+1 ) + 1
               KK = ISIZE
               ALPHA = -ONE
               BETA = ONE
               HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 ) .AND. ( KK.GE.1 )
               IF( HASWORK ) THEN
                  CALL PDGEMM( 'T', 'NoTranspose', MM, NN, KK, ALPHA,
     $                         WORK( IP_BMAT ), ISTART, ( IEND+1 ),
     $                         DESCB, WORK( IP_BMAT ), ISTART,
     $                         ( ISIZEB+1 ), DESCB, BETA,
     $                         WORK( IP_BMAT ), ( IEND+1 ),
     $                         ( ISIZEB+1 ), DESCB )
               ENDIF
   60       CONTINUE
   70       CONTINUE
* end do istart
         ENDIF
*
*        Write out B.
*
*         C(istartB:iendB, jstartB:jendB) <- B
*
         INFO = 0
         IIC = ( IC-1 ) + ISTARTB
         JJC = ( JC-1 ) + JSTARTB
         CALL DLAWRITE( IODEV, ISIZEB, JSIZEB, IIC, JJC,
     $                  WORK( IP_BMAT ), 1, 1, DESCB, INFO )
         CALL ASSERT( INFO.EQ.0,
     $                '** PFxUCHFACT: LAWRITE returns info = ', INFO )
         ISTARTB = IENDB + 1
         GOTO 30
      ENDIF
   80 CONTINUE
* end while (istartB <= m)
      WORK( 1 ) = INEED
      INFO = 0
      RETURN
      END
