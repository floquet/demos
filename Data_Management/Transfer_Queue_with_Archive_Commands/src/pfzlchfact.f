      SUBROUTINE PFZLCHFACT( IODEV, M, IC, JC, WORK, LWORK, INFO )
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
* This routine uses the lower triangular part of matrix.
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
      COMPLEX*16         WORK( * )
*     ..
*     .. Local Scalars ..
      LOGICAL            HASWORK, LWORKQUERY
      INTEGER            ANEED, BNEED, CSRCIO, IB, ICOLA, ICOLB,
     $                   ICONTXT, IEND, IENDA, IENDB, IFREE, IIC, INC,
     $                   INCA, INEED, IP_AMAT, IP_BMAT, IROWA, IROWB,
     $                   ISIZE, ISIZEA, ISIZEB, ISTART, ISTARTA,
     $                   ISTARTB, JB, JENDA, JENDB, JJC, JSIZEA, JSIZEB,
     $                   JSTARTA, JSTARTB, KK, LDA, LDB, LINFO, LOCPA,
     $                   LOCPB, LOCQ, LOCQA, LOCQB, LWORK0, LWORK1, MB,
     $                   MM, MMB, MYPCOL, MYPROW, N, NB, NCOLA, NCOLB,
     $                   NFREE, NN, NNB, NPCOL, NPROW, NROWA, NROWB, P0,
     $                   Q0, RSRCIO
      COMPLEX*16         ALPHA, BETA, ONE, ZERO
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
     $                   PFMAXSIZE, PXERBLA, PZGEMM, PZHERK, PZPOTRF,
     $                   PZTRSM, ZLAREAD, ZLAWRITE
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, INT, MAX, MIN
*     ..
*     .. Executable Statements ..
      N = M
      ONE = DCMPLX( DBLE( 1 ) )
      ZERO = DCMPLX( DBLE( 0 ) )
      CALL LAIO_INFO( IODEV, MM, NN, MMB, NNB, MB, NB, CSRCIO, RSRCIO,
     $                ICONTXT )
      CALL BLACS_GRIDINFO( ICONTXT, NPROW, NPCOL, MYPROW, MYPCOL )
*
* Determine storage requirement.
*
*
*   Storage for two m x nnb panels.
*
      P0 = MYPROW
      Q0 = MYPCOL
      LOCPB = NUMROC( M, MB, MYPROW, P0, NPROW )
      LOCQB = NUMROC( NNB, NB, MYPCOL, Q0, NPCOL )
      LWORK0 = MAX( LOCPB, 1 )*MAX( LOCQB, 1 )*2
*
*      Storage to hold entire matrix.
*
      LOCQ = NUMROC( N, NB, MYPCOL, Q0, NPCOL )
      LWORK1 = MAX( LOCPB, 1 )*MAX( LOCQ, 1 )
      INEED = MIN( LWORK1, LWORK0 )
      LWORKQUERY = ( LWORK.EQ.-1 )
      IF( LWORK.LT.INEED ) THEN
         INFO = -6
         WORK( 1 ) = INEED
         IF( .NOT.LWORKQUERY ) THEN
            CALL PXERBLA( ICONTXT, 'PFZLCHFACT', -INFO )
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
      NROWB = M
      NCOLB = -1
      CALL PFMAXSIZE( 'C', LWORK, NROWB, NCOLB, MB, NB, P0, Q0, ICONTXT,
     $                INFO )
      CALL ASSERT( INFO.EQ.0, '** PFxLCHFACT:  pfmaxsize returns ',
     $             INFO )
      JSTARTB = 1
      JENDB = MIN( N, JSTARTB+NCOLB-1 )
      IF( JENDB.NE.N ) THEN
*
*       Try to align on block boundary
*
         JJC = ( JC-1 ) + JENDB
         JJC = MAX( JC, INT( JJC / NNB )*NNB )
         JENDB = JJC - ( JC-1 )
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
      CALL ASSERT( INFO.EQ.0, '** PFxLCHFACT: descinit returns info ',
     $             INFO )
      ISTARTB = 1
      IENDB = M
      ISIZEB = IENDB - ISTARTB + 1
      JSIZEB = JENDB - JSTARTB + 1
*
*       Read in
*
*       B <- C(istartB:iendB, jstartB:jendB )
*
      INFO = 0
      IIC = ( IC-1 ) + ISTARTB
      JJC = ( JC-1 ) + JSTARTB
      CALL ZLAREAD( IODEV, ISIZEB, JSIZEB, IIC, JJC, WORK( IP_BMAT ), 1,
     $              1, DESCB, INFO )
      CALL ASSERT( INFO.EQ.0, '** PFxLCHFACT: LAREAD returns info = ',
     $             INFO )
      IF( USE_PTRSM ) THEN
*
*       Factor diagonal block
*
*        B(1:jsizeB,1:jsizeB) <- chol(B(1:jsizeB,1:jsizeB))
*
         LINFO = 0
         CALL PZPOTRF( 'Lower', JSIZEB, WORK( IP_BMAT ), 1, 1, DESCB,
     $                 LINFO )
         IF( LINFO.GE.1 ) THEN
*
*         Part of matrix is not positive definite.
*
            WORK( 1 ) = INEED
            INFO = ( JSTARTB-1 ) + LINFO
            RETURN
         ENDIF
*
*       Update the rest of B.
*
*   B( (jsizeB+1):isizeB, 1:jsizeB ) = ...
*        B( (jsizeB+1):isizeB, 1:jsizeB )/ (  B(1:jsizeB,1:jsizeB)');
*
*
         MM = ISIZEB - ( JSIZEB+1 ) + 1
         NN = JSIZEB
         ALPHA = ONE
         IB = ( JSIZEB+1 )
         JB = 1
         HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
         IF( HASWORK ) THEN
            CALL PZTRSM( 'Rightside', 'LowerTriangular', 'C',
     $                   'NonunitDiagonal', MM, NN, ALPHA,
     $                   WORK( IP_BMAT ), 1, 1, DESCB, WORK( IP_BMAT ),
     $                   IB, JB, DESCB )
         ENDIF
      ELSE
*
*        Perform cholesky factorization and
*        update together by unrolling the computation.
*        Seems to give better performance than
*        one call to P_TRSM.
*
*        This may be considered a 'Right-looking' algorithm
*        to factor the in-core panel.
*
*
*
*        factor B(1:isizeB, 1:jsizeB ).
*
         INC = NB
         DO 10 ISTART = 1, JSIZEB, INC
            IEND = MIN( JSIZEB, ISTART+INC-1 )
            ISIZE = IEND - ISTART + 1
*
*          Factor diagonal block
*
*          B(istart:iend,istart:iend).
*
            LINFO = 0
            CALL PZPOTRF( 'Lower', ISIZE, WORK( IP_BMAT ), ISTART,
     $                    ISTART, DESCB, LINFO )
            IF( LINFO.GE.1 ) THEN
*
*                Matrix may not be positive definite.
*
               INFO = ( JSTARTB-1 ) + ( ISTART-1 ) + LINFO
               WORK( 1 ) = INEED
               RETURN
            ENDIF
*
*       Update rest of row/column.
*
*        B( (iend+1):isizeB, istart:iend) =...
*           B( (iend+1):isizeB, istart:iend) /
*                   B(istart:iend,istart:iend)'
*
            MM = ISIZEB - ( IEND+1 ) + 1
            NN = ISIZE
            ALPHA = ONE
            HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
            IF( HASWORK ) THEN
               CALL PZTRSM( 'Rightside', 'LowerTriangular', 'C',
     $                      'NonunitDiagonal', MM, NN, ALPHA,
     $                      WORK( IP_BMAT ), ISTART, ISTART, DESCB,
     $                      WORK( IP_BMAT ), ( IEND+1 ), ISTART, DESCB )
            ENDIF
*
*  Symmetric update of diagonal block.
*
*  B( (iend+1):jsizeB, (iend+1):jsizeB ) = ...
*       B( (iend+1):jsizeB, (iend+1):jsizeB ) -
*          B( (iend+1):jsizeB, istart:iend ) *
*               B((iend+1):jsizeB, istart:iend )'
*
            NN = JSIZEB - ( IEND+1 ) + 1
            KK = ISIZE
            ALPHA = -ONE
            BETA = ONE
            HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
            IF( HASWORK ) THEN
               CALL PZHERK( 'Lower', 'NoTranspose', NN, KK,
     $                      DBLE( ALPHA ), WORK( IP_BMAT ), ( IEND+1 ),
     $                      ISTART, DESCB, DBLE( BETA ),
     $                      WORK( IP_BMAT ), ( IEND+1 ), ( IEND+1 ),
     $                      DESCB )
            ENDIF
*
* Update the rest
*
* B( (jsizeB+1):isizeB, (iend+1):jsizeB ) <-
*    B( (jsizeB+1):isizeB, (iend+1):jsizeB )
*           B((jsizeB+1):isizeB,  istart:iend )*
*               B((iend+1):jsizeB, istart:iend )'
*
            MM = ISIZEB - ( JSIZEB+1 ) + 1
            NN = JSIZEB - ( IEND+1 ) + 1
            KK = ISIZE
            ALPHA = -ONE
            BETA = ONE
            HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 ) .AND. ( KK.GE.1 )
            IF( HASWORK ) THEN
               CALL PZGEMM( 'NoTranspose', 'C', MM, NN, KK, ALPHA,
     $                      WORK( IP_BMAT ), ( JSIZEB+1 ), ISTART,
     $                      DESCB, WORK( IP_BMAT ), ( IEND+1 ), ISTART,
     $                      DESCB, BETA, WORK( IP_BMAT ), ( JSIZEB+1 ),
     $                      ( IEND+1 ), DESCB )
            ENDIF
   10    CONTINUE
   20    CONTINUE
* end do istart
      ENDIF
*
*  Write out
*  C(istartB:iendB, jstartB:jendB) = B;
*
      IIC = ( IC-1 ) + ISTARTB
      JJC = ( JC-1 ) + JSTARTB
      CALL ZLAWRITE( IODEV, ISIZEB, JSIZEB, IIC, JJC, WORK( IP_BMAT ),
     $               1, 1, DESCB, INFO )
      CALL ASSERT( INFO.EQ.0, '** PFxLCHFACT: LAWRITE returns info = ',
     $             INFO )
*
*  Factor the rest of matrix.
*
      JSTARTB = JENDB + 1
   30 CONTINUE
      IF( JSTARTB.LE.N ) THEN
         ISTARTB = JSTARTB
         IENDB = M
         ISIZEB = IENDB - ISTARTB + 1
         ISTARTA = ISTARTB
         IENDA = IENDB
         ISIZEA = IENDA - ISTARTA + 1
         IFREE = 1
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
         JSIZEA = MAX( 1, MIN( ( JSTARTB-1 ), NNB ) )
         LOCPA = NUMROC( ISIZEA, MB, MYPROW, P0, NPROW )
         LOCQA = NUMROC( JSIZEA, NB, MYPCOL, Q0, NPCOL )
         ANEED = LOCPA*LOCQA
         IP_AMAT = IFREE
         IFREE = IFREE + ANEED
         NFREE = LWORK - ANEED
         CALL ASSERT( NFREE.GE.1, '** PFxLCHFACT: nfree <= 0 ', NFREE )
*
*           The rest of storage goes to B.
*
         P0 = IROWB
         Q0 = ICOLB
         NCOLB = JSIZEB
         NROWB = -1
         NROWB = ISIZEB
         NCOLB = -1
         CALL PFMAXSIZE( 'C', NFREE, NROWB, NCOLB, MB, NB, P0, Q0,
     $                   ICONTXT, INFO )
         JENDB = MIN( N, JSTARTB+NCOLB-1 )
         IF( JENDB.NE.N ) THEN
            JJC = ( JC-1 ) + JENDB
            JJC = MAX( ( JC-1 )+JSTARTB, INT( JJC / NNB )*NNB )
            JENDB = JJC - ( JC-1 )
         ENDIF
         JSIZEB = JENDB - JSTARTB + 1
*
*            Double check storage.
*
         LOCPB = NUMROC( ISIZEB, MB, MYPROW, P0, NPROW )
         LOCQB = NUMROC( JSIZEB, NB, MYPCOL, Q0, NPCOL )
         BNEED = LOCPB*LOCQB
         CALL ASSERT( BNEED.LE.NFREE,
     $                '** PFxLCHFACT: need more storage, Bneed = ',
     $                BNEED )
         IP_BMAT = IFREE
         IFREE = IFREE + BNEED
         LDB = MAX( 1, LOCPB )
         CALL DESCINIT( DESCB, NROWB, NCOLB, MB, NB, P0, Q0, ICONTXT,
     $                  LDB, INFO )
         CALL ASSERT( INFO.EQ.0,
     $                '** PFxLCHFACT: descinit returns info =', INFO )
*
*         Read in
*         B <- C(istartB:iendB, jstartB:jendB )
*
         INFO = 0
         IIC = ( IC-1 ) + ISTARTB
         JJC = ( JC-1 ) + JSTARTB
         CALL ZLAREAD( IODEV, ISIZEB, JSIZEB, IIC, JJC, WORK( IP_BMAT ),
     $                 1, 1, DESCB, INFO )
         CALL ASSERT( INFO.EQ.0, '** PFxLCHFACT: LAREAD returns info = '
     $                , INFO )
         JSTARTA = 1
         INCA = NNB
   40    CONTINUE
         IF( JSTARTA.LE.( JSTARTB-1 ) ) THEN
            JENDA = MIN( ( JSTARTB-1 ), JSTARTA+INCA-1 )
            IF( JENDA.NE.( JSTARTB-1 ) ) THEN
               JJC = ( JC-1 ) + JENDA
               JJC = MAX( ( JC-1 )+JSTARTA, INT( JJC / NNB )*NNB )
               JENDA = JJC - ( JC-1 )
            ENDIF
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
            NCOLA = JSIZEA
            INFO = 0
            LOCPA = NUMROC( NROWA, MB, MYPROW, P0, NPROW )
            LDA = MAX( 1, LOCPA )
            CALL DESCINIT( DESCA, NROWA, NCOLA, MB, NB, P0, Q0, ICONTXT,
     $                     LDA, INFO )
            CALL ASSERT( INFO.EQ.0,
     $                   '** PFxLCHFACT: descinit returns info =',
     $                   INFO )
*
*                 Read in
*               A = C(istartA:iendA, jstartA:jendA );
*
            INFO = 0
            IIC = ( IC-1 ) + ISTARTA
            JJC = ( JC-1 ) + JSTARTA
            CALL ZLAREAD( IODEV, ISIZEA, JSIZEA, IIC, JJC,
     $                    WORK( IP_AMAT ), 1, 1, DESCA, INFO )
            CALL ASSERT( INFO.EQ.0,
     $                   '** PFxLCHFACT: LAREAD returns info ', INFO )
*
*    Symmetric update to diagonal block.
*
*    B(1:jsizeB,1:jsizeB) = B(1:jsizeB,1:jsizeB) - ...
*           A(1:jsizeB,1:jsizeA)*A(1:jsizeB,1:jsizeA)';
*
*
            ALPHA = -ONE
            BETA = ONE
            NN = JSIZEB
            KK = JSIZEA
            HASWORK = ( NN.GE.1 ) .AND. ( KK.GE.1 )
            IF( HASWORK ) THEN
               CALL PZHERK( 'Lower', 'NoTranspose', NN, KK,
     $                      DBLE( ALPHA ), WORK( IP_AMAT ), 1, 1, DESCA,
     $                      DBLE( BETA ), WORK( IP_BMAT ), 1, 1, DESCB )
            ENDIF
*
*    Update remaining part of B.
*
*   B( (jsizeB+1):isizeB, 1:jsizeB) =
*           B( (jsizeB+1):isizeB, 1:jsizeB) - ...
*           A( (jsizeB+1):isizeA, 1:jsizeA) *
*           A(1:jsizeB,1:jsizeA)';
*
            MM = ISIZEB - ( JSIZEB+1 ) + 1
            NN = JSIZEB
            KK = JSIZEA
            ALPHA = -ONE
            BETA = ONE
            HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 ) .AND. ( KK.GE.1 )
            IF( HASWORK ) THEN
               CALL PZGEMM( 'NoTranspose', 'C', MM, NN, KK, ALPHA,
     $                      WORK( IP_AMAT ), ( JSIZEB+1 ), 1, DESCA,
     $                      WORK( IP_AMAT ), 1, 1, DESCA, BETA,
     $                      WORK( IP_BMAT ), ( JSIZEB+1 ), 1, DESCB )
            ENDIF
            JSTARTA = JENDA + 1
            GOTO 40
         ENDIF
   50    CONTINUE
* end while
*
*        All updates are done.
*
         IF( USE_PTRSM ) THEN
*
*     Factor diagonal block.
*
*    B(1:jsizeB,1:jsizeB) = chol( B(1:jsizeB, 1:jsizeB) );
*
*
            LINFO = 0
            CALL PZPOTRF( 'Lower', JSIZEB, WORK( IP_BMAT ), 1, 1, DESCB,
     $                    LINFO )
            IF( LINFO.GE.1 ) THEN
*
*               Part of matrix is not positive definite.
*
               WORK( 1 ) = INEED
               INFO = ( JSTARTB-1 ) + LINFO
               RETURN
            ENDIF
*
*      Modify the rest of B.
*
*    B( (jsizeB+1):isizeB,  1:jsizeB) = ...
*       B( (jsizeB+1):isizeB,  1:jsizeB) /
*               (B(1:jsizeB,1:jsizeB)');
*
            MM = ISIZEB - ( JSIZEB+1 ) + 1
            NN = JSIZEB
            ALPHA = ONE
            HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
            IF( HASWORK ) THEN
               CALL PZTRSM( 'Rightside', 'LowerTriangular', 'C',
     $                      'NonunitDiagonal', MM, NN, ALPHA,
     $                      WORK( IP_BMAT ), 1, 1, DESCB,
     $                      WORK( IP_BMAT ), ( JSIZEB+1 ), 1, DESCB )
            ENDIF
         ELSE
*
*        Perform cholesky factorization and
*        update together by unrolling the computation.
*        Seems to give better performance than
*        one call to P_TRSM.
*
*        This may be considered a 'Right-looking' algorithm
*        to factor the in-core panel.
*
*
*
*        factor B(1:isizeB, 1:jsizeB ).
*
            INC = NB
            DO 60 ISTART = 1, JSIZEB, INC
               IEND = MIN( JSIZEB, ISTART+INC-1 )
               ISIZE = IEND - ISTART + 1
*
*          Factor diagonal block
*
*          B(istart:iend,istart:iend).
*
               LINFO = 0
               CALL PZPOTRF( 'Lower', ISIZE, WORK( IP_BMAT ), ISTART,
     $                       ISTART, DESCB, LINFO )
               IF( LINFO.GE.1 ) THEN
*
*                Matrix may not be positive definite.
*
                  INFO = ( JSTARTB-1 ) + ( ISTART-1 ) + LINFO
                  WORK( 1 ) = INEED
                  RETURN
               ENDIF
*
*       Update rest of row/column.
*
*        B( (iend+1):isizeB, istart:iend) =...
*           B( (iend+1):isizeB, istart:iend) /
*                   B(istart:iend,istart:iend)'
*
               MM = ISIZEB - ( IEND+1 ) + 1
               NN = ISIZE
               ALPHA = ONE
               HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
               IF( HASWORK ) THEN
                  CALL PZTRSM( 'Rightside', 'LowerTriangular', 'C',
     $                         'NonunitDiagonal', MM, NN, ALPHA,
     $                         WORK( IP_BMAT ), ISTART, ISTART, DESCB,
     $                         WORK( IP_BMAT ), ( IEND+1 ), ISTART,
     $                         DESCB )
               ENDIF
*
*  Symmetric update of diagonal block.
*
*  B( (iend+1):jsizeB, (iend+1):jsizeB ) = ...
*       B( (iend+1):jsizeB, (iend+1):jsizeB ) -
*          B( (iend+1):jsizeB, istart:iend ) *
*               B((iend+1):jsizeB, istart:iend )'
*
               NN = JSIZEB - ( IEND+1 ) + 1
               KK = ISIZE
               ALPHA = -ONE
               BETA = ONE
               HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 )
               IF( HASWORK ) THEN
                  CALL PZHERK( 'Lower', 'NoTranspose', NN, KK,
     $                         DBLE( ALPHA ), WORK( IP_BMAT ),
     $                         ( IEND+1 ), ISTART, DESCB, DBLE( BETA ),
     $                         WORK( IP_BMAT ), ( IEND+1 ), ( IEND+1 ),
     $                         DESCB )
               ENDIF
*
* Update the rest
*
* B( (jsizeB+1):isizeB, (iend+1):jsizeB ) <-
*    B( (jsizeB+1):isizeB, (iend+1):jsizeB )
*           B((jsizeB+1):isizeB,  istart:iend )*
*               B((iend+1):jsizeB, istart:iend )'
*
               MM = ISIZEB - ( JSIZEB+1 ) + 1
               NN = JSIZEB - ( IEND+1 ) + 1
               KK = ISIZE
               ALPHA = -ONE
               BETA = ONE
               HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 ) .AND. ( KK.GE.1 )
               IF( HASWORK ) THEN
                  CALL PZGEMM( 'NoTranspose', 'C', MM, NN, KK, ALPHA,
     $                         WORK( IP_BMAT ), ( JSIZEB+1 ), ISTART,
     $                         DESCB, WORK( IP_BMAT ), ( IEND+1 ),
     $                         ISTART, DESCB, BETA, WORK( IP_BMAT ),
     $                         ( JSIZEB+1 ), ( IEND+1 ), DESCB )
               ENDIF
   60       CONTINUE
   70       CONTINUE
* end do istart
         ENDIF
*
*     Write out B.
*
*    C(istartB:iendB, jstartB:jendB ) = B;
*
*
         INFO = 0
         IIC = ( IC-1 ) + ISTARTB
         JJC = ( JC-1 ) + JSTARTB
         CALL ZLAWRITE( IODEV, ISIZEB, JSIZEB, IIC, JJC,
     $                  WORK( IP_BMAT ), 1, 1, DESCB, INFO )
         CALL ASSERT( INFO.EQ.0,
     $                '** PFxLCHFACT: LAWRITE returns info = ', INFO )
         JSTARTB = JENDB + 1
         GOTO 30
      ENDIF
   80 CONTINUE
* end while
      WORK( 1 ) = INEED
      INFO = 0
      RETURN
      END
