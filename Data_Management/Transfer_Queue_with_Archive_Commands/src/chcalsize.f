      SUBROUTINE CHCALSIZE( DESCC, M, NNB, JSTARTB, LNFREE, GNCOLA,
     $                      GNCOLB )
* intent(out):: gncolA,gncolB
* calculate optimal sizes for A and B matrices
* used in out-of-core cholesky factorization
*     .. Parameters ..
      INTEGER            DLEN_
      PARAMETER          ( DLEN_ = 9 )
      INTEGER            CTXT_, MB_, NB_
      PARAMETER          ( CTXT_ = 2, MB_ = 5, NB_ = 6 )
      INTEGER            RSRC_, CSRC_
      PARAMETER          ( RSRC_ = 7, CSRC_ = 8 )
      LOGICAL            TRYALIGN
      PARAMETER          ( TRYALIGN = .true. )
*     ..
*     .. Scalar Arguments ..
      INTEGER            GNCOLA, GNCOLB, JSTARTB, LNFREE, M, NNB
*     ..
*     .. Array Arguments ..
      INTEGER            DESCC( DLEN_ )
*     ..
*     .. Local Scalars ..
      LOGICAL            ATEND, GOODENOUGH, NEEDPRINT
      CHARACTER          SCOPE, TOP
      INTEGER            ASIZE, BSIZE, CA, CDEST, COLSRC, CONTXT,
     $                   GNROWA, GNROWB, IEND, IGOOD, ISIZE, ISTART,
     $                   JENDB, JSIZEB, LDA, LDB, MB, MYID, MYPCOL,
     $                   MYPROW, NB, NCOLA, NCOLB, NFREE, NPCOL, NPROC,
     $                   NPROW, NROWA, NROWB, RA, RCFLAG, RDEST, ROWSRC
*     ..
*     .. External Functions ..
      INTEGER            NUMROC
      EXTERNAL           NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           ASSERT, BLACS_GRIDINFO, BLACS_PINFO,
     $                   CALCOLSIZE, IGAMN2D
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
      CALL BLACS_PINFO( MYID, NPROC )
      CALL BLACS_GRIDINFO( CONTXT, NPROW, NPCOL, MYPROW, MYPCOL )
      NEEDPRINT = ( 0.GE.3 ) .OR. ( ( 0.GE.2 ) .AND. ( MYID.EQ.0 ) )
      CONTXT = DESCC( CTXT_ )
      MB = DESCC( MB_ )
      NB = DESCC( NB_ )
      ROWSRC = DESCC( RSRC_ )
      COLSRC = DESCC( CSRC_ )
      NFREE = LNFREE
      SCOPE = 'A'
      TOP = ' '
      RA = -1
      CA = -1
      RCFLAG = -1
      RDEST = -1
      CDEST = -1
      CALL IGAMN2D( CONTXT, SCOPE, TOP, 1, 1, NFREE, 1, RA, CA, RCFLAG,
     $              RDEST, CDEST )
      ISTART = JSTARTB
      IEND = M
      ISIZE = IEND - ISTART + 1
      GNROWB = ISIZE
      NROWB = NUMROC( GNROWB, MB, MYPROW, ROWSRC, NPROW )
      LDB = NROWB
      GNROWA = M - JSTARTB + 1
      NROWA = NUMROC( GNROWA, MB, MYPROW, ROWSRC, NPROW )
      LDA = NROWA
      IF( JSTARTB.EQ.1 ) THEN
* special case,
* initial factorization,
* no updated needed so no need to allocate A
         GNCOLA = 0
         GNROWA = 0
         CALL CALCOLSIZE( NFREE, GNROWB, DESCC, GNCOLB )
         GNCOLB = MAX( 1, MIN( M, GNCOLB ) )
         SCOPE = 'A'
         TOP = ' '
         RA = -1
         CA = -1
         RCFLAG = -1
         RDEST = -1
         CDEST = -1
         CALL IGAMN2D( CONTXT, SCOPE, TOP, 1, 1, GNCOLB, 1, RA, CA,
     $                 RCFLAG, RDEST, CDEST )
         JENDB = JSTARTB + GNCOLB - 1
         ATEND = ( JENDB.EQ.M )
         IF( .NOT.ATEND ) THEN
            IF( TRYALIGN .AND. ( GNCOLB.GE.NNB ) ) THEN
               GNCOLB = GNCOLB - MOD( GNCOLB, NNB )
            ENDIF
            IF( ( GNCOLB.GE.NB ) ) THEN
               GNCOLB = GNCOLB - MOD( GNCOLB, NB )
            ENDIF
         ENDIF
         NCOLB = NUMROC( GNCOLB, NB, MYPCOL, COLSRC, NPCOL )
         BSIZE = LDB*NCOLB
         GOODENOUGH = ( BSIZE.LE.LNFREE )
         IF( GOODENOUGH ) THEN
            IGOOD = 1
         ELSE
            IGOOD = 0
         ENDIF
         SCOPE = 'A'
         TOP = ' '
         RA = -1
         CA = -1
         RCFLAG = -1
         RDEST = -1
         CDEST = -1
         CALL IGAMN2D( CONTXT, SCOPE, TOP, 1, 1, IGOOD, 1, RA, CA,
     $                 RCFLAG, RDEST, CDEST )
         GOODENOUGH = ( IGOOD.NE.0 )
         CALL ASSERT( GOODENOUGH, '** chcalsize: increase nfree by ',
     $                BSIZE-NFREE )
         RETURN
      ENDIF
* =============================
* first attempt,
* allocate A  to have sufficient space for I/O
* and then let B have the rest
* =============================
      GNCOLA = MIN( NNB, JSTARTB-1 )
      NCOLA = NUMROC( GNCOLA, NB, MYPCOL, COLSRC, NPCOL )
      ASIZE = LDA*NCOLA
* the rest for B
      CALL CALCOLSIZE( NFREE-ASIZE, GNROWB, DESCC, GNCOLB )
      GNCOLB = MAX( 1, MIN( M-JSTARTB+1, GNCOLB ) )
      SCOPE = 'A'
      TOP = ' '
      RA = -1
      CA = -1
      RCFLAG = -1
      RDEST = -1
      CDEST = -1
      CALL IGAMN2D( CONTXT, SCOPE, TOP, 1, 1, GNCOLB, 1, RA, CA, RCFLAG,
     $              RDEST, CDEST )
      JENDB = JSTARTB + GNCOLB - 1
      ATEND = ( JENDB.EQ.M )
      IF( .NOT.ATEND ) THEN
         IF( TRYALIGN .AND. ( GNCOLB.GE.NNB ) ) THEN
            GNCOLB = GNCOLB - MOD( GNCOLB, NNB )
         ENDIF
         IF( ( GNCOLB.GE.NB ) ) THEN
            GNCOLB = GNCOLB - MOD( GNCOLB, NB )
         ENDIF
      ENDIF
      NCOLB = NUMROC( GNCOLB, NB, MYPCOL, COLSRC, NPCOL )
      BSIZE = NCOLB*LDB
      JSIZEB = GNCOLB
      JENDB = JSTARTB + JSIZEB - 1
      GOODENOUGH = ( ASIZE+BSIZE.LE.LNFREE ) .AND.
     $             ( ( ATEND ) .OR. ( GNCOLB.GE.NB*NPCOL ) )
      IF( NEEDPRINT ) THEN
         WRITE( *, FMT = 9999 )MYID, GOODENOUGH, GNCOLA, GNCOLB, NFREE
 9999    FORMAT( 'myid ', I4, ' chcalsize: case1 ', ' igood ', L3,
     $         ' gncolA,gncolB ', 2( 1X, I6 ), ' nfree ', I7 )
         WRITE( *, FMT = 9998 )MYID, NROWA, NCOLA, NCOLB, ASIZE, BSIZE
 9998    FORMAT( ' myid ', I4, ' nrowA,ncolA,ncolB ', 3( 1X, I6 ),
     $         ' Asize,Bsize ', 2( 1X, I6 ) )
      ENDIF
* all processors must agree
      IF( GOODENOUGH ) THEN
         IGOOD = 1
      ELSE
         IGOOD = 0
      ENDIF
      SCOPE = 'A'
      TOP = ' '
      RA = -1
      CA = -1
      RCFLAG = -1
      RDEST = -1
      CDEST = -1
      CALL IGAMN2D( CONTXT, SCOPE, TOP, 1, 1, IGOOD, 1, RA, CA, RCFLAG,
     $              RDEST, CDEST )
      GOODENOUGH = ( IGOOD.NE.0 )
      IF( GOODENOUGH ) THEN
         RETURN
      ENDIF
* ========================================
* try gncolB = nb*npcol, gncolA = nb*npcol
* load balance, each processor has at least one block
* ========================================
      GNCOLB = MAX( 1, MIN( NB*NPCOL, M-JSTARTB+1 ) )
      NCOLB = NUMROC( GNCOLB, NB, MYPCOL, COLSRC, NPCOL )
      BSIZE = NCOLB*LDB
      JENDB = JSTARTB + GNCOLB - 1
      ATEND = ( JENDB.EQ.M )
      IF( .NOT.ATEND ) THEN
         IF( TRYALIGN .AND. ( GNCOLB.GE.NNB ) ) THEN
            GNCOLB = GNCOLB - MOD( GNCOLB, NNB )
         ENDIF
         IF( ( GNCOLB.GE.NB ) ) THEN
            GNCOLB = GNCOLB - MOD( GNCOLB, NB )
         ENDIF
      ENDIF
      GNCOLA = MIN( NB*NPCOL, JSTARTB-1 )
      NCOLA = NUMROC( GNCOLA, NB, MYPCOL, COLSRC, NPCOL )
      ASIZE = LDA*NCOLA
      GOODENOUGH = ( ASIZE+BSIZE.LE.LNFREE )
      IF( NEEDPRINT ) THEN
         WRITE( *, FMT = 9997 )MYID, GOODENOUGH, GNCOLA, GNCOLB, NFREE
 9997    FORMAT( 'myid ', I4, ' chcalsize: case2, ', ' igood ', L3,
     $         ' gncolA,gncolB ', 2( 1X, I6 ), ' nfree ', I7 )
         WRITE( *, FMT = 9998 )MYID, NROWA, NCOLA, NCOLB, ASIZE, BSIZE
      ENDIF
* all processors must agree
      IF( GOODENOUGH ) THEN
         IGOOD = 1
      ELSE
         IGOOD = 0
      ENDIF
      SCOPE = 'A'
      TOP = ' '
      RA = -1
      CA = -1
      RCFLAG = -1
      RDEST = -1
      CDEST = -1
      CALL IGAMN2D( CONTXT, SCOPE, TOP, 1, 1, IGOOD, 1, RA, CA, RCFLAG,
     $              RDEST, CDEST )
      GOODENOUGH = ( IGOOD.NE.0 )
      IF( GOODENOUGH ) THEN
         RETURN
      ENDIF
* ========================================
* gncolB is rather small,
* try allocating min amount of B for load balancing,
* then the rest for A
* ========================================
      GNCOLB = MAX( 1, MIN( NB*NPCOL, M-JSTARTB+1 ) )
      NCOLB = NUMROC( GNCOLB, NB, MYPCOL, COLSRC, NPCOL )
      BSIZE = NCOLB*LDB
      JENDB = JSTARTB + GNCOLB - 1
      ATEND = ( JENDB.EQ.M )
      CALL CALCOLSIZE( NFREE-BSIZE, GNROWA, DESCC, GNCOLA )
      SCOPE = 'A'
      TOP = ' '
      RA = -1
      CA = -1
      RCFLAG = -1
      RDEST = -1
      CDEST = -1
      CALL IGAMN2D( CONTXT, SCOPE, TOP, 1, 1, GNCOLA, 1, RA, CA, RCFLAG,
     $              RDEST, CDEST )
      IF( TRYALIGN .AND. ( GNCOLA.GE.NNB ) ) THEN
         GNCOLA = GNCOLA - MOD( GNCOLA, NNB )
      ENDIF
      IF( ( GNCOLA.GE.NB ) ) THEN
         GNCOLA = GNCOLA - MOD( GNCOLA, NB )
      ENDIF
      NCOLA = NUMROC( GNCOLA, NB, MYPCOL, COLSRC, NPCOL )
      ASIZE = LDA*NCOLA
      GOODENOUGH = ( ASIZE+BSIZE.LE.LNFREE )
* all processors must agree
      IF( NEEDPRINT ) THEN
         WRITE( *, FMT = 9996 )MYID, GOODENOUGH, GNCOLA, GNCOLB, NFREE
 9996    FORMAT( 'myid ', I4, ' chcalsize: case2b, ', ' igood ', L3,
     $         ' gncolA,gncolB ', 2( 1X, I6 ), ' nfree ', I7 )
         WRITE( *, FMT = 9998 )MYID, NROWA, NCOLA, NCOLB, ASIZE, BSIZE
      ENDIF
      IF( GOODENOUGH ) THEN
         IGOOD = 1
      ELSE
         IGOOD = 0
      ENDIF
      SCOPE = 'A'
      TOP = ' '
      RA = -1
      CA = -1
      RCFLAG = -1
      RDEST = -1
      CDEST = -1
      CALL IGAMN2D( CONTXT, SCOPE, TOP, 1, 1, IGOOD, 1, RA, CA, RCFLAG,
     $              RDEST, CDEST )
      GOODENOUGH = ( IGOOD.NE.0 )
      IF( GOODENOUGH ) THEN
         RETURN
      ENDIF
* memory is really tight,
* try splitting it 1/2 and 1/2 between A and B
      BSIZE = NFREE / 2
* what is my last column index
      CALL CALCOLSIZE( NFREE / 2, GNROWB, DESCC, GNCOLB )
      GNCOLB = MAX( 1, MIN( M-JENDB+1, GNCOLB ) )
      JENDB = JSTARTB + GNCOLB - 1
      ATEND = ( JENDB.EQ.M )
      IF( .NOT.ATEND ) THEN
         IF( TRYALIGN .AND. ( GNCOLB.GE.NNB ) ) THEN
            GNCOLB = GNCOLB - MOD( GNCOLB, NNB )
         ENDIF
         IF( ( GNCOLB.GE.NB ) ) THEN
            GNCOLB = GNCOLB - MOD( GNCOLB, NB )
         ENDIF
      ENDIF
      NCOLB = NUMROC( GNCOLB, NB, MYPCOL, COLSRC, NPCOL )
      BSIZE = LDB*NCOLB
      ASIZE = NFREE - BSIZE
      CALL CALCOLSIZE( ASIZE, GNROWA, DESCC, GNCOLA )
      GNCOLA = MAX( 1, MIN( JSTARTB-1, GNCOLA ) )
      SCOPE = 'A'
      TOP = ' '
      RA = -1
      CA = -1
      RCFLAG = -1
      RDEST = -1
      CDEST = -1
      CALL IGAMN2D( CONTXT, SCOPE, TOP, 1, 1, GNCOLA, 1, RA, CA, RCFLAG,
     $              RDEST, CDEST )
      NCOLA = NUMROC( GNCOLA, NB, MYPCOL, COLSRC, NPCOL )
      ASIZE = LDA*NCOLA
      GOODENOUGH = ( ASIZE+BSIZE.LE.LNFREE )
      IF( NEEDPRINT ) THEN
         WRITE( *, FMT = 9995 )MYID, GOODENOUGH, GNCOLA, GNCOLB, NFREE
 9995    FORMAT( ' myid ', I4, ' chcalsize: case3 ', ' igood ', L3,
     $         ' gncolA,gncolb ', 2( 1X, I6 ), ' nfree ', I7 )
         WRITE( *, FMT = 9998 )MYID, NROWA, NCOLA, NCOLB, ASIZE, BSIZE
      ENDIF
      IF( GOODENOUGH ) THEN
         IGOOD = 1
      ELSE
         IGOOD = 0
      ENDIF
      SCOPE = 'A'
      TOP = ' '
      RA = -1
      CA = -1
      RCFLAG = -1
      RDEST = -1
      CDEST = -1
      CALL IGAMN2D( CONTXT, SCOPE, TOP, 1, 1, IGOOD, 1, RA, CA, RCFLAG,
     $              RDEST, CDEST )
      GOODENOUGH = ( IGOOD.NE.0 )
      IF( GOODENOUGH ) THEN
         RETURN
      ENDIF
* getting desperate set gncolA = npcol, gncolB = npcol
      GNCOLA = NPCOL
      GNCOLB = NPCOL
      GNCOLA = MIN( GNCOLA, JSTARTB-1 )
      GNCOLB = MAX( 1, MIN( GNCOLB, M-JSTARTB+1 ) )
      JENDB = JSTARTB + GNCOLB - 1
      ATEND = ( JENDB.EQ.M )
      NCOLA = NUMROC( GNCOLA, NB, MYPCOL, COLSRC, NPCOL )
      ASIZE = LDA*NCOLA
      NCOLB = NUMROC( GNCOLB, NB, MYPCOL, COLSRC, NPCOL )
      BSIZE = NCOLB*LDB
      GOODENOUGH = ( ASIZE+BSIZE.LE.LNFREE )
      IF( NEEDPRINT ) THEN
         WRITE( *, FMT = 9994 )MYID, GOODENOUGH, GNCOLA, GNCOLB, NFREE
 9994    FORMAT( 'myid ', I4, ' chcalsize: case4  ', ' igood ', L3,
     $         ' gncolA,gncolB ', 2( 1X, I6 ), ' nfree ', I7 )
         WRITE( *, FMT = 9998 )MYID, NROWA, NCOLA, NCOLB, ASIZE, BSIZE
      ENDIF
      IF( GOODENOUGH ) THEN
         IGOOD = 1
      ELSE
         IGOOD = 0
      ENDIF
      SCOPE = 'A'
      TOP = ' '
      RA = -1
      CA = -1
      RCFLAG = -1
      RDEST = -1
      CDEST = -1
      CALL IGAMN2D( CONTXT, SCOPE, TOP, 1, 1, IGOOD, 1, RA, CA, RCFLAG,
     $              RDEST, CDEST )
      GOODENOUGH = ( IGOOD.NE.0 )
      IF( GOODENOUGH ) THEN
         RETURN
      ENDIF
* last attempt, gncolA = 1, gncolB = 1
      GNCOLA = 1
      GNCOLB = 1
      GNCOLA = MIN( GNCOLA, JSTARTB-1 )
      GNCOLB = MIN( GNCOLB, M-JSTARTB+1 )
      NCOLA = NUMROC( GNCOLA, NB, MYPCOL, COLSRC, NPCOL )
      ASIZE = LDA*NCOLA
      NCOLB = NUMROC( GNCOLB, NB, MYPCOL, COLSRC, NPCOL )
      BSIZE = NCOLB*LDB
      GOODENOUGH = ( ASIZE+BSIZE.LE.LNFREE )
      IF( NEEDPRINT ) THEN
         WRITE( *, FMT = 9993 )MYID, GOODENOUGH, GNCOLA, GNCOLB, NFREE
 9993    FORMAT( 'myid ', I4, ' chcalsize: case5 last attempt ',
     $         ' igood ', G3.3, ' gncolA,gncolB ', 2( 1X, I6 ),
     $         ' nfree ', I7 )
         WRITE( *, FMT = 9998 )MYID, NROWA, NCOLA, NCOLB, ASIZE, BSIZE
      ENDIF
      IF( GOODENOUGH ) THEN
         IGOOD = 1
      ELSE
         IGOOD = 0
      ENDIF
      SCOPE = 'A'
      TOP = ' '
      RA = -1
      CA = -1
      RCFLAG = -1
      RDEST = -1
      CDEST = -1
      CALL IGAMN2D( CONTXT, SCOPE, TOP, 1, 1, IGOOD, 1, RA, CA, RCFLAG,
     $              RDEST, CDEST )
      GOODENOUGH = ( IGOOD.NE.0 )
      IF( GOODENOUGH ) THEN
         RETURN
      ENDIF
* insufficient memory
      CALL ASSERT( ASIZE+BSIZE.LE.NFREE,
     $             '**  chcalsize: increase buffer size by at least ',
     $             ASIZE+BSIZE-NFREE+1 )
      RETURN
      END
