      SUBROUTINE PFDCOPYTAU( M, IA, TAUA, DESCA, IB, TAUB, DESCB )
*
*
*  -- ScaLAPACK auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     Oct 10, 1996
*
*
*   Purpose:
*   ========
*
*   The 'tau' vector is tied to the distribution of
*   a ScaLAPACK matrix.
*
*   PFCOPYTAU copies the tau vector tied to matrix A,
*   to another vector consistent with matrix B.
*
*
*     .. Parameters ..
      INTEGER            DLEN_
      PARAMETER          ( DLEN_ = 9 )
      INTEGER            CTXT_, MB_, NB_
      PARAMETER          ( CTXT_ = 2, MB_ = 5, NB_ = 6 )
      INTEGER            RSRC_, CSRC_
      PARAMETER          ( RSRC_ = 7, CSRC_ = 8 )
      INTEGER            NDIM
      PARAMETER          ( NDIM = 256 )
      LOGICAL            USEMAX
      PARAMETER          ( USEMAX = .false. )
*     ..
*     .. Scalar Arguments ..
      INTEGER            IA, IB, M
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( DLEN_ ), DESCB( DLEN_ )
      DOUBLE PRECISION   TAUA( * ), TAUB( * )
*     ..
*     .. Local Scalars ..
      LOGICAL            ISMINE, ISVALID, LOCALCOPY, SAMELEN,
     $                   SAMEOFFSET, SAMEPROC
      CHARACTER          SCOPE, TOP
      INTEGER            CA, CDEST, CONTXTA, CONTXTB, CSRCA, CSRCB,
     $                   GINDX, GINDXA, GINDXB, I, IBEND, IBSTART,
     $                   IDEST, IDEST_END, IEND, IPROCA, IPROCB, IROW,
     $                   ISIZE, ISRC, ISRC_END, ISTART, JCOL, LENA,
     $                   LENB, LINDX, LOCQA, LOCQB, MYPCOLA, MYPCOLB,
     $                   MYPROWA, MYPROWB, NB_A, NB_B, NPCOLA, NPCOLB,
     $                   NPROWA, NPROWB, RA, RCFLAG, RDEST, RSRCA, RSRCB
      DOUBLE PRECISION   DZERO
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   GTAU( NDIM )
*     ..
*     .. External Functions ..
      INTEGER            INDXG2L, INDXG2P, INDXL2G, NUMROC
      EXTERNAL           INDXG2L, INDXG2P, INDXL2G, NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DFILL, DGAMX2D, DGSUM2D,
     $                   INFOG1L
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
      DZERO = DBLE( 0 )
      CONTXTA = DESCA( CTXT_ )
      CONTXTB = DESCB( CTXT_ )
      CALL BLACS_GRIDINFO( CONTXTA, NPROWA, NPCOLA, MYPROWA, MYPCOLA )
      CALL BLACS_GRIDINFO( CONTXTB, NPROWB, NPCOLB, MYPROWB, MYPCOLB )
      CSRCA = DESCA( CSRC_ )
      RSRCA = DESCA( RSRC_ )
      CSRCB = DESCB( CSRC_ )
      RSRCB = DESCB( RSRC_ )
*
*  Check for special case where
*  everything is aligned.
*
      NB_A = DESCA( NB_ )
      NB_B = DESCB( NB_ )
      IPROCA = INDXG2P( IA, NB_A, MYPCOLA, CSRCA, NPCOLA )
      IPROCB = INDXG2P( IB, NB_B, MYPCOLB, CSRCB, NPCOLB )
      SAMEPROC = ( IPROCA.EQ.IPROCB ) .AND. ( RSRCA.EQ.RSRCB ) .AND.
     $           ( CSRCA.EQ.CSRCB ) .AND. ( NPROWA.EQ.NPROWB ) .AND.
     $           ( NPCOLA.EQ.NPCOLB )
      SAMEOFFSET = ( NB_A.EQ.NB_B ) .AND.
     $             ( MOD( IA, NB_A ).EQ.MOD( IB, NB_B ) )
      CALL INFOG1L( IA, NB_A, NPCOLA, MYPCOLA, CSRCA, ISRC, IPROCA )
      CALL INFOG1L( IB, NB_B, NPCOLB, MYPCOLB, CSRCB, IDEST, IPROCB )
      CALL INFOG1L( IA+M-1, NB_A, NPCOLA, MYPCOLA, CSRCA, ISRC_END,
     $              IPROCA )
      CALL INFOG1L( IB+M-1, NB_B, NPCOLB, MYPCOLB, CSRCB, IDEST_END,
     $              IPROCB )
      IF( IPROCA.NE.MYPCOLA ) THEN
         ISRC_END = ISRC_END - 1
      ENDIF
      IF( IPROCB.NE.MYPCOLB ) THEN
         IDEST_END = IDEST_END - 1
      ENDIF
      LENA = MAX( 0, ISRC_END-ISRC+1 )
      LENB = MAX( 0, IDEST_END-IDEST+1 )
      SAMELEN = ( LENA.EQ.LENB )
      LOCALCOPY = ( ( SAMEPROC .AND. SAMEOFFSET ) .AND. SAMELEN )
*
*  debug only
*
      LOCALCOPY = .false.
      IF( LOCALCOPY ) THEN
*
*         Everything is aligned and local,
*         so just perform copy.
*
         IF( LENA.GE.1 ) THEN
            LOCQA = NUMROC( IA+M-1, NB_A, MYPCOLA, CSRCA, NPCOLA )
            LOCQB = NUMROC( IB+M-1, NB_B, MYPCOLB, RSRCB, NPCOLB )
            ISVALID = ( 1.LE.ISRC ) .AND. ( ISRC.LE.ISRC_END ) .AND.
     $                ( ISRC_END.LE.LOCQA ) .AND. ( 1.LE.IDEST ) .AND.
     $                ( IDEST.LE.IDEST_END ) .AND.
     $                ( IDEST_END.LE.LOCQB )
            IF( .NOT.ISVALID ) THEN
               WRITE( *, FMT = 9999 )MYPROWA, MYPCOLA, M, IA, IB,
     $            LOCQA, ISRC, ISRC_END, LOCQB, IDEST, IDEST_END
 9999          FORMAT( ' myprowA,mypcolA ', 2( 1X, I5 ), ' m,ia,ib ',
     $               3( 1X, I5 ), ' LocqA, isrc,isrc_end ', 3( 1X, I5 ),
     $               ' LocqB, idest,idest_end ', 3( 1X, I5 ) )
               STOP '** error ** '
            ENDIF
            DO 10 I = 1, LENA
               GINDXA = INDXL2G( ( ISRC-1 )+I, NB_A, MYPCOLA, CSRCA,
     $                  NPCOLA )
               GINDXB = INDXL2G( ( IDEST-1 )+I, NB_B, MYPCOLB, RSRCB,
     $                  NPCOLB )
               ISVALID = ( IA.LE.GINDXA ) .AND.
     $                   ( GINDXA.LE.IA+M-1 ) .AND.
     $                   ( IB.LE.GINDXB ) .AND. ( GINDXB.LE.IB+M-1 )
               IF( .NOT.ISVALID ) THEN
                  WRITE( *, FMT = 9999 )MYPROWA, MYPCOLA, M, IA, IB,
     $               LOCQA, ISRC, ISRC_END, LOCQB, IDEST, IDEST_END
                  WRITE( *, FMT = * )'i,gindxA,gindxB ', I, GINDXA,
     $               GINDXB
                  STOP '** error ** '
               ENDIF
               TAUB( ( IDEST-1 )+I ) = TAUA( ( ISRC-1 )+I )
   10       CONTINUE
   20       CONTINUE
         ENDIF
      ELSE
         DO 70 ISTART = IA, IA + M - 1, NDIM
            IEND = MIN( ISTART+NDIM-1, IA+M-1 )
            ISIZE = IEND - ISTART + 1
            CALL DFILL( NDIM, DZERO, GTAU( 1 ), 1 )
            DO 30 GINDX = ISTART, IEND
* unique contribution by processor
*            holding the diagonal entry,
* the unique entry is essential for dsum to work
               JCOL = INDXG2P( GINDX, DESCA( NB_ ), MYPCOLA, CSRCA,
     $                NPCOLA )
               IROW = INDXG2P( GINDX, DESCA( MB_ ), MYPROWA, RSRCA,
     $                NPROWA )
               ISMINE = ( ( IROW.EQ.MYPROWA ) .AND.
     $                  ( JCOL.EQ.MYPCOLA ) )
               IF( ISMINE ) THEN
                  LINDX = INDXG2L( GINDX, DESCA( NB_ ), MYPCOLA, CSRCA,
     $                    NPCOLA )
                  GTAU( GINDX-ISTART+1 ) = TAUA( LINDX )
               ENDIF
   30       CONTINUE
   40       CONTINUE
* Note tau is between 1 and 2
            SCOPE = 'A'
            TOP = ' '
            RDEST = -1
            CDEST = -1
            RA = 0
            CA = 0
            RCFLAG = -1
            IF( USEMAX ) THEN
               CALL DGAMX2D( CONTXTA, SCOPE, TOP, ISIZE, 1, GTAU( 1 ),
     $                       ISIZE, RA, CA, RCFLAG, RDEST, CDEST )
            ELSE
* use summation, this is ok since only
* the processor holding the diagonal entry
* contributes to gtau(*)
               CALL DGSUM2D( CONTXTA, SCOPE, TOP, ISIZE, 1, GTAU( 1 ),
     $                       ISIZE, RDEST, CDEST )
            ENDIF
* copy from gtau back to local storage
            IBSTART = IB + ( ISTART-IA )
            IBEND = IB + ( IEND-IA )
            DO 50 GINDX = IBSTART, IBEND
               JCOL = INDXG2P( GINDX, DESCB( NB_ ), MYPCOLB, CSRCB,
     $                NPCOLB )
               ISMINE = ( JCOL.EQ.MYPCOLB )
               IF( ISMINE ) THEN
                  LINDX = INDXG2L( GINDX, DESCB( NB_ ), MYPCOLB, CSRCB,
     $                    NPCOLB )
                  TAUB( LINDX ) = GTAU( GINDX-IBSTART+1 )
               ENDIF
   50       CONTINUE
   60       CONTINUE
   70    CONTINUE
   80    CONTINUE
* end do istart
      ENDIF
      RETURN
      END
