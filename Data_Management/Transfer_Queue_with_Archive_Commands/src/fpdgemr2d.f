      SUBROUTINE FPDGEMR2D( M, N, AMAT, IA, JA, DESCA, BMAT, IB, JB,
     $                      DESCB )
*     .. Parameters ..
      INTEGER            DLEN_
      PARAMETER          ( DLEN_ = 9 )
      INTEGER            CTXT_, M_, N_, MB_, NB_
      PARAMETER          ( CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6 )
      INTEGER            RSRC_, CSRC_, LLD_
      PARAMETER          ( RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      LOGICAL            TRYCOPY, TRYSHIFT
      PARAMETER          ( TRYCOPY = .true., TRYSHIFT = .false. )
*     ..
*     .. Scalar Arguments ..
      INTEGER            IA, IB, JA, JB, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( DLEN_ ), DESCB( DLEN_ )
      DOUBLE PRECISION   AMAT( * ), BMAT( * )
*     ..
*     .. Local Scalars ..
      LOGICAL            ISALIGNED, ISOK, ISVALID, USEINLINE
      INTEGER            CONTXTA, CONTXTB, CSRCA, CSRCB, I, ICOL, IDXA,
     $                   IDXB, II, INFO, IROW, ISIZEA, ISIZEB, J, JJ,
     $                   JSIZEA, JSIZEB, LCINDX, LDA, LDB, LRINDX,
     $                   MAXIA, MAXIB, MAXJA, MAXJB, MBA, MBB, MINIA,
     $                   MINIB, MINJA, MINJB, MYPCOLA, MYPCOLB, MYPROWA,
     $                   MYPROWB, NBA, NBB, NCOLA, NCOLB, NPCOLA,
     $                   NPCOLB, NPROWA, NPROWB, NROWA, NROWB, RSRCA,
     $                   RSRCB
      DOUBLE PRECISION   AIJ, ONE, ZERO
*     ..
*     .. External Functions ..
      LOGICAL            CHECKALIGN
      INTEGER            NUMROC
      EXTERNAL           CHECKALIGN, NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           ASSERT, BLACS_GRIDINFO, DCOPY, INFOG2L,
     $                   PDGEMR2DO, PDSHIFT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MOD
*     ..
*     .. Executable Statements ..
      ONE = DBLE( 1 )
      ZERO = DBLE( 0 )
      LDA = DESCA( LLD_ )
      LDB = DESCB( LLD_ )
      IF( TRYCOPY ) THEN
         ISALIGNED = CHECKALIGN( M, N, IA, JA, DESCA, IB, JB, DESCB )
      ELSE
         ISALIGNED = .false.
* treat as unaligned copy, use redist
      ENDIF
      IF( ISALIGNED ) THEN
         CONTXTA = DESCA( CTXT_ )
         CALL BLACS_GRIDINFO( CONTXTA, NPROWA, NPCOLA, MYPROWA,
     $                        MYPCOLA )
         NROWA = NUMROC( DESCA( M_ ), DESCA( MB_ ), MYPROWA,
     $           DESCA( RSRC_ ), NPROWA )
         NCOLA = NUMROC( DESCA( N_ ), DESCA( NB_ ), MYPCOLA,
     $           DESCA( CSRC_ ), NPCOLA )
         MBA = DESCA( MB_ )
         NBA = DESCA( NB_ )
         RSRCA = DESCA( RSRC_ )
         CSRCA = DESCA( CSRC_ )
         MINIA = NROWA + 1
         MAXIA = 0
         MINJA = NCOLA + 1
         MAXJA = 0
         CONTXTB = DESCB( CTXT_ )
         CALL BLACS_GRIDINFO( CONTXTB, NPROWB, NPCOLB, MYPROWB,
     $                        MYPCOLB )
         NROWB = NUMROC( DESCB( M_ ), DESCB( MB_ ), MYPROWB,
     $           DESCB( RSRC_ ), NPROWB )
         NCOLB = NUMROC( DESCB( N_ ), DESCB( NB_ ), MYPCOLB,
     $           DESCB( CSRC_ ), NPCOLB )
         MBB = DESCB( MB_ )
         NBB = DESCB( NB_ )
         RSRCB = DESCB( RSRC_ )
         CSRCB = DESCB( CSRC_ )
         MINIB = NROWB + 1
         MAXIB = 0
         MINJB = NCOLB + 1
         MAXJB = 0
*
*       Compute local extent
*
         CALL INFOG2L( IA, JA, DESCA, NPROWA, NPCOLA, MYPROWA, MYPCOLA,
     $                 LRINDX, LCINDX, IROW, ICOL )
         MINIA = LRINDX
         MINJA = LCINDX
         CALL INFOG2L( IA+M-1, JA+N-1, DESCA, NPROWA, NPCOLA, MYPROWA,
     $                 MYPCOLA, LRINDX, LCINDX, IROW, ICOL )
         IF( IROW.NE.MYPROWA ) THEN
            LRINDX = LRINDX - 1
         ENDIF
         IF( ICOL.NE.MYPCOLA ) THEN
            LCINDX = LCINDX - 1
         ENDIF
         MAXIA = LRINDX
         MAXJA = LCINDX
         CALL INFOG2L( IB, JB, DESCB, NPROWB, NPCOLB, MYPROWB, MYPCOLB,
     $                 LRINDX, LCINDX, IROW, ICOL )
         MINIB = LRINDX
         MINJB = LCINDX
         CALL INFOG2L( IB+M-1, JB+N-1, DESCB, NPROWB, NPCOLB, MYPROWB,
     $                 MYPCOLB, LRINDX, LCINDX, IROW, ICOL )
         IF( IROW.NE.MYPROWB ) THEN
            LRINDX = LRINDX - 1
         ENDIF
         IF( ICOL.NE.MYPCOLB ) THEN
            LCINDX = LCINDX - 1
         ENDIF
         MAXIB = LRINDX
         MAXJB = LCINDX
* perform memory to memory copy
         ISIZEA = ( MAXIA-MINIA ) + 1
         JSIZEA = ( MAXJA-MINJA ) + 1
         ISIZEB = ( MAXIB-MINIB ) + 1
         JSIZEB = ( MAXJB-MINJB ) + 1
         IF( ISIZEA.LE.0 ) THEN
            ISIZEA = 0
         ENDIF
         IF( JSIZEA.LE.0 ) THEN
            JSIZEA = 0
         ENDIF
         IF( ISIZEB.LE.0 ) THEN
            ISIZEB = 0
         ENDIF
         IF( JSIZEB.LE.0 ) THEN
            JSIZEB = 0
         ENDIF
         ISVALID = ( ISIZEA.EQ.ISIZEB ) .AND. ( JSIZEA.EQ.JSIZEB )
         IF( .NOT.ISVALID ) THEN
            WRITE( *, FMT = 9999 )IA, JA, IB, JB, MYPROWA, MYPCOLA,
     $         MYPROWB, MYPCOLB, MINIA, MAXIA, MINJA, MAXJA, MINIA,
     $         MAXIA, MINJA, MAXJA, ISIZEA, JSIZEA, ISIZEB, JSIZEB
 9999       FORMAT( ' ia,ja ', 2( 1X, I6 ), ' ib,jb ', 2( 1X, I6 ),
     $            ' myprowA,mypcolA ', 2( 1X, I4 ), ' myprowB,mypcolB ',
     $            2( 1X, I4 ), / ' miniA,maxiA ', 2( 1X, I6 ),
     $            ' minjA,maxjA ', 2( 1X, I6 ), ' miniB,maxiB ',
     $            2( 1X, I6 ), ' minjB,maxjB ', 2( 1X, I6 ),
     $            / ' isizeA,jsizeA ', 2( 1X, I6 ), ' isizeB,jsizeB ',
     $            2( 1X, I6 ) )
            STOP '** error in FPGEMR2D ** '
         ENDIF
* may have nothing to do on this processor
         IF( ( ISIZEA.GE.1 ) .AND. ( JSIZEA.GE.1 ) ) THEN
            USEINLINE = ( ISIZEA.LE.20 )
            IF( USEINLINE ) THEN
               DO 30 J = 1, JSIZEA
                  DO 10 I = 1, ISIZEA
                     II = MINIA + ( I-1 )
                     JJ = MINJA + ( J-1 )
                     AIJ = AMAT( II+( JJ-1 )*LDA )
                     II = MINIB + ( I-1 )
                     JJ = MINJB + ( J-1 )
                     BMAT( II+( JJ-1 )*LDB ) = AIJ
   10             CONTINUE
   20             CONTINUE
   30          CONTINUE
   40          CONTINUE
            ELSE
               I = 1
               DO 50 J = 1, JSIZEA
                  II = MINIA + ( I-1 )
                  JJ = MINJA + ( J-1 )
                  IDXA = II + ( JJ-1 )*LDA
                  II = MINIB + ( I-1 )
                  JJ = MINJB + ( J-1 )
                  IDXB = II + ( JJ-1 )*LDB
                  CALL DCOPY( ISIZEA, AMAT( IDXA ), 1, BMAT( IDXB ), 1 )
   50          CONTINUE
   60          CONTINUE
            ENDIF
         ENDIF
      ELSE
*
*         Check if we can do a shift.
*
         IF( TRYSHIFT ) THEN
            ISOK = ( DESCA( CTXT_ ).EQ.DESCB( CTXT_ ) ) .AND.
     $             ( DESCA( MB_ ).EQ.DESCB( MB_ ) ) .AND.
     $             ( DESCA( NB_ ).EQ.DESCB( NB_ ) ) .AND.
     $             ( MOD( IA-1, DESCA( MB_ ) ).EQ.
     $             MOD( IB-1, DESCB( MB_ ) ) ) .AND.
     $             ( MOD( JA-1, DESCA( NB_ ) ).EQ.
     $             MOD( JB-1, DESCB( NB_ ) ) )
         ELSE
            ISOK = .false.
         ENDIF
         IF( ISOK ) THEN
            INFO = 0
            CALL PDSHIFT( M, N, AMAT, IA, JA, DESCA, BMAT, IB, JB,
     $                    DESCB, INFO )
            CALL ASSERT( INFO.EQ.0,
     $                   '** FPxGEMR2D: PxSHIFT returns info= ', INFO )
         ELSE
            CALL PDGEMR2DO( M, N, AMAT, IA, JA, DESCA, BMAT, IB, JB,
     $                      DESCB )
         ENDIF
      ENDIF
      RETURN
      END
