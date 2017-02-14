* DO 30 loop corrected.  Was calculating GINDX = JA + ( I-1 )
* and GINDX = JB + ( I-1 ).  Changed to GINDX = JA + ( J-1 )
* and GINDX = JB + ( J-1 ).
* SBC, 4/13/04

      LOGICAL          FUNCTION CHECKALIGN( M, N, IA, JA, DESCA, IB, JB,
     $                 DESCB )
*     .. Parameters ..
      INTEGER            DLEN_
      PARAMETER          ( DLEN_ = 9 )
      INTEGER            CTXT_, MB_, NB_
      PARAMETER          ( CTXT_ = 2, MB_ = 5, NB_ = 6 )
      INTEGER            RSRC_, CSRC_
      PARAMETER          ( RSRC_ = 7, CSRC_ = 8 )
*     ..
*     .. Scalar Arguments ..
      INTEGER            IA, IB, JA, JB, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( DLEN_ ), DESCB( DLEN_ )
*     ..
*     .. Local Scalars ..
      LOGICAL            ISALIGNED, PROCSHAPEOK, SAMELENG, SAMEOFFSET,
     $                   SAMEPROC, SAMESHAPE, SAMEWIDTH
      INTEGER            CONTXTA, CONTXTB, CSRCA, CSRCB, GINDX, I, IHIA,
     $                   IHIB, ILOA, ILOB, IPCOLA, IPCOLB, IPROWA,
     $                   IPROWB, J, JHIA, JHIB, JLOA, JLOB, LCINDXA,
     $                   LCINDXB, LENA, LENB, LRINDXA, LRINDXB, MYIDA,
     $                   MYIDB, MYPCOLA, MYPCOLB, MYPROWA, MYPROWB,
     $                   NPCOLA, NPCOLB, NPROWA, NPROWB, PIDA, PIDB,
     $                   RSRCA, RSRCB, WIDTHA, WIDTHB
*     ..
*     .. External Functions ..
      INTEGER            BLACS_PNUM, INDXG2P
      EXTERNAL           BLACS_PNUM, INDXG2P
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, INFOG2L
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MOD
*     ..
*     .. Executable Statements ..
      CONTXTA = DESCA( CTXT_ )
      CONTXTB = DESCB( CTXT_ )
      CALL BLACS_GRIDINFO( CONTXTA, NPROWA, NPCOLA, MYPROWA, MYPCOLA )
      CALL BLACS_GRIDINFO( CONTXTB, NPROWB, NPCOLB, MYPROWB, MYPCOLB )
      CALL INFOG2L( IA, JA, DESCA, NPROWA, NPCOLA, MYPROWA, MYPCOLA,
     $              LRINDXA, LCINDXA, RSRCA, CSRCA )
      ILOA = LRINDXA
      JLOA = LCINDXA
      CALL INFOG2L( IB, JB, DESCB, NPROWB, NPCOLB, MYPROWB, MYPCOLB,
     $              LRINDXB, LCINDXB, RSRCB, CSRCB )
      ILOB = LRINDXB
      JLOB = LCINDXB
      CALL INFOG2L( ( IA-1 )+M, ( JA-1 )+N, DESCA, NPROWA, NPCOLA,
     $              MYPROWA, MYPCOLA, LRINDXA, LCINDXA, RSRCA, CSRCA )
      IF( MYPROWA.NE.RSRCA ) THEN
         LRINDXA = LRINDXA - 1
      ENDIF
      IF( MYPCOLA.NE.CSRCA ) THEN
         LCINDXA = LCINDXA - 1
      ENDIF
      IHIA = LRINDXA
      JHIA = LCINDXA
      CALL INFOG2L( ( IB-1 )+M, ( JB-1 )+N, DESCB, NPROWB, NPCOLB,
     $              MYPROWB, MYPCOLB, LRINDXB, LCINDXB, RSRCB, CSRCB )
      IF( MYPROWB.NE.RSRCB ) THEN
         LRINDXB = LRINDXB - 1
      ENDIF
      IF( MYPCOLB.NE.CSRCB ) THEN
         LCINDXB = LCINDXB - 1
      ENDIF
      IHIB = LRINDXB
      JHIB = LCINDXB
      LENA = MAX( 0, IHIA-ILOA+1 )
      LENB = MAX( 0, IHIB-ILOB+1 )
      WIDTHA = MAX( 0, JHIA-JLOA+1 )
      WIDTHB = MAX( 0, JHIB-JLOB+1 )
      SAMELENG = ( LENA.EQ.LENB )
      SAMEWIDTH = ( WIDTHA.EQ.WIDTHB )
      MYIDA = BLACS_PNUM( CONTXTA, MYPROWA, MYPCOLA )
      MYIDB = BLACS_PNUM( CONTXTB, MYPROWB, MYPCOLB )
      PROCSHAPEOK = ( ( NPROWA.EQ.NPROWB ) .AND. ( NPCOLA.EQ.NPCOLB ) )
      SAMEPROC = ( ( MYPROWA.EQ.MYPROWB ) .AND.
     $           ( MYPCOLA.EQ.MYPCOLB ) .AND. ( MYIDA.EQ.MYIDB ) )
      SAMESHAPE = ( ( DESCA( MB_ ).EQ.DESCB( MB_ ) ) .AND.
     $            ( DESCA( NB_ ).EQ.DESCB( NB_ ) ) ) .AND.
     $            ( SAMELENG .AND. SAMEWIDTH )
      SAMEOFFSET = ( ( MOD( IA, DESCA( MB_ ) ).EQ.MOD( IB,
     $             DESCB( MB_ ) ) ) .AND. ( MOD( JA,
     $             DESCA( NB_ ) ).EQ.MOD( JB, DESCB( NB_ ) ) ) )
      ISALIGNED = ( PROCSHAPEOK .AND. SAMEPROC .AND. SAMESHAPE .AND.
     $            SAMEOFFSET )
      IF( .NOT.ISALIGNED ) THEN
         CHECKALIGN = ( .false. )
         RETURN
      ENDIF
* extra check
      GINDX = JA
      IPCOLA = INDXG2P( GINDX, DESCA( NB_ ), MYPCOLA, DESCA( CSRC_ ),
     $         NPCOLA )
      GINDX = JB
      IPCOLB = INDXG2P( GINDX, DESCB( NB_ ), MYPCOLB, DESCB( CSRC_ ),
     $         NPCOLB )
      DO 10 I = 1, M
         GINDX = IA + ( I-1 )
         IPROWA = INDXG2P( GINDX, DESCA( MB_ ), MYPROWA, DESCA( RSRC_ ),
     $            NPROWA )
         GINDX = IB + ( I-1 )
         IPROWB = INDXG2P( GINDX, DESCB( MB_ ), MYPROWB, DESCB( RSRC_ ),
     $            NPROWB )
         PIDA = BLACS_PNUM( CONTXTA, IPROWA, IPCOLA )
         PIDB = BLACS_PNUM( CONTXTB, IPROWB, IPCOLB )
         ISALIGNED = ISALIGNED .AND. ( PIDA.EQ.PIDB )
         IF( .NOT.ISALIGNED ) THEN
            CHECKALIGN = ( .false. )
            RETURN
         ENDIF
   10 CONTINUE
   20 CONTINUE
* end do i
      IPROWA = INDXG2P( IA, DESCA( MB_ ), MYPROWA, DESCA( RSRC_ ),
     $         NPROWA )
      IPROWB = INDXG2P( IB, DESCB( MB_ ), MYPROWB, DESCB( RSRC_ ),
     $         NPROWB )
      DO 30 J = 1, N
         GINDX = JA + ( J-1 )
         IPCOLA = INDXG2P( GINDX, DESCA( NB_ ), MYPCOLA, DESCA( CSRC_ ),
     $            NPCOLA )
         GINDX = JB + ( J-1 )
         IPCOLB = INDXG2P( GINDX, DESCB( NB_ ), MYPCOLB, DESCB( CSRC_ ),
     $            NPCOLB )
         PIDA = BLACS_PNUM( CONTXTA, IPROWA, IPCOLA )
         PIDB = BLACS_PNUM( CONTXTB, IPROWB, IPCOLB )
         ISALIGNED = ISALIGNED .AND. ( PIDA.EQ.PIDB )
         IF( .NOT.ISALIGNED ) THEN
            CHECKALIGN = ( .false. )
            RETURN
         ENDIF
   30 CONTINUE
   40 CONTINUE
* end do j
* passed all tests
      CHECKALIGN = ( .true. )
      RETURN
      END
