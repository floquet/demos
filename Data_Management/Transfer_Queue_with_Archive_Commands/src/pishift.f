      SUBROUTINE PISHIFT( M, N, A, IA, JA, DESCA, B, IB, JB, DESCB,
     $                    INFO )
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
* Perform an aligned shift or copy of a M by N submatrix
* from  A(ia:ia+m-1, ja:ja+n-1) to B( ib:ib+m-1,jb:jb+n-1).
*
* The routine assumes  the following are all satisfied:
*
* descA(CTXT_) == descB(CTXT_)
* descA(MB_) ==  descB(MB_)
* descA(NB_) ==  descB(NB_)
* mod(ia-1,MB) == mod(ib-1,MB)
* mod(ja-1,NB) == mod(jb-1,NB)
*
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DT_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1 )
      INTEGER            CTXT_, MB_, NB_
      PARAMETER          ( CTXT_ = 2, MB_ = 5, NB_ = 6 )
      INTEGER            RSRC_, CSRC_, LLD_
      PARAMETER          ( RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Scalar Arguments ..
      INTEGER            IA, IB, INFO, JA, JB, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            A( * ), B( * ), DESCA( DLEN_ ), DESCB( DLEN_ )
*     ..
*     .. Local Scalars ..
      LOGICAL            ALL_LOCAL, ISVALID
      INTEGER            DESTID, ICOLA, ICOLB, ICTXT, IDEST, IDEST_COL,
     $                   IDEST_ROW, IDUMMY, IROWA, IROWB, ISRC,
     $                   ISRC_COL, ISRC_ROW, J, LCINDXA, LCINDXB,
     $                   LRINDXA, LRINDXB, MB, MMA, MMB, MYID, MYPCOL,
     $                   MYPROW, NB, NNA, NNB, NPCOL, NPROW, SAVEDT
*     ..
*     .. External Functions ..
      INTEGER            BLACS_PNUM, NUMROC2
      EXTERNAL           BLACS_PNUM, NUMROC2
*     ..
*     .. External Subroutines ..
      EXTERNAL           ASSERT, BLACS_GRIDINFO, CHK1MAT, ICOPY,
     $                   IGERV2D, IGESD2D, INFOG2L, PCHK1MAT, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
* Perform error checking.
*
      INFO = 0
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYPROW,
     $                     MYPCOL )
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 6*100+CTXT_ )
         RETURN
      ENDIF
      ICTXT = DESCB( CTXT_ )
      CALL BLACS_GRIDINFO( DESCB( CTXT_ ), NPROW, NPCOL, MYPROW,
     $                     MYPCOL )
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 10*100+CTXT_ )
         RETURN
      ENDIF
      SAVEDT = DESCA( DT_ )
      DESCA( DT_ ) = BLOCK_CYCLIC_2D
      CALL CHK1MAT( M, 1, N, 2, IA, JA, DESCA, 6, INFO )
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PISHIFT', 6 )
         RETURN
      ENDIF
      IDUMMY = 1
      CALL PCHK1MAT( M, 1, N, 2, IA, JA, DESCA, 6, IDUMMY, IDUMMY,
     $               IDUMMY, INFO )
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PISHIFT', 6 )
         RETURN
      ENDIF
      DESCA( DT_ ) = SAVEDT
      SAVEDT = DESCB( DT_ )
      DESCB( DT_ ) = BLOCK_CYCLIC_2D
      CALL CHK1MAT( M, 1, N, 2, IB, JB, DESCB, 10, INFO )
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCB( CTXT_ ), 'PISHIFT', 10 )
         RETURN
      ENDIF
      IDUMMY = 1
      CALL PCHK1MAT( M, 1, N, 2, IB, JB, DESCB, 10, IDUMMY, IDUMMY,
     $               IDUMMY, INFO )
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCB( CTXT_ ), 'PISHIFT', 10 )
         RETURN
      ENDIF
      DESCB( DT_ ) = SAVEDT
      MB = DESCA( MB_ )
      NB = DESCA( NB_ )
      IF( DESCA( CTXT_ ).NE.DESCB( CTXT_ ) ) THEN
         INFO = -( 6*100+CTXT_ )
         CALL PXERBLA( DESCA( CTXT_ ), 'PISHIFT', 6 )
      ENDIF
      ISVALID = ( 1.LE.MB ) .AND. ( DESCA( MB_ ).EQ.DESCB( MB_ ) )
      IF( .NOT.ISVALID ) THEN
         INFO = -( 6*100+MB_ )
         CALL PXERBLA( DESCA( CTXT_ ), 'PISHIFT', 6 )
         RETURN
      ENDIF
      ISVALID = ( 1.LE.NB ) .AND. ( DESCA( NB_ ).EQ.DESCB( NB_ ) )
      IF( .NOT.ISVALID ) THEN
         INFO = -( 6*100+NB_ )
         CALL PXERBLA( DESCA( CTXT_ ), 'PISHIFT', 6 )
         RETURN
      ENDIF
      IF( MOD( IA-1, MB ).NE.MOD( IB-1, MB ) ) THEN
         INFO = -4
         CALL PXERBLA( DESCA( CTXT_ ), 'PISHIFT', ( -INFO ) )
         RETURN
      ENDIF
      IF( MOD( JA-1, NB ).NE.MOD( JB-1, MB ) ) THEN
         INFO = -5
         CALL PXERBLA( DESCA( CTXT_ ), 'PISHIFT', ( -INFO ) )
         RETURN
      ENDIF
*
*  All seems to check out ok.
*
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYPROW, MYPCOL,
     $              LRINDXA, LCINDXA, IROWA, ICOLA )
      CALL INFOG2L( IB, JB, DESCB, NPROW, NPCOL, MYPROW, MYPCOL,
     $              LRINDXB, LCINDXB, IROWB, ICOLB )
      MMA = NUMROC2( M, IA, MB, MYPROW, DESCA( RSRC_ ), NPROW )
      NNA = NUMROC2( N, JA, NB, MYPCOL, DESCA( CSRC_ ), NPCOL )
      MMB = NUMROC2( M, IB, MB, MYPROW, DESCB( RSRC_ ), NPROW )
      NNB = NUMROC2( N, JB, NB, MYPCOL, DESCB( CSRC_ ), NPCOL )
      CALL ASSERT( MMA.GE.0, '** PxSHIFT: invalid mmA = ', MMA )
      CALL ASSERT( MMB.GE.0, '** PxSHIFT: invalid mmB = ', MMB )
      CALL ASSERT( NNA.GE.0, '** PxSHIFT: invalid nnA = ', NNA )
      CALL ASSERT( NNB.GE.0, '** PxSHIFT: invalid nnB = ', NNB )
      ALL_LOCAL = ( ( IROWA.EQ.IROWB ) .AND. ( ICOLA.EQ.ICOLB ) )
      IF( ALL_LOCAL ) THEN
*
*          No communication is required.
*          just use copy.
*
         CALL ASSERT( MMA.EQ.MMB, '** PSHIFT: mmA != mmB ', MMA )
         CALL ASSERT( NNA.EQ.NNB, '** PSHIFT: nnA != nnB ', NNA )
         ISRC = LRINDXA + ( LCINDXA-1 )*DESCA( LLD_ )
         IDEST = LRINDXB + ( LCINDXB-1 )*DESCB( LLD_ )
         DO 10 J = 1, NNA
            CALL ICOPY( MMA, A( ISRC ), 1, B( IDEST ), 1 )
            ISRC = ISRC + DESCA( LLD_ )
            IDEST = IDEST + DESCB( LLD_ )
   10    CONTINUE
   20    CONTINUE
      ELSE
*
*               Need communication.
*               Assume sufficient send buffer space to avoid deadlock.
*
         IDEST_ROW = MOD( 2*NPROW+( IROWB-IROWA )+MYPROW, NPROW )
         IDEST_COL = MOD( 2*NPCOL+( ICOLB-ICOLA )+MYPCOL, NPCOL )
         ISRC_ROW = MOD( 2*NPROW+MYPROW-( IROWB-IROWA ), NPROW )
         ISRC_COL = MOD( 2*NPCOL+MYPCOL-( ICOLB-ICOLA ), NPCOL )
         DESTID = BLACS_PNUM( DESCA( CTXT_ ), IDEST_ROW, IDEST_COL )
         MYID = BLACS_PNUM( DESCA( CTXT_ ), MYPROW, MYPCOL )
         ISRC = LRINDXA + ( LCINDXA-1 )*DESCA( LLD_ )
         IDEST = LRINDXB + ( LCINDXB-1 )*DESCB( LLD_ )
         IF( ( MMA.GE.1 ) .AND. ( NNA.GE.1 ) ) THEN
            CALL IGESD2D( DESCA( CTXT_ ), MMA, NNA, A( ISRC ),
     $                    DESCA( LLD_ ), IDEST_ROW, IDEST_COL )
         ENDIF
         IF( ( MMB.GE.1 ) .AND. ( NNB.GE.1 ) ) THEN
            CALL IGERV2D( DESCB( CTXT_ ), MMB, NNB, B( IDEST ),
     $                    DESCB( LLD_ ), ISRC_ROW, ISRC_COL )
         ENDIF
      ENDIF
      INFO = 0
      RETURN
      END
