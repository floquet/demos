      SUBROUTINE PFSLAPRNT( M, N, A, IA, JA, DESCA, IRPRNT, ICPRNT,
     $                      CMATNM, NOUT, WORK )
*
*
*  -- ScaLAPACK auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     Oct 10, 1996
*
*
*
* Purpose:
* ========
* PFxLAPRNT prints  out an out of core matrix
* similar in function to PxLAPRNT()
*
* Note:
* =====
*
* This routine uses PxELGET for simplicity and may be
* very inefficient for print out large submatrices.
*
*
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DT_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1 )
      INTEGER            CTXT_, MB_, NB_
      PARAMETER          ( CTXT_ = 2, MB_ = 5, NB_ = 6 )
      INTEGER            RSRC_, CSRC_
      PARAMETER          ( RSRC_ = 7, CSRC_ = 8 )
      INTEGER            IODEV_, SIZE_, DISK_BLOCK_CYCLIC_2D
      PARAMETER          ( IODEV_ = 10, SIZE_ = 11,
     $                   DISK_BLOCK_CYCLIC_2D = 601 )
      INTEGER            FLEN_
      PARAMETER          ( FLEN_ = 12 )
*     ..
*     .. Scalar Arguments ..
      CHARACTER*( * )    CMATNM
      INTEGER            IA, ICPRNT, IRPRNT, JA, M, N, NOUT
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( FLEN_ )
      REAL               A( * ), WORK( * )
*     ..
*     .. Local Scalars ..
      LOGICAL            HASWORK
      INTEGER            ANEED, ASIZE, CSRC, IAOFF, ICONTXT, ICTXT,
     $                   IDUMMY, IEND, IIA, INFO, IOCSRC, IODEV, IORSRC,
     $                   ISIZE, ISTART, JAOFF, JEND, JJA, JSIZE, JSTART,
     $                   LDA, LOCP, LOCQ, MB, MM, MMB, MYPCOL, MYPROW,
     $                   NB, NN, NNB, NPCOL, NPROW, RSRC, SAVEDT
*     ..
*     .. Local Arrays ..
      INTEGER            DESCATMP( DLEN_ )
*     ..
*     .. External Functions ..
      INTEGER            INDXG2P, NUMROC
      EXTERNAL           INDXG2P, NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           ASSERT, BLACS_GRIDINFO, CHK1MAT, DESCINIT,
     $                   LAIO_INFO, PCHK1MAT, PFSLAPRNT2, PXERBLA,
     $                   SLAREAD
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYPROW,
     $                     MYPCOL )
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 6*100+CTXT_ )
         RETURN
      ENDIF
      SAVEDT = DESCA( DT_ )
      DESCA( DT_ ) = BLOCK_CYCLIC_2D
      CALL CHK1MAT( M, 1, N, 2, IA, JA, DESCA, 6, INFO )
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PFSLAPRNT', 6 )
         RETURN
      ENDIF
      IDUMMY = 1
      CALL PCHK1MAT( M, 1, N, 2, IA, JA, DESCA, 6, IDUMMY, IDUMMY,
     $               IDUMMY, INFO )
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PFSLAPRNT', 6 )
         RETURN
      ENDIF
      DESCA( DT_ ) = SAVEDT
      IF( DESCA( DT_ ).NE.DISK_BLOCK_CYCLIC_2D ) THEN
         INFO = -( 6*100+DT_ )
         CALL PXERBLA( DESCA( CTXT_ ), 'PFSLAPRNT', 6 )
         RETURN
      ENDIF
      INFO = 0
      CALL LAIO_INFO( DESCA( IODEV_ ), MM, NN, MMB, NNB, MB, NB, CSRC,
     $                RSRC, ICONTXT )
      IF( MB.NE.DESCA( MB_ ) ) THEN
         INFO = -( 6*100+MB_ )
      ENDIF
      IF( NB.NE.DESCA( NB_ ) ) THEN
         INFO = -( 6*100+NB_ )
      ENDIF
      IF( ICONTXT.NE.DESCA( CTXT_ ) ) THEN
         INFO = -( 6*100+CTXT_ )
      ENDIF
      IF( DESCA( RSRC_ ).NE.RSRC ) THEN
         INFO = -( 6*100+RSRC_ )
      ENDIF
      IF( DESCA( CSRC_ ).NE.CSRC ) THEN
         INFO = -( 6*100+CSRC_ )
      ENDIF
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PFSLAPRNT', 6 )
         RETURN
      ENDIF
      IODEV = DESCA( IODEV_ )
      MB = DESCA( MB_ )
      NB = DESCA( NB_ )
*
*  Make temporary buffer A(1,1) aligned
*  for local I/O operations.
*
      IORSRC = RSRC
      IOCSRC = CSRC
      RSRC = INDXG2P( IA, MB, MYPROW, IORSRC, NPROW )
      CSRC = INDXG2P( JA, NB, MYPCOL, IOCSRC, NPCOL )
      LOCP = NUMROC( MMB, MB, MYPROW, RSRC, NPROW )
      LOCQ = NUMROC( NNB, NB, MYPCOL, CSRC, NPCOL )
      LDA = MAX( 1, LOCP )
      ANEED = LOCP*LOCQ
      ASIZE = DESCA( SIZE_ )
      IF( ASIZE.LT.ANEED ) THEN
         A( 1 ) = ANEED
         CALL PXERBLA( DESCA( CTXT_ ), 'PFSLAPRNT', 6 )
         RETURN
      ENDIF
*
*   Bring in mmb x nnb patches.
*
      CALL DESCINIT( DESCATMP, MMB, NNB, MB, MB, RSRC, CSRC,
     $               DESCA( CTXT_ ), LDA, INFO )
      CALL ASSERT( INFO.EQ.0, '** PFLAPRNT: descinit returns info ',
     $             INFO )
      DO 30 JSTART = 1, N, NNB
         JEND = MIN( N, JSTART+NNB-1 )
         JSIZE = JEND - JSTART + 1
         DO 10 ISTART = 1, M, MMB
            IEND = MIN( M, ISTART+MMB-1 )
            ISIZE = IEND - ISTART + 1
            HASWORK = ( ISTART.GE.1 ) .AND. ( JSTART.GE.1 )
            IF( HASWORK ) THEN
               INFO = 0
               IIA = ( IA-1 ) + ISTART
               JJA = ( JA-1 ) + JSTART
               CALL SLAREAD( IODEV, ISIZE, JSIZE, IIA, JJA, A, 1, 1,
     $                       DESCATMP, INFO )
               CALL ASSERT( INFO.EQ.0,
     $                '** PFLAPRNT: internal error,LAREAD returns info '
     $                      , INFO )
               IAOFF = IIA - 1
               JAOFF = JJA - 1
               CALL PFSLAPRNT2( ISIZE, JSIZE, IAOFF, JAOFF, A, 1, 1,
     $                          DESCATMP, IRPRNT, ICPRNT, CMATNM, NOUT,
     $                          WORK )
            ENDIF
   10    CONTINUE
   20    CONTINUE
   30 CONTINUE
   40 CONTINUE
      A( 1 ) = ANEED
      RETURN
      END
