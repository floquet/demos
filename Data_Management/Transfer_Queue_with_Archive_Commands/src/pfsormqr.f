      SUBROUTINE PFSORMQR( SIDE, TRANS, M, N, K, A, IA, JA, DESCA, TAU,
     $                     C, IC, JC, DESCC, WORK, LWORK, INFO )
*
*  Purpose
*  =======
*
*  PFORMQR overwrites the general real M-by-N distributed matrix
*  sub( C ) = C(IC:IC+M-1,JC:JC+N-1) with
*
*                      SIDE = 'L'            SIDE = 'R'
*  TRANS = 'N':      Q * sub( C )          sub( C ) * Q
*  TRANS = 'C':      Q**T * sub( C )       sub( C ) * Q**T
*
*  where Q is a orthogonal(unitary) distributed matrix defined as the
*  product of k elementary reflectors
*
*        Q = H(1) H(2) . . . H(k)
*
*  as returned by PDGEQRF. Q is of order M if SIDE = 'L' and of order N
*  if SIDE = 'R'.
*
*
*  The elementary reflectors are computed by PFGEQRF and stored
*  in an out-of-core matrix.
*
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DT_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1 )
      INTEGER            CTXT_, M_, N_, MB_, NB_
      PARAMETER          ( CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6 )
      INTEGER            RSRC_, CSRC_, LLD_
      PARAMETER          ( RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      INTEGER            IODEV_, SIZE_, DISK_BLOCK_CYCLIC_2D
      PARAMETER          ( IODEV_ = 10, SIZE_ = 11,
     $                   DISK_BLOCK_CYCLIC_2D = 601 )
      INTEGER            FLEN_
      PARAMETER          ( FLEN_ = 12 )
      LOGICAL            TRY_VARIABLE
      PARAMETER          ( TRY_VARIABLE = .false. )
*     ..
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            IA, IC, INFO, JA, JC, K, LWORK, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( FLEN_ ), DESCC( DLEN_ )
      REAL               A( * ), C( * ), TAU( * ), WORK( * )
*     ..
*     .. Local Scalars ..
      LOGICAL            ASIZEQUERY, FORWARD, HASWORK, ISTRANS, ISVALID,
     $                   LEFT, LWORKQUERY, NOTRANS, RIGHT
      INTEGER            ANEED, ASIZE, CSRC, IACOL, IAROW, ICOFF,
     $                   ICONTXT, ICPOS, ICROW, ICTXT, IDUMMY, IEND,
     $                   IFREE, IODEV, IPTAU, IROFF, ISIZE, ISTART,
     $                   JCPOS, JEND, JSIZE, JSTART, KK, LDA, LDATMP,
     $                   LINFO, LOCP, LOCQ, MB, MM, MMB, MQRNEED,
     $                   MYPCOL, MYPROW, M_A, NB, NFREE, NN, NNB, NPCOL,
     $                   NPROW, N_A, P0, Q0, RSRC, SAVEDT, TAUNEED
*     ..
*     .. Local Arrays ..
      INTEGER            DESCAIN( DLEN_ ), DESCATMP( DLEN_ )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            INDXG2P, NUMROC
      EXTERNAL           LSAME, INDXG2P, NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           ASSERT, BLACS_GRIDINFO, CHK1MAT, DESCINIT,
     $                   LAIO_INFO, PCHK1MAT, PFMAXSIZE, PFSCOPYTAU,
     $                   PSORMQR, PXERBLA, SLAREAD
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
      INFO = 0
      TAUNEED = 0
      ANEED = 0
      MQRNEED = 0
      LWORKQUERY = ( LWORK.EQ.-1 )
* error checking
      ICTXT = DESCC( CTXT_ )
      CALL BLACS_GRIDINFO( DESCC( CTXT_ ), NPROW, NPCOL, MYPROW,
     $                     MYPCOL )
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 14*100+CTXT_ )
         RETURN
      ENDIF
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYPROW,
     $                     MYPCOL )
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 9*100+CTXT_ )
         RETURN
      ENDIF
      HASWORK = ( M.GE.1 ) .AND. ( N.GE.1 ) .AND. ( K.GE.1 )
      IF( .NOT.HASWORK ) THEN
         RETURN
      ENDIF
      LEFT = LSAME( SIDE, 'L' )
      RIGHT = LSAME( SIDE, 'R' )
      ISVALID = ( LEFT .OR. RIGHT )
      IF( .NOT.ISVALID ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PFSORMQR', 1 )
         INFO = -1
         RETURN
      ENDIF
      ISTRANS = ( LSAME( TRANS, 'T' ) .OR. LSAME( TRANS, 'C' ) )
      NOTRANS = LSAME( TRANS, 'N' )
      ISVALID = ( ISTRANS .OR. NOTRANS )
      IF( .NOT.ISVALID ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PFSORMQR', 2 )
         INFO = -2
         RETURN
      ENDIF
      SAVEDT = DESCA( DT_ )
      DESCA( DT_ ) = BLOCK_CYCLIC_2D
      CALL CHK1MAT( M, 3, K, 5, IA, JA, DESCA, 9, INFO )
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PFSORMQR', 9 )
         RETURN
      ENDIF
      IDUMMY = 1
      CALL PCHK1MAT( M, 3, K, 5, IA, JA, DESCA, 9, IDUMMY, IDUMMY,
     $               IDUMMY, INFO )
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PFSORMQR', 9 )
         RETURN
      ENDIF
      DESCA( DT_ ) = SAVEDT
      SAVEDT = DESCC( DT_ )
      DESCC( DT_ ) = BLOCK_CYCLIC_2D
      CALL CHK1MAT( M, 3, N, 4, IC, JC, DESCC, 14, INFO )
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCC( CTXT_ ), 'PFSORMQR', 14 )
         RETURN
      ENDIF
      IDUMMY = 1
      CALL PCHK1MAT( M, 3, N, 4, IC, JC, DESCC, 14, IDUMMY, IDUMMY,
     $               IDUMMY, INFO )
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCC( CTXT_ ), 'PFSORMQR', 14 )
         RETURN
      ENDIF
      DESCC( DT_ ) = SAVEDT
      IODEV = DESCA( IODEV_ )
      ASIZE = DESCA( SIZE_ )
      ASIZEQUERY = ( ASIZE.EQ.-1 )
      IF( DESCA( DT_ ).NE.DISK_BLOCK_CYCLIC_2D ) THEN
         INFO = -( 9*100+DT_ )
         CALL PXERBLA( DESCA( CTXT_ ), 'PFSORMQR', 9 )
         RETURN
      ENDIF
      INFO = 0
      CALL LAIO_INFO( DESCA( IODEV_ ), MM, NN, MMB, NNB, MB, NB, CSRC,
     $                RSRC, ICONTXT )
      IF( MB.NE.DESCA( MB_ ) ) THEN
         INFO = -( 9*100+MB_ )
      ENDIF
      IF( NB.NE.DESCA( NB_ ) ) THEN
         INFO = -( 9*100+NB_ )
      ENDIF
      IF( ICONTXT.NE.DESCA( CTXT_ ) ) THEN
         INFO = -( 9*100+CTXT_ )
      ENDIF
      IF( DESCA( RSRC_ ).NE.RSRC ) THEN
         INFO = -( 9*100+RSRC_ )
      ENDIF
      IF( DESCA( CSRC_ ).NE.CSRC ) THEN
         INFO = -( 9*100+CSRC_ )
      ENDIF
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PFSORMQR', 9 )
         RETURN
      ENDIF
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYPROW,
     $                     MYPCOL )
      M_A = M
*
*     Note the use of p0,q0 to overestimate storage.
*
      P0 = MYPROW
      Q0 = MYPCOL
      LOCP = NUMROC( M_A, DESCA( MB_ ), MYPROW, P0, NPROW )
      LOCQ = NUMROC( MIN( K, NNB ), DESCA( NB_ ), MYPCOL, Q0, NPCOL )
      ANEED = LOCP*LOCQ
      IF( ANEED.GT.ASIZE ) THEN
         A( 1 ) = ANEED
         INFO = -( 9*100+SIZE_ )
         IF( .NOT.ASIZEQUERY ) THEN
            CALL PXERBLA( DESCA( CTXT_ ), 'PFSORMQR', 9 )
         ENDIF
         RETURN
      ENDIF
*
*  Allocate temp work storage
*
      N_A = -1
      CALL PFMAXSIZE( 'Column', ASIZE, M, N_A, DESCA( MB_ ),
     $                DESCA( NB_ ), DESCA( RSRC_ ), DESCA( CSRC_ ),
     $                DESCA( CTXT_ ), INFO )
      CALL ASSERT( INFO.EQ.0,
     $             '** PFORMQR: internal error, pfmaxsize returns info '
     $             , INFO )
      IF( N_A.NE.K ) THEN
         IF( N_A.GE.NNB ) THEN
            N_A = N_A - MOD( N_A, NNB )
         ENDIF
      ENDIF
      LOCQ = NUMROC( N_A, NB, MYPCOL, CSRC, NPCOL )
      IFREE = 1
      IPTAU = IFREE
      IFREE = IFREE + LOCQ
      NFREE = LWORK - IFREE + 1
      TAUNEED = LWORK - NFREE
      IF( NFREE.LE.0 ) THEN
         WORK( 1 ) = TAUNEED
         INFO = -16
         IF( .NOT.LWORKQUERY ) THEN
            CALL PXERBLA( DESCA( CTXT_ ), 'PFSORMQR', 16 )
         ENDIF
         RETURN
      ENDIF
      LINFO = 0
      CALL DESCINIT( DESCAIN, DESCA( M_ ), DESCA( N_ ), DESCA( MB_ ),
     $               DESCA( NB_ ), DESCA( RSRC_ ), DESCA( CSRC_ ),
     $               DESCA( CTXT_ ), DESCA( LLD_ ), LINFO )
      IF( LINFO.NE.0 ) THEN
         INFO = -( 9*100+ABS( LINFO ) )
         CALL PXERBLA( DESCA( CTXT_ ), 'PFSORMQR', 9 )
         RETURN
      ENDIF
      FORWARD = ( LEFT .AND. ISTRANS ) .OR. ( RIGHT .AND. NOTRANS )
      IF( FORWARD ) THEN
*
*  Try to use variable width panels
*
         JSTART = JA
   10    CONTINUE
         IF( JSTART.LE.JA+MIN( M, K )-1 ) THEN
            ISTART = ( JSTART-JA ) + IA
            IEND = IA + M - 1
            ISIZE = IEND - ISTART + 1
            M_A = ISIZE
            LOCP = NUMROC( M_A, MB, MYPROW, RSRC, NPROW )
            LDA = MAX( 1, LOCP )
            IF( TRY_VARIABLE ) THEN
               CALL PFMAXSIZE( 'Column', ASIZE, M_A, N_A, DESCA( MB_ ),
     $                         DESCA( NB_ ), DESCA( RSRC_ ),
     $                         DESCA( CSRC_ ), DESCA( CTXT_ ), INFO )
               CALL ASSERT( INFO.EQ.0,
     $             '** PFORMQR: internal error, pfmaxsize returns info '
     $                      , INFO )
            ENDIF
*
*            Want either
*               jend extends to end of matrix
*               jend = ja + min(m,k)-1
*
*            or
*               jend falls on a record boundary
*
*               mod(jend,nnb) == 0
*
            JEND = JSTART + N_A - 1
            JEND = MIN( JEND, JA+MIN( M, K )-1 )
            IF( JEND.NE.JA+MIN( M, K )-1 ) THEN
               IF( JEND.GE.NNB ) THEN
                  JEND = JEND - MOD( JEND, NNB )
               ENDIF
               IF( JEND.GE.NB ) THEN
                  JEND = JEND - MOD( JEND, NB )
               ENDIF
               JEND = MAX( JSTART, JEND )
            ENDIF
            JSIZE = JEND - JSTART + 1
            ISVALID = ( 1.LE.JSIZE ) .AND. ( JSIZE.LE.N_A )
            CALL ASSERT( ISVALID, '**  PFORMQR: internal error, jsize ',
     $                   JSIZE )
            IROFF = 1 + MOD( ( ISTART-1 ), MB )
            ICOFF = 1 + MOD( ( JSTART-1 ), NB )
            MM = ISIZE
            NN = N
            KK = JSIZE
            HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 ) .AND. ( KK.GE.1 )
            IF( HASWORK ) THEN
               ICPOS = ( ISTART-IA ) + IC
               JCPOS = JC
*
*                Satisfy PORMQR alignment constraint.
*
               IAROW = INDXG2P( ISTART, MB, MYPROW, RSRC, NPROW )
               IACOL = INDXG2P( JSTART, MB, MYPCOL, CSRC, NPCOL )
               ICROW = INDXG2P( ICPOS, DESCC( MB_ ), MYPROW,
     $                 DESCC( RSRC_ ), NPROW )
               IF( LSAME( SIDE, 'L' ) ) THEN
*
*                       Need iArow == iCrow
*
                  IAROW = ICROW
               ENDIF
               LOCP = NUMROC( ISIZE, DESCA( MB_ ), MYPROW, IAROW,
     $                NPROW )
               LOCQ = NUMROC( JSIZE, DESCA( NB_ ), MYPCOL, IACOL,
     $                NPCOL )
               LDATMP = MAX( 1, LOCP )
               CALL ASSERT( LOCP*LOCQ.LE.DESCA( SIZE_ ),
     $                      '**   PFORMQR: internal error, Asize = ',
     $                      DESCA( SIZE_ ) )
               CALL DESCINIT( DESCATMP, ISIZE, JSIZE, DESCA( MB_ ),
     $                        DESCA( NB_ ), IAROW, IACOL,
     $                        DESCA( CTXT_ ), LDATMP, INFO )
               CALL ASSERT( INFO.EQ.0,
     $                      '** PFOMQR: descinit returns info ', INFO )
               CALL SLAREAD( DESCA( IODEV_ ), ISIZE, JSIZE, ISTART,
     $                       JSTART, A, IROFF, ICOFF, DESCATMP, INFO )
               CALL PFSCOPYTAU( JSIZE, JSTART, TAU, DESCAIN, ICOFF,
     $                          WORK( IPTAU ), DESCATMP )
               CALL PSORMQR( SIDE, TRANS, MM, NN, KK, A, IROFF, ICOFF,
     $                       DESCATMP, WORK( IPTAU ), C, ICPOS, JCPOS,
     $                       DESCC, WORK( IFREE ), NFREE, INFO )
               MQRNEED = MAX( MQRNEED, INT( WORK( IFREE ) ) )
               IF( INFO.EQ.-16 ) THEN
*
*                       insufficient work space
*
                  INFO = -16
                  WORK( 1 ) = TAUNEED + MQRNEED
                  RETURN
               ENDIF
            ENDIF
            JSTART = JEND + 1
            GOTO 10
         ENDIF
   20    CONTINUE
      ELSE
*
*  Backward traversal, more complicated,
*  so use fixed width pannels.
*
         M_A = M
         LOCP = NUMROC( M_A, MB, MYPROW, RSRC, NPROW )
         LDA = MAX( 1, LOCP )
         CALL PFMAXSIZE( 'Column', ASIZE, M_A, N_A, DESCA( MB_ ),
     $                   DESCA( NB_ ), DESCA( RSRC_ ), DESCA( CSRC_ ),
     $                   DESCA( CTXT_ ), INFO )
         CALL ASSERT( INFO.EQ.0,
     $             '** PFORMQR: internal error, pfmaxsize returns info '
     $                , INFO )
         JEND = JA + MIN( M, K ) - 1
   30    CONTINUE
         IF( JEND.GE.JA ) THEN
            JSTART = MAX( JA, JEND-N_A+1 )
*
*               Want either jstart is at beginning of matrix
*               or fall on a record boundary
*               mod(jstart,nnb) == 1
*
*
            IF( JSTART.NE.JA ) THEN
               JSTART = MAX( JEND-N_A+1, INT( JSTART / NNB )*NNB+1 )
            ENDIF
            JSIZE = JEND - JSTART + 1
            ISTART = ( JSTART-JA ) + IA
            IEND = IA + M - 1
            ISIZE = IEND - ISTART + 1
            IROFF = 1 + MOD( ( ISTART-1 ), MB )
            ICOFF = 1 + MOD( ( JSTART-1 ), NB )
            MM = ISIZE
            KK = JSIZE
            NN = N
            HASWORK = ( MM.GE.1 ) .AND. ( KK.GE.1 ) .AND. ( NN.GE.1 )
            IF( HASWORK ) THEN
               ICPOS = ( ISTART-IA ) + IC
               JCPOS = JC
*
*                 Satisfy PORMQR alignment constraint.
*
               IAROW = INDXG2P( ISTART, MB, MYPROW, RSRC, NPROW )
               IACOL = INDXG2P( JSTART, NB, MYPCOL, CSRC, NPCOL )
               ICROW = INDXG2P( ICPOS, DESCC( MB_ ), MYPROW,
     $                 DESCC( RSRC_ ), NPROW )
               IF( LSAME( SIDE, 'L' ) ) THEN
* Need iArow == iCrow
                  IAROW = ICROW
               ENDIF
               LOCP = NUMROC( ISIZE, DESCA( MB_ ), MYPROW, IAROW,
     $                NPROW )
               LOCQ = NUMROC( JSIZE, DESCA( NB_ ), MYPCOL, IACOL,
     $                NPCOL )
               LDATMP = MAX( 1, LOCP )
               CALL ASSERT( LOCP*LOCQ.LE.DESCA( SIZE_ ),
     $                      '** PFORMQR: internal error, Asize = ',
     $                      DESCA( SIZE_ ) )
               CALL DESCINIT( DESCATMP, ISIZE, JSIZE, DESCA( MB_ ),
     $                        DESCA( NB_ ), IAROW, IACOL,
     $                        DESCA( CTXT_ ), LDATMP, INFO )
               CALL ASSERT( INFO.EQ.0,
     $                      '** PFOMQR: descinit returns info ', INFO )
               CALL SLAREAD( DESCA( IODEV_ ), ISIZE, JSIZE, ISTART,
     $                       JSTART, A, IROFF, ICOFF, DESCATMP, INFO )
               CALL PFSCOPYTAU( JSIZE, JSTART, TAU, DESCAIN, ICOFF,
     $                          WORK( IPTAU ), DESCATMP )
               CALL PSORMQR( SIDE, TRANS, MM, NN, KK, A, IROFF, ICOFF,
     $                       DESCATMP, WORK( IPTAU ), C, ICPOS, JCPOS,
     $                       DESCC, WORK( IFREE ), NFREE, INFO )
               MQRNEED = MAX( MQRNEED, INT( WORK( IFREE ) ) )
               IF( INFO.EQ.-16 ) THEN
* insufficient work space
                  INFO = -16
                  WORK( 1 ) = TAUNEED + MQRNEED
                  RETURN
               ENDIF
            ENDIF
            JEND = JSTART - 1
            GOTO 30
         ENDIF
   40    CONTINUE
      ENDIF
      A( 1 ) = ANEED
      WORK( 1 ) = TAUNEED + MQRNEED
      INFO = 0
      RETURN
      END
