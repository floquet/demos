      SUBROUTINE PFCGEMM( TRANSA, TRANSB, M, N, K, ALPHA, A, IA, JA,
     $                    DESCA, B, IB, JB, DESCB, BETA, C, IC, JC,
     $                    DESCC, INFO )
*
*
*  -- ScaLAPACK auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     Oct 10, 1996
*
*
*   Purpose:
*   =======
*
*   PFvGEMM performs a matrix-matrix update similar to PvGEMM where
*   at most one of A,B, or C is an out-of-core matrix.
*
*   The Goal is to bring in parts of out-of-core matrix only ONCE.
*   This may require bringing in complete row or complete column
*   instead of  rectangular patches.
*
*   See documentation on PBLAS routine PxGEMM()
*   for more details.
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
*     ..
*     .. Scalar Arguments ..
      CHARACTER          TRANSA, TRANSB
      INTEGER            IA, IB, IC, INFO, JA, JB, JC, K, M, N
      COMPLEX            ALPHA, BETA
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * ), DESCC( * )
      COMPLEX            A( * ), B( * ), C( * )
*     ..
*     .. Local Scalars ..
      LOGICAL            ALL_INCORE, AONDISK, ASIZEQUERY, BONDISK,
     $                   BSIZEQUERY, CONDISK, CSIZEQUERY, HASWORK,
     $                   ISTRANSA, ISTRANSB, ISVALID, NOTRANSA, NOTRANSB
      INTEGER            ANEED, ASIZE, BNEED, BSIZE, CNEED, CSIZE, CSRC,
     $                   ICOLA, ICOLB, ICOLC, ICONTXT, ICTXT, IDUMMY,
     $                   IEND, IERR, II, IIA, IIB, IIC, IIHI, IILO,
     $                   IINC, IODEV, IOFF, IROWA, IROWB, IROWC, ISIZE,
     $                   ISTART, JEND, JINC, JJ, JJA, JJB, JJC, JJHI,
     $                   JJLO, JOFF, JSIZE, JSTART, KK, LDA, LDB, LDC,
     $                   LOCP, LOCQ, MB, MM, MMB, MYPCOL, MYPROW, M_A,
     $                   M_B, NB, NDISK, NN, NNB, NPCOL, NPROW, N_A,
     $                   N_B, P0, Q0, RSRC, SAVEDT
      COMPLEX            LBETA, ONE
*     ..
*     .. Local Arrays ..
      INTEGER            DESCATMP( DLEN_ ), DESCBTMP( DLEN_ ),
     $                   DESCCTMP( DLEN_ )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL, INDXG2P, NUMROC
      EXTERNAL           LSAME, ICEIL, INDXG2P, NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           ASSERT, BLACS_GRIDINFO, CHK1MAT, CLAREAD,
     $                   CLAWRITE, DESCINIT, IGSUM2D, LAIO_INFO, PCGEMM,
     $                   PCHK1MAT, PCSCAL, PFMAX2SIZE, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, INT, MAX, MIN, MOD, REAL
*     ..
*     .. Executable Statements ..
*
*       Perform error checking
*
      ONE = CMPLX( REAL( 1 ) )
      INFO = 0
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYPROW,
     $                     MYPCOL )
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 10*100+CTXT_ )
         RETURN
      ENDIF
      ICTXT = DESCB( CTXT_ )
      CALL BLACS_GRIDINFO( DESCB( CTXT_ ), NPROW, NPCOL, MYPROW,
     $                     MYPCOL )
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 14*100+CTXT_ )
         RETURN
      ENDIF
      ICTXT = DESCC( CTXT_ )
      CALL BLACS_GRIDINFO( DESCC( CTXT_ ), NPROW, NPCOL, MYPROW,
     $                     MYPCOL )
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 16*100+CTXT_ )
         RETURN
      ENDIF
      ISVALID = ( LSAME( TRANSA, 'T' ) .OR. LSAME( TRANSA, 'C' ) .OR.
     $          LSAME( TRANSA, 'N' ) )
      IF( .NOT.ISVALID ) THEN
         INFO = -1
         CALL PXERBLA( DESCA( CTXT_ ), 'PFCGEMM', 1 )
         RETURN
      ENDIF
      ISVALID = ( LSAME( TRANSB, 'T' ) .OR. LSAME( TRANSB, 'C' ) .OR.
     $          LSAME( TRANSB, 'N' ) )
      IF( .NOT.ISVALID ) THEN
         INFO = -2
         CALL PXERBLA( DESCB( CTXT_ ), 'PFCGEMM', 2 )
         RETURN
      ENDIF
      NOTRANSA = LSAME( TRANSA, 'N' )
      ISTRANSA = ( .NOT.NOTRANSA )
      NOTRANSB = LSAME( TRANSB, 'N' )
      ISTRANSB = ( .NOT.NOTRANSB )
      IF( ISTRANSA ) THEN
*
*       A is k x m
*
         SAVEDT = DESCA( DT_ )
         DESCA( DT_ ) = BLOCK_CYCLIC_2D
         CALL CHK1MAT( K, 5, M, 3, IA, JA, DESCA, 10, INFO )
         IF( INFO.NE.0 ) THEN
            CALL PXERBLA( DESCA( CTXT_ ), 'PFCGEMM', 10 )
            RETURN
         ENDIF
         IDUMMY = 1
         CALL PCHK1MAT( K, 5, M, 3, IA, JA, DESCA, 10, IDUMMY, IDUMMY,
     $                  IDUMMY, INFO )
         IF( INFO.NE.0 ) THEN
            CALL PXERBLA( DESCA( CTXT_ ), 'PFCGEMM', 10 )
            RETURN
         ENDIF
         DESCA( DT_ ) = SAVEDT
      ELSE
*
*       A is m x k.
*
         SAVEDT = DESCA( DT_ )
         DESCA( DT_ ) = BLOCK_CYCLIC_2D
         CALL CHK1MAT( M, 3, K, 5, IA, JA, DESCA, 10, INFO )
         IF( INFO.NE.0 ) THEN
            CALL PXERBLA( DESCA( CTXT_ ), 'PFCGEMM', 10 )
            RETURN
         ENDIF
         IDUMMY = 1
         CALL PCHK1MAT( M, 3, K, 5, IA, JA, DESCA, 10, IDUMMY, IDUMMY,
     $                  IDUMMY, INFO )
         IF( INFO.NE.0 ) THEN
            CALL PXERBLA( DESCA( CTXT_ ), 'PFCGEMM', 10 )
            RETURN
         ENDIF
         DESCA( DT_ ) = SAVEDT
      ENDIF
      IF( ISTRANSB ) THEN
*
*       B is n x k
*
         SAVEDT = DESCB( DT_ )
         DESCB( DT_ ) = BLOCK_CYCLIC_2D
         CALL CHK1MAT( N, 4, K, 5, IB, JB, DESCB, 14, INFO )
         IF( INFO.NE.0 ) THEN
            CALL PXERBLA( DESCB( CTXT_ ), 'PFCGEMM', 14 )
            RETURN
         ENDIF
         IDUMMY = 1
         CALL PCHK1MAT( N, 4, K, 5, IB, JB, DESCB, 14, IDUMMY, IDUMMY,
     $                  IDUMMY, INFO )
         IF( INFO.NE.0 ) THEN
            CALL PXERBLA( DESCB( CTXT_ ), 'PFCGEMM', 14 )
            RETURN
         ENDIF
         DESCB( DT_ ) = SAVEDT
      ELSE
*
*       B is k x n
*
         SAVEDT = DESCB( DT_ )
         DESCB( DT_ ) = BLOCK_CYCLIC_2D
         CALL CHK1MAT( K, 5, N, 4, IB, JB, DESCB, 14, INFO )
         IF( INFO.NE.0 ) THEN
            CALL PXERBLA( DESCB( CTXT_ ), 'PFCGEMM', 14 )
            RETURN
         ENDIF
         IDUMMY = 1
         CALL PCHK1MAT( K, 5, N, 4, IB, JB, DESCB, 14, IDUMMY, IDUMMY,
     $                  IDUMMY, INFO )
         IF( INFO.NE.0 ) THEN
            CALL PXERBLA( DESCB( CTXT_ ), 'PFCGEMM', 14 )
            RETURN
         ENDIF
         DESCB( DT_ ) = SAVEDT
      ENDIF
*
*       C is m x n
*
      SAVEDT = DESCC( DT_ )
      DESCC( DT_ ) = BLOCK_CYCLIC_2D
      CALL CHK1MAT( M, 3, N, 4, IC, JC, DESCC, 16, INFO )
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCC( CTXT_ ), 'PFCGEMM', 16 )
         RETURN
      ENDIF
      IDUMMY = 1
      CALL PCHK1MAT( M, 3, N, 4, IC, JC, DESCC, 16, IDUMMY, IDUMMY,
     $               IDUMMY, INFO )
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCC( CTXT_ ), 'PFCGEMM', 16 )
         RETURN
      ENDIF
      DESCC( DT_ ) = SAVEDT
      AONDISK = ( DESCA( DT_ ).EQ.DISK_BLOCK_CYCLIC_2D )
      BONDISK = ( DESCB( DT_ ).EQ.DISK_BLOCK_CYCLIC_2D )
      CONDISK = ( DESCC( DT_ ).EQ.DISK_BLOCK_CYCLIC_2D )
*
*    Can support at most one out-of-core matrix.
*
      NDISK = 0
      IF( AONDISK ) THEN
         NDISK = NDISK + 1
      ENDIF
      IF( BONDISK ) THEN
         NDISK = NDISK + 1
      ENDIF
      IF( CONDISK ) THEN
         NDISK = NDISK + 1
      ENDIF
      IF( NDISK.EQ.0 ) THEN
*
*          All matrices are in core.
*
         CALL PCGEMM( TRANSA, TRANSB, M, N, K, ALPHA, A, IA, JA, DESCA,
     $                B, IB, JB, DESCB, BETA, C, IC, JC, DESCC )
         RETURN
      ENDIF
      IF( NDISK.NE.1 ) THEN
         IF( AONDISK .AND. BONDISK ) THEN
            INFO = -14
         ENDIF
         IF( BONDISK .AND. CONDISK ) THEN
            INFO = -16
         ENDIF
         CALL PXERBLA( DESCA( CTXT_ ), 'PFCGEMM', -INFO )
         RETURN
      ENDIF
      IF( CONDISK ) THEN
         CSIZE = DESCC( SIZE_ )
         CSIZEQUERY = ( CSIZE.EQ.-1 )
         IODEV = DESCC( IODEV_ )
         IF( DESCC( DT_ ).NE.DISK_BLOCK_CYCLIC_2D ) THEN
            INFO = -( 16*100+DT_ )
            CALL PXERBLA( DESCC( CTXT_ ), 'PFCGEMM', 16 )
            RETURN
         ENDIF
         INFO = 0
         CALL LAIO_INFO( DESCC( IODEV_ ), MM, NN, MMB, NNB, MB, NB,
     $                   CSRC, RSRC, ICONTXT )
         IF( MB.NE.DESCC( MB_ ) ) THEN
            INFO = -( 16*100+MB_ )
         ENDIF
         IF( NB.NE.DESCC( NB_ ) ) THEN
            INFO = -( 16*100+NB_ )
         ENDIF
         IF( ICONTXT.NE.DESCC( CTXT_ ) ) THEN
            INFO = -( 16*100+CTXT_ )
         ENDIF
         IF( DESCC( RSRC_ ).NE.RSRC ) THEN
            INFO = -( 16*100+RSRC_ )
         ENDIF
         IF( DESCC( CSRC_ ).NE.CSRC ) THEN
            INFO = -( 16*100+CSRC_ )
         ENDIF
         IF( INFO.NE.0 ) THEN
            CALL PXERBLA( DESCC( CTXT_ ), 'PFCGEMM', 16 )
            RETURN
         ENDIF
         CALL LAIO_INFO( IODEV, MM, NN, MMB, NNB, MB, NB, CSRC, RSRC,
     $                   ICONTXT )
         CALL BLACS_GRIDINFO( ICONTXT, NPROW, NPCOL, MYPROW, MYPCOL )
*
*        Determine dimensions of temporary in-core block
*
         P0 = MYPROW
         Q0 = MYPCOL
         IINC = MIN( M, MMB )
         JINC = MIN( N, NNB )
         LOCP = NUMROC( IINC, MB, MYPROW, P0, NPROW )
         LOCQ = NUMROC( JINC, NB, MYPCOL, Q0, NPCOL )
         CNEED = LOCP*LOCQ
         IF( CSIZE.LT.CNEED ) THEN
            C( 1 ) = CNEED
            INFO = -( 16*100+SIZE_ )
            IF( .NOT.CSIZEQUERY ) THEN
               CALL PXERBLA( DESCC( CTXT_ ), 'PFCGEMM', 16 )
            ENDIF
            RETURN
         ENDIF
         CALL PFMAX2SIZE( CSIZE, M, N, IINC, JINC, IODEV )
         IF( IINC.NE.M ) THEN
            IF( IINC.GE.MMB ) THEN
               IINC = IINC - MOD( IINC, MMB )
            ENDIF
            IF( IINC.GE.MB ) THEN
               IINC = IINC - MOD( IINC, MB )
            ENDIF
         ENDIF
         IF( JINC.NE.N ) THEN
            IF( JINC.GE.NNB ) THEN
               JINC = JINC - MOD( JINC, NNB )
            ENDIF
            IF( JINC.GE.NB ) THEN
               JINC = JINC - MOD( JINC, NB )
            ENDIF
         ENDIF
         ISVALID = ( IINC.GE.MIN( M, MMB ) ) .AND.
     $             ( JINC.GE.MIN( N, NNB ) )
         CALL ASSERT( ISVALID, '** PFxGEMM: internal error, Cneed ',
     $                CNEED )
*
*        Will everything fit?
*        Note the use of p0,q0 to overestimage storage.
*
         P0 = MYPROW
         Q0 = MYPCOL
         LOCP = NUMROC( M, MB, MYPROW, P0, NPROW )
         LOCQ = NUMROC( N, NB, MYPCOL, Q0, NPCOL )
         IERR = 0
         IF( LOCP*LOCQ.GT.CSIZE ) THEN
            IERR = 1
         ENDIF
         CALL IGSUM2D( ICONTXT, 'All', ' ', 1, 1, IERR, 1, -1, -1 )
         ALL_INCORE = ( IERR.EQ.0 )
         IF( ALL_INCORE ) THEN
*
*               Alignment considerations.
*
            IROWC = INDXG2P( IC, MB, MYPROW, RSRC, NPROW )
            ICOLC = INDXG2P( JC, NB, MYPCOL, CSRC, NPCOL )
            P0 = IROWC
            Q0 = ICOLC
            IROWA = INDXG2P( IA, DESCA( MB_ ), MYPROW, DESCA( RSRC_ ),
     $              NPROW )
            ICOLA = INDXG2P( JA, DESCA( NB_ ), MYPCOL, DESCA( CSRC_ ),
     $              NPCOL )
            IROWB = INDXG2P( IB, DESCB( MB_ ), MYPROW, DESCB( RSRC_ ),
     $              NPROW )
            ICOLB = INDXG2P( JB, DESCB( NB_ ), MYPCOL, DESCA( CSRC_ ),
     $              NPCOL )
            IF( NOTRANSA .AND. NOTRANSB ) THEN
               IF( IROWA.NE.IROWC ) THEN
                  Q0 = ICOLB
               ENDIF
            ELSE
               IF( ISTRANSA .AND. NOTRANSB ) THEN
                  IF( IROWA.NE.IROWB ) THEN
                     Q0 = ICOLB
                  ENDIF
               ELSE
                  IF( NOTRANSA .AND. ISTRANSB ) THEN
                     IF( ICOLA.NE.ICOLB ) THEN
                        P0 = IROWA
                     ENDIF
                  ELSE
                     IF( ISTRANSA .AND. ISTRANSB ) THEN
*
*                       No alignment restriction.
*
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
            LOCP = NUMROC( M, MB, MYPROW, P0, NPROW )
            LOCQ = NUMROC( N, NB, MYPCOL, Q0, NPCOL )
            CALL ASSERT( LOCP*LOCQ.LE.CSIZE,
     $                   '** PFxGEMM: internal error, Csize = ', CSIZE )
            LDC = MAX( 1, LOCP )
            CALL DESCINIT( DESCCTMP, M, N, MB, NB, P0, Q0, ICONTXT, LDC,
     $                     INFO )
            CALL ASSERT( INFO.EQ.0,
     $                   '** PFxGEMM: descinit returns info = ', INFO )
            INFO = 0
            CALL CLAREAD( IODEV, M, N, IC, JC, C, 1, 1, DESCCTMP, INFO )
            CALL ASSERT( INFO.EQ.0, '** PFxGEMM: LAREAD returns info = '
     $                   , INFO )
            CALL PCGEMM( TRANSA, TRANSB, M, N, K, ALPHA, A, IA, JA,
     $                   DESCA, B, IB, JB, DESCB, BETA, C, 1, 1,
     $                   DESCCTMP )
            INFO = 0
            CALL CLAWRITE( IODEV, M, N, IC, JC, C, 1, 1, DESCCTMP,
     $                     INFO )
            CALL ASSERT( INFO.EQ.0,
     $                   '** PFxGEMM: LAWRITE returns info = ', INFO )
            C( 1 ) = CNEED
            RETURN
         ENDIF
         JJLO = MAX( 0, INT( JC / JINC ) )
         JJHI = ICEIL( JC+N-1, JINC )
         IILO = MAX( 0, INT( IC / IINC ) )
         IIHI = ICEIL( IA+M-1, IINC )
         DO 30 JJ = JJLO, JJHI
            JSTART = MAX( JC, JJ*JINC+1 )
            JEND = MIN( JC+( N-1 ), ( JJ+1 )*JINC )
            JSIZE = JEND - JSTART + 1
            DO 10 II = IILO, IIHI
               ISTART = MAX( IC, II*IINC+1 )
               IEND = MIN( IC+( M-1 ), ( II+1 )*IINC )
               ISIZE = IEND - ISTART + 1
*
*               Use (ioff,joff) to adjust for alignment
*
               IOFF = 1 + MOD( IINC+( ISTART-1 ), IINC )
               JOFF = 1 + MOD( JINC+( JSTART-1 ), JINC )
               IIC = ISTART
               JJC = JSTART
               MM = ISIZE
               NN = JSIZE
               KK = K
               IF( ISTRANSA ) THEN
                  IIA = IA
                  JJA = JA + ( ISTART-IC )
               ELSE
                  IIA = IA + ( ISTART-IC )
                  JJA = JA
               ENDIF
               IF( ISTRANSB ) THEN
                  IIB = IB + ( JSTART-JC )
                  JJB = JB
               ELSE
                  IIB = IB
                  JJB = JB + ( JSTART-JC )
               ENDIF
               HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 ) .AND. ( KK.GE.1 )
               IF( HASWORK ) THEN
*
*                    Alignment considerations.
*
                  IROWC = INDXG2P( IIC, MB, MYPROW, RSRC, NPROW )
                  ICOLC = INDXG2P( JJC, NB, MYPCOL, CSRC, NPCOL )
                  P0 = IROWC
                  Q0 = ICOLC
                  IROWA = INDXG2P( IIA, DESCA( MB_ ), MYPROW,
     $                    DESCA( RSRC_ ), NPROW )
                  ICOLA = INDXG2P( JJA, DESCA( NB_ ), MYPCOL,
     $                    DESCA( CSRC_ ), NPCOL )
                  IROWB = INDXG2P( IIB, DESCB( MB_ ), MYPROW,
     $                    DESCB( RSRC_ ), NPROW )
                  ICOLB = INDXG2P( JJB, DESCB( NB_ ), MYPCOL,
     $                    DESCB( CSRC_ ), NPCOL )
                  IF( NOTRANSA .AND. NOTRANSB ) THEN
                     IF( IROWC.NE.IROWA ) THEN
                        Q0 = ICOLB
                     ENDIF
                  ELSE
                     IF( ISTRANSA .AND. NOTRANSB ) THEN
                        IF( IROWA.NE.IROWB ) THEN
                           Q0 = ICOLB
                        ENDIF
                     ELSE
                        IF( NOTRANSA .AND. ISTRANSA ) THEN
                           IF( ICOLA.NE.ICOLB ) THEN
                              P0 = IROWA
                           ENDIF
                        ELSE
                           IF( ISTRANSA .AND. ISTRANSB ) THEN
*
*                       No alignment restriction.
*
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
                  LOCP = NUMROC( IINC, MB, MYPROW, P0, NPROW )
                  LOCQ = NUMROC( JINC, NB, MYPCOL, Q0, NPCOL )
                  LDC = MAX( 1, LOCP )
                  CALL ASSERT( LOCP*LOCQ.LE.CSIZE,
     $                         '** PFxGEMM: internal error, Csize = ',
     $                         CSIZE )
                  CALL DESCINIT( DESCCTMP, IINC, JINC, MB, NB, P0, Q0,
     $                           ICONTXT, LDC, INFO )
                  CALL ASSERT( INFO.EQ.0,
     $                         '** PFxGEMM: descinit returns info = ',
     $                         INFO )
                  CALL CLAREAD( IODEV, ISIZE, JSIZE, IIC, JJC, C, IOFF,
     $                          JOFF, DESCCTMP, INFO )
                  CALL ASSERT( INFO.EQ.0,
     $                         '** PFxGEMM: LAREAD returns info ',
     $                         INFO )
                  CALL PCGEMM( TRANSA, TRANSB, MM, NN, KK, ALPHA, A,
     $                         IIA, JJA, DESCA, B, IIB, JJB, DESCB,
     $                         BETA, C, IOFF, JOFF, DESCC )
                  CALL CLAWRITE( IODEV, ISIZE, JSIZE, IIC, JJC, C, IOFF,
     $                           JOFF, DESCCTMP, INFO )
               ENDIF
* end if (haswork)
   10       CONTINUE
   20       CONTINUE
* end do ii
   30    CONTINUE
   40    CONTINUE
* end do jj
         C( 1 ) = CNEED
      ELSE
         IF( BONDISK ) THEN
*
*       Bring in parts of B and perform matrix-matrix product
*
            BSIZE = DESCB( SIZE_ )
            BSIZEQUERY = ( BSIZE.EQ.-1 )
            IODEV = DESCB( IODEV_ )
            IF( DESCB( DT_ ).NE.DISK_BLOCK_CYCLIC_2D ) THEN
               INFO = -( 14*100+DT_ )
               CALL PXERBLA( DESCB( CTXT_ ), 'PFCGEMM', 14 )
               RETURN
            ENDIF
            INFO = 0
            CALL LAIO_INFO( DESCB( IODEV_ ), MM, NN, MMB, NNB, MB, NB,
     $                      CSRC, RSRC, ICONTXT )
            IF( MB.NE.DESCB( MB_ ) ) THEN
               INFO = -( 14*100+MB_ )
            ENDIF
            IF( NB.NE.DESCB( NB_ ) ) THEN
               INFO = -( 14*100+NB_ )
            ENDIF
            IF( ICONTXT.NE.DESCB( CTXT_ ) ) THEN
               INFO = -( 14*100+CTXT_ )
            ENDIF
            IF( DESCB( RSRC_ ).NE.RSRC ) THEN
               INFO = -( 14*100+RSRC_ )
            ENDIF
            IF( DESCB( CSRC_ ).NE.CSRC ) THEN
               INFO = -( 14*100+CSRC_ )
            ENDIF
            IF( INFO.NE.0 ) THEN
               CALL PXERBLA( DESCB( CTXT_ ), 'PFCGEMM', 14 )
               RETURN
            ENDIF
            CALL BLACS_GRIDINFO( ICONTXT, NPROW, NPCOL, MYPROW, MYPCOL )
*
*        Determine storage requirements.
*
            IF( ISTRANSB ) THEN
               M_B = N
               N_B = K
            ELSE
               M_B = K
               N_B = N
            ENDIF
*
*         Allocate buffer.
*         Note the use of p0,q0 to overestimate storage.
*
            P0 = MYPROW
            Q0 = MYPCOL
            LOCP = NUMROC( MIN( MMB, M_B ), MB, MYPROW, P0, NPROW )
            LOCQ = NUMROC( MIN( NNB, N_B ), NB, MYPCOL, Q0, NPCOL )
            BNEED = LOCP*LOCQ
            CALL PFMAX2SIZE( DESCB( SIZE_ ), M_B, N_B, IINC, JINC,
     $                       DESCB( IODEV_ ) )
            ISVALID = ( BNEED.LE.BSIZE ) .AND.
     $                ( IINC.GE.MIN( M_B, MMB ) ) .AND.
     $                ( JINC.GE.MIN( N_B, NNB ) )
            IF( .NOT.ISVALID ) THEN
               B( 1 ) = BNEED
               INFO = -( 14*100+SIZE_ )
               IF( .NOT.BSIZEQUERY ) THEN
                  CALL PXERBLA( ICONTXT, 'PFCGEMM', 14 )
               ENDIF
               RETURN
            ENDIF
*
*        Will everything fit in core?
*        Note the use of p0,q0
*
            IROWB = INDXG2P( IB, DESCB( MB_ ), MYPROW, RSRC, NPROW )
            ICOLB = INDXG2P( JB, DESCB( NB_ ), MYPCOL, CSRC, NPCOL )
            P0 = MYPROW
            Q0 = MYPCOL
            LOCP = NUMROC( M_B, MB, MYPROW, P0, NPROW )
            LOCQ = NUMROC( N_B, NB, MYPCOL, Q0, NPCOL )
            IERR = 0
            IF( LOCP*LOCQ.GT.BSIZE ) THEN
               IERR = 1
            ENDIF
            CALL IGSUM2D( ICONTXT, 'A', ' ', 1, 1, IERR, 1, -1, -1 )
            ALL_INCORE = ( IERR.EQ.0 )
            IF( ALL_INCORE ) THEN
*
*                 Alignment consideration.
*
               IIB = IB
               JJB = JB
               IIA = IA
               JJA = JA
               IIC = IC
               JJC = JC
               IINC = M_B
               JINC = N_B
               IROWB = INDXG2P( IIB, MB, MYPROW, RSRC, NPROW )
               ICOLB = INDXG2P( JJB, NB, MYPCOL, CSRC, NPCOL )
               P0 = IROWB
               Q0 = ICOLB
               IROWA = INDXG2P( IIA, DESCA( MB_ ), MYPROW,
     $                 DESCA( RSRC_ ), NPROW )
               ICOLA = INDXG2P( JJA, DESCA( NB_ ), MYPCOL,
     $                 DESCA( CSRC_ ), NPCOL )
               IROWC = INDXG2P( IIC, DESCC( MB_ ), MYPROW,
     $                 DESCC( RSRC_ ), NPROW )
               ICOLC = INDXG2P( JJC, DESCC( NB_ ), MYPCOL,
     $                 DESCC( CSRC_ ), NPCOL )
               IF( NOTRANSA .AND. NOTRANSB ) THEN
                  IF( IROWA.NE.IROWC ) THEN
                     Q0 = ICOLC
                  ENDIF
               ELSE
                  IF( ISTRANSA .AND. NOTRANSB ) THEN
                     IF( IROWA.NE.IROWB ) THEN
                        Q0 = ICOLC
                     ENDIF
                  ELSE
                     IF( NOTRANSA .AND. ISTRANSB ) THEN
                        IF( IROWA.NE.IROWC ) THEN
                           Q0 = ICOLA
                        ENDIF
                     ELSE
                        IF( ISTRANSA .AND. ISTRANSB ) THEN
*
*                        No processor alignment restriction.
*
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
               LOCP = NUMROC( IINC, MB, MYPROW, P0, NPROW )
               LOCQ = NUMROC( JINC, NB, MYPCOL, Q0, NPCOL )
               LDB = MAX( 1, LOCP )
               CALL ASSERT( LOCP*LOCQ.LE.DESCB( SIZE_ ),
     $                      '** PFxGEMM: internal error Bsize = ',
     $                      BSIZE )
               CALL DESCINIT( DESCBTMP, IINC, JINC, MB, NB, P0, Q0,
     $                        ICONTXT, LDB, INFO )
               CALL ASSERT( INFO.EQ.0,
     $                      '** PFxGEMM: descinit returns info ', INFO )
               INFO = 0
               CALL DESCINIT( DESCBTMP, M_B, N_B, DESCB( MB_ ),
     $                        DESCB( NB_ ), P0, Q0, ICONTXT, LDB, INFO )
               CALL ASSERT( INFO.EQ.0,
     $                      '** PFxGEMM: descinit returns info = ',
     $                      INFO )
*
*                Read everything in and perform in-core computation.
*
               CALL CLAREAD( IODEV, M_B, N_B, IIB, JJB, B, 1, 1,
     $                       DESCBTMP, INFO )
               CALL ASSERT( INFO.EQ.0,
     $                      '** PFxGEMM: LAREAD returns info = ', INFO )
               CALL PCGEMM( TRANSA, TRANSB, M, N, K, ALPHA, A, IIA, JJA,
     $                      DESCA, B, 1, 1, DESCB, BETA, C, IIC, JJC,
     $                      DESCC )
               B( 1 ) = BNEED
               RETURN
            ENDIF
*
*           Prescale C.
*
            IF( BETA.NE.ONE ) THEN
               DO 50 JJ = 1, N
                  CALL PCSCAL( M, BETA, C, IC, ( JC-1 )+JJ, DESCC, 1 )
   50          CONTINUE
   60          CONTINUE
            ENDIF
            LBETA = ONE
            JJLO = MAX( 0, JB / JINC )
            JJHI = ICEIL( JB+N_B-1, JINC )
            IILO = MAX( 0, IB / IINC )
            IIHI = ICEIL( IB+M_B-1, IINC )
            DO 90 JJ = JJLO, JJHI
               JSTART = MAX( JB, JJ*JINC+1 )
               JEND = MIN( JB+( N_B-1 ), ( JJ+1 )*JINC )
               JSIZE = JEND - JSTART + 1
               DO 70 II = IILO, IIHI
                  ISTART = MAX( IB, II*IINC+1 )
                  IEND = MIN( IB+( M_B-1 ), ( II+1 )*IINC )
                  ISIZE = IEND - ISTART + 1
*
*               Use (ioff,joff) for alignment
*
                  IOFF = 1 + MOD( IINC+( ISTART-1 ), IINC )
                  JOFF = 1 + MOD( JINC+( JSTART-1 ), JINC )
                  IIB = ISTART
                  JJB = JSTART
                  IF( ISTRANSB ) THEN
                     MM = M
                     NN = ISIZE
                     KK = JSIZE
                     IIA = IA
                     JJA = JA + ( JJB-JB )
                     IIC = IC
                     JJC = JC + ( IIB-IB )
                  ELSE
                     MM = M
                     NN = JSIZE
                     KK = ISIZE
                     IIA = IA
                     JJA = JA + ( IIB-IB )
                     IIC = JSIZE
                     JJC = JC + ( JJB-JB )
                  ENDIF
                  HASWORK = ( MM.GE.1 ) .AND. ( NN.GE.1 ) .AND.
     $                      ( KK.GE.1 )
                  IF( HASWORK ) THEN
*
*                       Alignment considerations.
*
                     IROWB = INDXG2P( IIB, MB, MYPROW, RSRC, NPROW )
                     ICOLB = INDXG2P( JJB, NB, MYPCOL, CSRC, NPCOL )
                     P0 = IROWB
                     Q0 = ICOLB
                     IROWA = INDXG2P( IIA, DESCA( MB_ ), MYPROW,
     $                       DESCA( RSRC_ ), NPROW )
                     ICOLA = INDXG2P( JJA, DESCA( NB_ ), MYPCOL,
     $                       DESCA( CSRC_ ), NPCOL )
                     IROWC = INDXG2P( IIC, DESCC( MB_ ), MYPROW,
     $                       DESCC( RSRC_ ), NPROW )
                     ICOLC = INDXG2P( JJC, DESCC( NB_ ), MYPCOL,
     $                       DESCC( CSRC_ ), NPCOL )
                     IF( NOTRANSA .AND. NOTRANSB ) THEN
                        IF( IROWA.NE.IROWC ) THEN
                           Q0 = ICOLC
                        ENDIF
                     ELSE
                        IF( ISTRANSA .AND. NOTRANSB ) THEN
                           IF( IROWA.NE.IROWB ) THEN
                              Q0 = ICOLC
                           ENDIF
                        ELSE
                           IF( NOTRANSA .AND. ISTRANSB ) THEN
                              IF( IROWA.NE.IROWC ) THEN
                                 Q0 = ICOLA
                              ENDIF
                           ELSE
                              IF( ISTRANSA .AND. ISTRANSB ) THEN
*
*                       No processor alignment restriction.
*
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDIF
                     LOCP = NUMROC( IINC, MB, MYPROW, P0, NPROW )
                     LOCQ = NUMROC( JINC, NB, MYPCOL, Q0, NPCOL )
                     LDB = MAX( 1, LOCP )
                     CALL ASSERT( LOCP*LOCQ.LE.DESCB( SIZE_ ),
     $                            '** PFxGEMM: internal error Bsize = ',
     $                            BSIZE )
                     CALL DESCINIT( DESCBTMP, IINC, JINC, MB, NB, P0,
     $                              Q0, ICONTXT, LDB, INFO )
                     CALL ASSERT( INFO.EQ.0,
     $                            '** PFxGEMM: descinit returns info ',
     $                            INFO )
                     CALL CLAREAD( IODEV, ISIZE, JSIZE, IIB, JJB, B,
     $                             IOFF, JOFF, DESCBTMP, INFO )
                     CALL ASSERT( INFO.EQ.0,
     $                            '** PFxGEMM: LAREAD returns info ',
     $                            INFO )
                     IIA = IA
                     JJA = JA
                     CALL PCGEMM( TRANSA, TRANSB, MM, NN, KK, ALPHA, A,
     $                            IIA, JJA, DESCA, B, IOFF, JOFF,
     $                            DESCBTMP, LBETA, C, IIC, JJC, DESCC )
                  ENDIF
* end if (haswork)
   70          CONTINUE
   80          CONTINUE
* end do ii
   90       CONTINUE
  100       CONTINUE
* end do jj
            B( 1 ) = BNEED
         ELSE
            IF( AONDISK ) THEN
               ASIZE = DESCA( SIZE_ )
               ASIZEQUERY = ( ASIZE.EQ.-1 )
               IODEV = DESCA( IODEV_ )
               IF( DESCA( DT_ ).NE.DISK_BLOCK_CYCLIC_2D ) THEN
                  INFO = -( 10*100+DT_ )
                  CALL PXERBLA( DESCA( CTXT_ ), 'PFCGEMM', 10 )
                  RETURN
               ENDIF
               INFO = 0
               CALL LAIO_INFO( DESCA( IODEV_ ), MM, NN, MMB, NNB, MB,
     $                         NB, CSRC, RSRC, ICONTXT )
               IF( MB.NE.DESCA( MB_ ) ) THEN
                  INFO = -( 10*100+MB_ )
               ENDIF
               IF( NB.NE.DESCA( NB_ ) ) THEN
                  INFO = -( 10*100+NB_ )
               ENDIF
               IF( ICONTXT.NE.DESCA( CTXT_ ) ) THEN
                  INFO = -( 10*100+CTXT_ )
               ENDIF
               IF( DESCA( RSRC_ ).NE.RSRC ) THEN
                  INFO = -( 10*100+RSRC_ )
               ENDIF
               IF( DESCA( CSRC_ ).NE.CSRC ) THEN
                  INFO = -( 10*100+CSRC_ )
               ENDIF
               IF( INFO.NE.0 ) THEN
                  CALL PXERBLA( DESCA( CTXT_ ), 'PFCGEMM', 10 )
                  RETURN
               ENDIF
               CALL BLACS_GRIDINFO( ICONTXT, NPROW, NPCOL, MYPROW,
     $                              MYPCOL )
               IF( ISTRANSA ) THEN
                  M_A = K
                  N_A = M
               ELSE
                  M_A = M
                  N_A = K
               ENDIF
*
*                Allocate temporary buffers.
*                Note the use of p0,q0 to overestimate storage.
*
               P0 = MYPROW
               Q0 = MYPCOL
               IINC = MIN( M_A, MMB )
               JINC = MIN( N_A, NNB )
               LOCP = NUMROC( IINC, MB, MYPROW, P0, NPROW )
               LOCQ = NUMROC( JINC, NB, MYPCOL, Q0, NPCOL )
               ANEED = LOCP*LOCQ
               CALL PFMAX2SIZE( DESCA( SIZE_ ), M_A, N_A, IINC, JINC,
     $                          DESCA( IODEV_ ) )
               ISVALID = ( ASIZE.GE.ANEED ) .AND.
     $                   ( IINC.GE.MIN( M_A, MMB ) ) .AND.
     $                   ( JINC.GE.MIN( N_A, NNB ) )
               IF( .NOT.ISVALID ) THEN
*
*                Insufficient temporary storage
*                return min requirement in A(1)
*
*
                  A( 1 ) = ANEED
                  INFO = -( 10*100+SIZE_ )
                  IF( .NOT.ASIZEQUERY ) THEN
                     CALL PXERBLA( ICONTXT, 'PFCGEMM', 10 )
                  ENDIF
                  RETURN
               ENDIF
*
*         Will everything fit?
*         Note the use of p0,q0 to over estimate storage.
*
               IROWA = INDXG2P( IA, MB, MYPROW, RSRC, NPROW )
               ICOLA = INDXG2P( JA, NB, MYPCOL, CSRC, NPCOL )
               P0 = MYPROW
               Q0 = MYPCOL
               LOCP = NUMROC( M_A, MB, MYPROW, P0, NPROW )
               LOCQ = NUMROC( N_A, NB, MYPCOL, Q0, NPCOL )
               IERR = 0
               IF( LOCP*LOCQ.GT.ASIZE ) THEN
                  IERR = 1
               ENDIF
               CALL IGSUM2D( ICONTXT, 'All', ' ', 1, 1, IERR, 1, -1,
     $                       -1 )
               ALL_INCORE = ( IERR.EQ.0 )
               IF( ALL_INCORE ) THEN
*
*                Perform in-core operation.
*
*
*                Alignment consideration.
*
                  P0 = IROWA
                  Q0 = ICOLA
                  IROWB = INDXG2P( IB, DESCB( MB_ ), MYPROW,
     $                    DESCB( RSRC_ ), NPROW )
                  ICOLB = INDXG2P( JB, DESCC( NB_ ), MYPCOL,
     $                    DESCB( CSRC_ ), NPCOL )
                  IROWC = INDXG2P( IC, DESCC( MB_ ), MYPROW,
     $                    DESCC( RSRC_ ), NPROW )
                  ICOLC = INDXG2P( JC, DESCC( NB_ ), MYPCOL,
     $                    DESCC( CSRC_ ), NPCOL )
                  IF( NOTRANSA .AND. NOTRANSB ) THEN
                     IF( ICOLB.NE.ICOLC ) THEN
                        P0 = IROWC
                     ENDIF
                  ELSE
                     IF( ( .NOT.NOTRANSA ) .AND. NOTRANSB ) THEN
                        IF( ICOLB.NE.ICOLC ) THEN
                           P0 = IROWB
                        ENDIF
                     ELSE
                        IF( NOTRANSA .AND. ( .NOT.NOTRANSB ) ) THEN
                           IF( ICOLA.NE.ICOLB ) THEN
                              P0 = IROWC
                           ENDIF
                        ELSE
                           IF( ( .NOT.NOTRANSA ) .AND.
     $                         ( .NOT.NOTRANSB ) ) THEN
*
*                        No restrictions.
*
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
                  LOCP = NUMROC( M_A, MB, MYPROW, P0, NPROW )
                  LOCQ = NUMROC( N_A, NB, MYPCOL, Q0, NPCOL )
                  CALL ASSERT( LOCP*LOCQ.LE.ASIZE,
     $                         '** PFxGEMM: internal error, Asize = ',
     $                         ASIZE )
                  LDA = MAX( 1, LOCP )
                  INFO = 0
                  CALL DESCINIT( DESCATMP, M_A, N_A, MB, NB, P0, Q0,
     $                           ICONTXT, LDA, INFO )
                  CALL ASSERT( INFO.EQ.0,
     $                         '** PFxGEMM: descinit returns info = ',
     $                         INFO )
                  CALL CLAREAD( IODEV, M_A, N_A, IA, JA, A, 1, 1,
     $                          DESCATMP, INFO )
                  CALL ASSERT( INFO.EQ.0,
     $                         '** PFxGEMM: LAREAD returns info = ',
     $                         INFO )
                  CALL PCGEMM( TRANSA, TRANSB, M, N, K, ALPHA, A, 1, 1,
     $                         DESCATMP, B, IB, JB, DESCB, BETA, C, IC,
     $                         JC, DESCC )
                  A( 1 ) = ANEED
                  RETURN
               ENDIF
* end if (all_incore)
*
*         Prescale C.
*
               IF( BETA.NE.ONE ) THEN
                  DO 110 JJ = 1, N
                     CALL PCSCAL( M, BETA, C, IC, ( JC-1 )+JJ, DESCC,
     $                            1 )
  110             CONTINUE
  120             CONTINUE
               ENDIF
               LBETA = ONE
               JJLO = MAX( 0, INT( JA / JINC ) )
               JJHI = ICEIL( JA+N_A-1, JINC )
               IILO = MAX( 0, INT( IA / IINC ) )
               IIHI = ICEIL( IA+M_A-1, IINC )
               DO 150 JJ = JJLO, JJHI
                  JSTART = MAX( JA, JJ*JINC+1 )
                  JEND = MIN( JA+( N_A-1 ), ( JJ+1 )*JINC )
                  JSIZE = JEND - JSTART + 1
                  DO 130 II = IILO, IIHI
                     ISTART = MAX( IA, II*IINC+1 )
                     IEND = MIN( IA+( M_A-1 ), ( II+1 )*IINC )
                     ISIZE = IEND - ISTART + 1
                     CALL ASSERT( IINC.GE.1, '** pfgemm: iinc < 1 ',
     $                            IINC )
                     CALL ASSERT( JINC.GE.1, '** pfgemm: jinc < 1 ',
     $                            JINC )
                     IOFF = 1 + MOD( IINC+( ISTART-1 ), IINC )
                     JOFF = 1 + MOD( JINC+( JSTART-1 ), JINC )
                     IIA = ISTART
                     JJA = JSTART
                     IF( NOTRANSA ) THEN
                        MM = ISIZE
                        NN = N
                        KK = JSIZE
                        IIC = IC + ( IIA-IA )
                        JJC = JC
                        IIB = IB + ( JJA-JA )
                        JJB = JB
                     ELSE
                        MM = JSIZE
                        NN = N
                        KK = ISIZE
                        IIC = IC + ( JJA-JA )
                        JJC = JC
                        IIB = IB + ( IIA-IA )
                        JJB = JB
                     ENDIF
                     HASWORK = ( ( MM.GE.1 ) .AND. ( NN.GE.1 ) .AND.
     $                         ( KK.GE.1 ) )
                     IF( HASWORK ) THEN
*
*                       Alignment consideration.
*
                        IROWA = INDXG2P( IIA, MB, MYPROW, RSRC, NPROW )
                        ICOLA = INDXG2P( JJA, NB, MYPCOL, CSRC, NPCOL )
                        P0 = IROWA
                        Q0 = ICOLA
                        IROWB = INDXG2P( IIB, DESCB( MB_ ), MYPROW,
     $                          DESCB( RSRC_ ), NPROW )
                        ICOLB = INDXG2P( JJB, DESCB( NB_ ), MYPCOL,
     $                          DESCB( CSRC_ ), NPCOL )
                        IROWC = INDXG2P( IIC, DESCC( MB_ ), MYPROW,
     $                          DESCC( RSRC_ ), NPROW )
                        ICOLC = INDXG2P( JJC, DESCC( NB_ ), MYPCOL,
     $                          DESCC( CSRC_ ), NPCOL )
                        IF( NOTRANSA .AND. NOTRANSB ) THEN
                           IF( ICOLB.NE.ICOLC ) THEN
                              P0 = IROWC
                           ENDIF
                        ELSE
                           IF( ISTRANSA .AND. NOTRANSB ) THEN
                              IF( ICOLB.NE.ICOLC ) THEN
                                 P0 = IROWB
                              ENDIF
                           ELSE
                              IF( NOTRANSA .AND. ISTRANSB ) THEN
                                 IF( ICOLA.NE.ICOLB ) THEN
                                    P0 = IROWC
                                 ENDIF
                              ELSE
                                 IF( ISTRANSA .AND. ISTRANSB ) THEN
*
*                        No alignment restrictions.
*
                                 ENDIF
                              ENDIF
                           ENDIF
                        ENDIF
                        LOCP = NUMROC( IINC, MB, MYPROW, P0, NPROW )
                        LOCQ = NUMROC( JINC, NB, MYPCOL, Q0, NPCOL )
                        CALL ASSERT( LOCP*LOCQ.LE.ASIZE,
     $                             '** PFxGEMM: internal error Asize = '
     $                               , ASIZE )
                        LDA = MAX( 1, LOCP )
                        CALL DESCINIT( DESCATMP, IINC, JINC, MB, NB, P0,
     $                                 Q0, ICONTXT, LDA, INFO )
                        CALL ASSERT( INFO.EQ.0,
     $                            '** PFxGEMM: descinit returns info = '
     $                               , INFO )
                        CALL CLAREAD( IODEV, ISIZE, JSIZE, IIA, JJA, A,
     $                                IOFF, JOFF, DESCATMP, INFO )
                        CALL ASSERT( INFO.EQ.0,
     $                               '** PFxGEMM: LAREAD returns info ',
     $                               INFO )
                        CALL PCGEMM( TRANSA, TRANSB, MM, NN, KK, ALPHA,
     $                               A, IOFF, JOFF, DESCATMP, B, IIB,
     $                               JJB, DESCB, LBETA, C, IIC, JJC,
     $                               DESCC )
                     ENDIF
* end if (haswork)
  130             CONTINUE
  140             CONTINUE
* end do ii
  150          CONTINUE
  160          CONTINUE
* end do jj
               A( 1 ) = ANEED
            ENDIF
         ENDIF
      ENDIF
      INFO = 0
      RETURN
      END
