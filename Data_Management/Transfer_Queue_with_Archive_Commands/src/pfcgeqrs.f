      SUBROUTINE PFCGEQRS( M, N, NRHS, A, IA, JA, DESCA, TAU, X, IX, JX,
     $                     DESCX, WORK, LWORK, INFO )
*
*
*  -- ScaLAPACK routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     Oct 10, 1996
*
*
*
*  Purpose:
*  =======
*
*  PFxGEQRS solves a system of distributed linear equations
*
*       sub(A)*X = sub(B)
*
*  where sub(A) is factored using the QR factorization
*  computed by PFxGEQRF.
*
*  sub(A) denotes A( ia:(ia+m-1),  ja:(ja+n-1) )
*
*  Note:
*  =====
*
*  work(Lwork) is used to store intermediate 'tau' vector and
*  also required for PxUNMQR or PxORMQR.
*  Since a variable width panel is used, the amount of storage
*  required in work(lwork) is dependent on the widest
*  possible panel.
*
*  The value of work(1) on exit of this routine is the
*  amount of lwork actually used.
*
*  A very conservative estimate is
*
*  lwork = tau_need + ormqr_need
*
*  p0 = myrow, q0= mycol
*  tau_need = numroc( n, nb_A, mypcol, p0, npcol)
*
*  ormqr_need = nb_A*nb_A + nb_A*max( (nb_A-1)/2, (NqB0 + MpB0))
*
*  MpB0 = iceil(numroc(m, nb_X,myprow,p0,nprow),nb_X)*nb_X
*  NpB0 = iceil(numroc(n, nb_X,mypcol,q0,nprow),mb_X)*mb_X
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
*     ..
*     .. Scalar Arguments ..
      INTEGER            IA, INFO, IX, JA, JX, LWORK, M, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( FLEN_ ), DESCX( DLEN_ )
      COMPLEX            A( * ), TAU( * ), WORK( * ), X( * )
*     ..
*     .. Local Scalars ..
      LOGICAL            ASIZEQUERY, LWORKQUERY
      INTEGER            ANEED, ASIZE, CSRC, ICONTXT, ICTXT, IDUMMY,
     $                   INEED, IODEV, LDD, LOCP, LOCQ, LWORK0, LWORK1,
     $                   MB, MM, MMB, MP0, MYPCOL, MYPROW, NB, NCOL, NN,
     $                   NNB, NPCOL, NPROW, NQ0, P0, Q0, RSRC, SAVEDT
*     ..
*     .. Local Arrays ..
      INTEGER            DESC1( DLEN_ ), DESC2( DLEN_ )
*     ..
*     .. External Functions ..
      INTEGER            NUMROC
      EXTERNAL           NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           ASSERT, BLACS_GRIDINFO, CHK1MAT, DESCINIT,
     $                   LAIO_INFO, PCHK1MAT, PFCQRSOLVE, PFMAXSIZE,
     $                   PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
      INFO = 0
      ICTXT = DESCA( CTXT_ )
      ICTXT = ICTXT
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYPROW, MYPCOL )
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 7*100+CTXT_ )
         RETURN
      ENDIF
      SAVEDT = DESCA( DT_ )
      DESCA( DT_ ) = BLOCK_CYCLIC_2D
      CALL CHK1MAT( M, 1, N, 2, IA, JA, DESCA, 7, INFO )
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PFCGEQRS', 7 )
         RETURN
      ENDIF
      IDUMMY = 1
      CALL PCHK1MAT( M, 1, N, 2, IA, JA, DESCA, 7, IDUMMY, IDUMMY,
     $               IDUMMY, INFO )
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PFCGEQRS', 7 )
         RETURN
      ENDIF
      DESCA( DT_ ) = SAVEDT
      IF( DESCA( DT_ ).NE.DISK_BLOCK_CYCLIC_2D ) THEN
         INFO = -( 7*100+DT_ )
         CALL PXERBLA( DESCA( CTXT_ ), 'PFCGEQRS', 7 )
         RETURN
      ENDIF
      INFO = 0
      CALL LAIO_INFO( DESCA( IODEV_ ), MM, NN, MMB, NNB, MB, NB, CSRC,
     $                RSRC, ICONTXT )
      IF( MB.NE.DESCA( MB_ ) ) THEN
         INFO = -( 7*100+MB_ )
      ENDIF
      IF( NB.NE.DESCA( NB_ ) ) THEN
         INFO = -( 7*100+NB_ )
      ENDIF
      IF( ICONTXT.NE.DESCA( CTXT_ ) ) THEN
         INFO = -( 7*100+CTXT_ )
      ENDIF
      IF( DESCA( RSRC_ ).NE.RSRC ) THEN
         INFO = -( 7*100+RSRC_ )
      ENDIF
      IF( DESCA( CSRC_ ).NE.CSRC ) THEN
         INFO = -( 7*100+CSRC_ )
      ENDIF
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PFCGEQRS', 7 )
         RETURN
      ENDIF
      IODEV = DESCA( IODEV_ )
      ASIZE = DESCA( SIZE_ )
      ASIZEQUERY = ( ASIZE.EQ.-1 )
      LWORKQUERY = ( LWORK.EQ.-1 )
*
*   Note the use of p0,q0 to overestimate storage.
*
      P0 = MYPROW
      Q0 = MYPCOL
      LOCP = NUMROC( M, MB, MYPROW, P0, NPROW )
      LOCQ = NUMROC( MIN( N, NNB ), NB, MYPCOL, Q0, NPCOL )
      ANEED = LOCP*LOCQ
*
*   Over-estimate amount of lwork storage.
*   lwork0 storage for PxORMQR().
*
      MB = DESCA( MB_ )
      NB = DESCA( NB_ )
      MP0 = MAX( 1, NUMROC( M+MB, MB, MYPROW, P0, NPROW ) ) + MB
      NQ0 = MAX( 1, NUMROC( N+NB, NB, MYPCOL, Q0, NPCOL ) ) + NB
      LWORK0 = MAX( NB*NB, ( NQ0+MP0 )*NB ) + NB*NB
      LWORK1 = MAX( 1, NUMROC( N, NB, MYPCOL, Q0, NPCOL ) )
      INEED = LWORK0 + 2*LWORK1
      IF( ANEED.GT.ASIZE ) THEN
         A( 1 ) = ANEED
         WORK( 1 ) = INEED
         INFO = -( 7*100+SIZE_ )
         IF( .NOT.ASIZEQUERY ) THEN
            CALL PXERBLA( ICTXT, 'PFCGEQRS', 7 )
         ENDIF
         RETURN
      ENDIF
      IF( LWORKQUERY ) THEN
         A( 1 ) = ANEED
         WORK( 1 ) = INEED
         INFO = -14
         RETURN
      ENDIF
      NCOL = -1
      CALL PFMAXSIZE( 'C', ASIZE, M, NCOL, MB, NB, RSRC, CSRC, ICTXT,
     $                INFO )
      NCOL = MAX( 1, MIN( N, NCOL ) )
      LOCP = NUMROC( M, MB, MYPROW, RSRC, NPROW )
      LOCQ = NUMROC( NCOL, NB, MYPCOL, CSRC, NPCOL )
      LDD = MAX( 1, LOCP )
      CALL DESCINIT( DESC1, M, NCOL, MB, NB, RSRC, CSRC, ICTXT, LDD,
     $               INFO )
      IF( INFO.NE.0 ) THEN
         WRITE( *, FMT = 9999 )MYPROW, MYPCOL, M, NCOL, MB, NB, RSRC,
     $      CSRC, ICTXT, LDD, INFO
 9999    FORMAT( '** PFGEQRS: ', ' myprow,mypcol ', 2( 1X, I4 ),
     $         ' m,ncol ', 2( 1X, I6 ), ' rsrc,csrc ', 2( 1X, I4 ),
     $         ' ictxt,ldd,info ', 3( 1X, I6 ) )
      ENDIF
      CALL ASSERT( INFO.EQ.0, '** PFGEQRS: descinit returns info ',
     $             INFO )
      CALL DESCINIT( DESC2, DESCA( M_ ), DESCA( N_ ), DESCA( MB_ ),
     $               DESCA( NB_ ), DESCA( RSRC_ ), DESCA( CSRC_ ),
     $               DESCA( CTXT_ ), DESCA( LLD_ ), INFO )
      CALL ASSERT( INFO.EQ.0, '** PFGEQRS: descinit returns info ',
     $             INFO )
      CALL PFCQRSOLVE( IODEV, M, N, NRHS, IA, JA, DESC2, TAU, A, DESC1,
     $                 ASIZE, X, IX, JX, DESCX, WORK, LWORK )
      INEED = WORK( 1 )
      A( 1 ) = ANEED
      IF( INEED.GT.LWORK ) THEN
         INFO = -14
         IF( .NOT.LWORKQUERY ) THEN
            CALL PXERBLA( ICTXT, 'PFCGEQRS', 14 )
         ENDIF
         RETURN
      ENDIF
      INFO = 0
      RETURN
      END
