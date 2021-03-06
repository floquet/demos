      SUBROUTINE PFZQRSOLVE( IODEV, M, N, NRHS, IC, JC, DESCC, TAUC, A,
     $                       DESCA, ASIZE, B, IB, JB, DESCB, WORK,
     $                       LWORK )
*
* Purpose:
* ========
*
* PFxQRSOLVE is an Internal auxiliary routine that
* uses QR factors generated by PFxGEQRF
* to solve linear equations.
*
* Matrix B gets overwritten with the solution
* sub(C) matrix  is m by n, B is m by nrhs
*
*
*     .. Parameters ..
      INTEGER            DLEN_, DT_
      PARAMETER          ( DLEN_ = 9, DT_ = 1 )
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
      INTEGER            ASIZE, IB, IC, IODEV, JB, JC, LWORK, M, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( DLEN_ ), DESCB( DLEN_ ), DESCC( DLEN_ )
      COMPLEX*16         A( * ), B( * ), TAUC( * ), WORK( * )
*     ..
*     .. Local Scalars ..
      LOGICAL            ASIZEQUERY, ISVALID, LWORKQUERY
      INTEGER            ANEED, IFREE, INEED, INFO, KK, LOCP, LOCQ, MM,
     $                   MYID, MYPCOL, MYPROW, NFREE, NN, NPCOL, NPROC,
     $                   NPROW, P0, Q0
      COMPLEX*16         ALPHA, ONE
*     ..
*     .. Local Arrays ..
      INTEGER            DESC( FLEN_ )
*     ..
*     .. External Functions ..
      INTEGER            NUMROC
      EXTERNAL           NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           ASSERT, BLACS_GRIDINFO, BLACS_PINFO, DESCINIT,
     $                   PFZLATRSM, PFZUNMQR, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, MAX, MIN
*     ..
*     .. Executable Statements ..
      LWORKQUERY = ( LWORK.EQ.-1 )
      ONE = DCMPLX( DBLE( 1 ) )
* check parameters
      CALL BLACS_PINFO( MYID, NPROC )
      ISVALID = ( 1.LE.IC ) .AND. ( 1.LE.JC ) .AND.
     $          ( IC+M-1.LE.DESCC( M_ ) ) .AND.
     $          ( JC+N-1.LE.DESCC( N_ ) ) .AND. ( 1.LE.IB ) .AND.
     $          ( 1.LE.JB ) .AND. ( IB+M-1.LE.DESCB( M_ ) ) .AND.
     $          ( JB+NRHS-1.LE.DESCB( N_ ) ) .AND.
     $          ( M.LE.DESCA( M_ ) ) .AND. ( 1.LE.DESCA( N_ ) )
      CALL ASSERT( ISVALID, ' ** PFQRSOLVE: invalid parameters, myid= ',
     $             MYID )
      IFREE = 1
      NFREE = LWORK
* perform Householder reduction
      MM = M
      NN = NRHS
      KK = MIN( M, N )
      INEED = MAX( 1, DESCA( N_ ) )
      CALL DESCINIT( DESC, DESCC( M_ ), DESCC( N_ ), DESCC( MB_ ),
     $               DESCC( NB_ ), DESCC( RSRC_ ), DESCC( CSRC_ ),
     $               DESCC( CTXT_ ), DESCC( LLD_ ), INFO )
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYPROW,
     $                     MYPCOL )
      P0 = MYPROW
      Q0 = MYPCOL
      LOCP = NUMROC( DESCA( M_ ), DESCA( MB_ ), MYPROW, P0, NPROW )
      LOCQ = NUMROC( DESCA( N_ ), DESCA( NB_ ), MYPCOL, Q0, NPCOL )
      ANEED = LOCP*LOCQ
      ASIZEQUERY = ( ASIZE.EQ.-1 )
      IF( ANEED.GT.ASIZE ) THEN
         A( 1 ) = ANEED
         WORK( 1 ) = INEED
         INFO = -11
         IF( .NOT.ASIZEQUERY ) THEN
            CALL PXERBLA( DESCA( CTXT_ ), 'PFxQRSOLVE', 11 )
         ENDIF
         RETURN
      ENDIF
      IF( INEED.GT.LWORK ) THEN
         WORK( 1 ) = INEED
         A( 1 ) = ANEED
         INFO = -17
         IF( .NOT.LWORKQUERY ) THEN
            CALL PXERBLA( DESCA( CTXT_ ), 'PFxQRSOLVE', 17 )
         ENDIF
      ENDIF
      DESC( DT_ ) = DISK_BLOCK_CYCLIC_2D
      DESC( IODEV_ ) = IODEV
      DESC( SIZE_ ) = ASIZE
      CALL PFZUNMQR( 'Left', 'C', MM, NN, KK, A, IC, JC, DESC, TAUC, B,
     $               IB, JB, DESCB, WORK( IFREE ), NFREE, INFO )
      CALL ASSERT( INFO.EQ.0, '** PFQRSOLVE: PFORMQR returns info ',
     $             INFO )
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYPROW,
     $                     MYPCOL )
      LOCP = NUMROC( DESCA( M_ ), DESCA( MB_ ), MYPROW, DESCA( RSRC_ ),
     $       NPROW )
      LOCQ = NUMROC( DESCA( N_ ), DESCA( NB_ ), MYPCOL, DESCA( CSRC_ ),
     $       NPCOL )
      ANEED = LOCP*LOCQ
      IF( ANEED.GT.ASIZE ) THEN
         A( 1 ) = ANEED
         WORK( 1 ) = INEED
         INFO = -11
         IF( .NOT.ASIZEQUERY ) THEN
            CALL PXERBLA( DESCA( CTXT_ ), 'PFxQRSOLVE', 11 )
         ENDIF
         RETURN
      ENDIF
      MM = M
      NN = NRHS
      ALPHA = ONE
      CALL PFZLATRSM( 'L', 'U', 'N', 'N', IODEV, MM, NN, IC, JC, ALPHA,
     $                B, IB, JB, DESCB, A, ASIZE, INFO )
      RETURN
      END
