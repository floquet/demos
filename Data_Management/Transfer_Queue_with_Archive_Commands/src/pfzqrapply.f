      SUBROUTINE PFZQRAPPLY( IODEV, M, N, K, IC, JC, DESCC, TAUC, A,
     $                       DESCA, TAUA, TAUA_DIM, IB, JB, B, DESCB,
     $                       WORK, LWORK )
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
* Apply Householder transformations to sub(B), where
* sub(B) is B( ib:(ib+m-1), jb:(jb+n-1) )
*
* The Householder transformations are computed by PFGEQRF
* and store in an out-of-core array
*
*
* ===============================================
* sub(B) is m by n,   apply k householder vectors
* ===============================================
*     .. Parameters ..
      INTEGER            DLEN_
      PARAMETER          ( DLEN_ = 9 )
      INTEGER            M_, N_, NB_
      PARAMETER          ( M_ = 3, N_ = 4, NB_ = 6 )
*     ..
*     .. Scalar Arguments ..
      INTEGER            IB, IC, IODEV, JB, JC, K, LWORK, M, N, TAUA_DIM
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( DLEN_ ), DESCB( DLEN_ ), DESCC( DLEN_ )
      COMPLEX*16         A( * ), B( * ), TAUA( TAUA_DIM ), TAUC( * ),
     $                   WORK( LWORK )
*     ..
*     .. Local Scalars ..
      LOGICAL            HASWORK, ISVALID
      CHARACTER          SIDE, TRANS
      INTEGER            GTAU_DIM, IA, IBPOS, ICPOS, IENDA, IFREE, INFO,
     $                   IP_GTAU, ISIZEA, ISTARTA, JA, JBPOS, JCPOS,
     $                   JENDA, JINC, JSIZEA, JSIZEB, JSTARTA, KK, MM,
     $                   MYID, NFREE, NN, NPROC
      COMPLEX*16         DZERO
*     ..
*     .. External Subroutines ..
      EXTERNAL           ASSERT, BLACS_PINFO, PFZCOPYTAU, PZUNMQR,
     $                   ZFILL, ZLAREAD
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
      DZERO = DCMPLX( DBLE( 0 ) )
* check parameter list
      CALL BLACS_PINFO( MYID, NPROC )
      ISVALID = ( 1.LE.IC ) .AND. ( 1.LE.JC ) .AND. ( 1.LE.IB ) .AND.
     $          ( 1.LE.JB ) .AND. ( M.GE.1 ) .AND. ( N.GE.1 ) .AND.
     $          ( IB+M-1.LE.DESCB( M_ ) ) .AND.
     $          ( JB+N-1.LE.DESCB( N_ ) ) .AND.
     $          ( M.LE.DESCA( M_ ) ) .AND. ( DESCA( N_ ).GE.1 )
      CALL ASSERT( ISVALID, '** qrapply: invalid parameters, myid = ',
     $             MYID )
* ============================
* allocate storage for gtau(*)
* ============================
      IFREE = 1
      GTAU_DIM = MAX( 1, MAX( DESCA( N_ ), DESCB( N_ ) ) )
      IP_GTAU = IFREE
      IFREE = IFREE + GTAU_DIM
      NFREE = LWORK - IFREE + 1
      CALL ASSERT( NFREE.GE.0, '** qrfact: increase lwork by ',
     $             ABS( NFREE )+2 )
      JSIZEB = N
      IF( DESCA( N_ ).GE.K ) THEN
         JINC = K
      ELSE
         JINC = DESCA( N_ )
         IF( JINC.GE.DESCA( NB_ ) ) THEN
            JINC = JINC - MOD( JINC, DESCA( NB_ ) )
         ENDIF
      ENDIF
      DO 10 JSTARTA = 1, K, JINC
         JENDA = MIN( K, JSTARTA+JINC-1 )
         JSIZEA = JENDA - JSTARTA + 1
         ISTARTA = JSTARTA
         IENDA = M
         ISIZEA = IENDA - ISTARTA + 1
* read in part of C matrix
         IA = ISTARTA
* not sure whether it should be ia = 1
         JA = 1
         ICPOS = ( IC-1 ) + ISTARTA
         JCPOS = ( JC-1 ) + JSTARTA
         INFO = 0
         CALL ZLAREAD( IODEV, ISIZEA, JSIZEA, ICPOS, JCPOS, A, IA, JA,
     $                 DESCA, INFO )
         CALL ASSERT( INFO.EQ.0,
     $                '** qrapply: LAREAD into A returns info ', INFO )
* copy tau
         CALL ASSERT( ( GTAU_DIM.GE.JSIZEA ),
     $                '** qrapply: increase gtau_dim to ', JSIZEA+1 )
         CALL ZFILL( TAUA_DIM, DZERO, TAUA( 1 ), 1 )
         JA = 1
         CALL PFZCOPYTAU( JSIZEA, JCPOS, TAUC, DESCC, JA, TAUA( 1 ),
     $                    DESCA )
*
*       Apply householder transformations
*       to only a part of matrix B
*
         HASWORK = ( ISIZEA.GE.1 ) .AND. ( JSIZEB.GE.1 ) .AND.
     $             ( JSIZEA.GE.1 )
         IF( HASWORK ) THEN
            IBPOS = ( IB-1 ) + IA
            JBPOS = ( JB-1 ) + 1
            JA = 1
            MM = ISIZEA
            NN = JSIZEB
            KK = JSIZEA
            INFO = 0
            SIDE = 'L'
            TRANS = 'C'
            CALL PZUNMQR( SIDE, TRANS, MM, NN, KK, A, IA, JA, DESCA,
     $                    TAUA( 1 ), B, IBPOS, JBPOS, DESCB,
     $                    WORK( IFREE ), NFREE, INFO )
            CALL ASSERT( INFO.EQ.0, '** qrapply: PORMQR returns info ',
     $                   INFO )
         ENDIF
   10 CONTINUE
   20 CONTINUE
* end do jstartA
      RETURN
      END
