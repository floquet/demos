      SUBROUTINE CALCOLSIZE( NFREE, GNROWA, DESCA, GNCOLA )
* give a fixed amount of storage,
* nfree entries on each processor
* figure out the largest scalapack array that can fit
* the number of rows (gnrowA) is fixed
* intent( out ) :: gncolA
* HUGE is 2^30
* if k = divup(n,m),
* then (k-1)*m < n <= k*m
* it can be computed as
*  if (mod(n,m) == 0) { k = int(n/m); }
*  else { k = int(n/m) + 1; }
* inline statement function
*     .. Parameters ..
      INTEGER            DLEN_
      PARAMETER          ( DLEN_ = 9 )
      INTEGER            CTXT_, MB_, NB_
      PARAMETER          ( CTXT_ = 2, MB_ = 5, NB_ = 6 )
      INTEGER            RSRC_, CSRC_
      PARAMETER          ( RSRC_ = 7, CSRC_ = 8 )
      INTEGER            HUGE
      PARAMETER          ( HUGE = 1073741824 )
*     ..
*     .. Scalar Arguments ..
      INTEGER            GNCOLA, GNROWA, NFREE
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( DLEN_ )
*     ..
*     .. Local Scalars ..
      LOGICAL            ISVALID
      CHARACTER          SCOPE, TOP
      INTEGER            CA, CDEST, CONTXT, LDA, M, MYID, MYPCOL,
     $                   MYPROW, N, NCOL, NPCOL, NPROC, NPROW, NROW, RA,
     $                   RCFLAG, RDEST
*     ..
*     .. External Functions ..
      INTEGER            NUMROC, NUMROCINV
      EXTERNAL           NUMROC, NUMROCINV
*     ..
*     .. External Subroutines ..
      EXTERNAL           ASSERT, BLACS_GRIDINFO, BLACS_PINFO, IGAMN2D
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT, MAX, MIN, MOD
*     ..
*     .. Statement Functions ..
      INTEGER            DIVUP
*     ..
*     .. Statement Function definitions ..
      DIVUP( N, M ) = ( INT( ( N ) / ( M ) )+
     $                MIN( 1, MOD( ( N ), ( M ) ) ) )
*     ..
*     .. Executable Statements ..
      CONTXT = DESCA( CTXT_ )
      CALL BLACS_PINFO( MYID, NPROC )
      CALL BLACS_GRIDINFO( CONTXT, NPROW, NPCOL, MYPROW, MYPCOL )
      NROW = NUMROC( GNROWA, DESCA( MB_ ), MYPROW, DESCA( RSRC_ ),
     $       NPROW )
      LDA = MAX( 1, NROW )
      NCOL = INT( NFREE / LDA )
* compute the largest global number of columns
* that is acceptable
      IF( NROW.GE.1 ) THEN
* use a 'reverse' numroc
* such that ncol == numroc( gncolA...) and
* numroc( 1+gncolA ...) > ncol
         GNCOLA = NUMROCINV( NCOL, DESCA( NB_ ), MYPCOL, DESCA( CSRC_ ),
     $            NPCOL )
      ELSE
* this processor does not participate
* so should not be a limiting factor
* it is independent of how much lwork space it has
         GNCOLA = HUGE
      ENDIF
      SCOPE = 'A'
      TOP = ' '
      RA = -1
      CA = -1
      RCFLAG = -1
      RDEST = -1
      CDEST = -1
      CALL IGAMN2D( CONTXT, SCOPE, TOP, 1, 1, GNCOLA, 1, RA, CA, RCFLAG,
     $              RDEST, CDEST )
      CALL ASSERT( GNCOLA.NE.HUGE,
     $             '** calcolsize: internal error, myid ', MYID )
* double check
      NCOL = NUMROC( GNCOLA, DESCA( NB_ ), MYPCOL, DESCA( CSRC_ ),
     $       NPCOL )
      ISVALID = ( NCOL*NROW.LE.NFREE )
      IF( .NOT.ISVALID ) THEN
         WRITE( *, FMT = 9999 )NFREE, NROW, NCOL, GNCOLA, MYPROW,
     $      MYPCOL
 9999    FORMAT( '** calcolsize: nfree,nrow,ncol ', 3( 1X, I7 ),
     $         / ' gncolA, myprow,mypcol ', 3( 1X, I7 ) )
      ENDIF
      CALL ASSERT( ISVALID, '** calcolsize: internal error ', GNCOLA )
      RETURN
      END
