      SUBROUTINE PFMAXSIZE( ROWORCOL, LWORK, M, N, MB, NB, P0, Q0,
     $                      ICONTXT, INFO )
*
*
*  -- ScaLAPACK auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     Oct 10, 1996
*
*
* Purpose:
* =======
*
* PFMAXSIZE computes  largest scalapack matrix that fit within
* user supplied temporary work space.
*
* if roworcol .eq. 'C', then m is given, the largest n is computed
* if roworcol .eq. 'R', then n is given, the largest m is computed
*
* Most of the work is performed in NUMROCINV()
*
*
* Arguments
* =========
*
* roworcol      (global input) character
*               Specifies whether to compute m or n
*               = 'R': computes m, assume given n
*               = 'C': computes n, assume given m
*
*
* m             (global input/output) integer
*               The number of rows of desired matrix
*
* n             (global input/output) integer
*               The number of columns of desired matrix
*
*
* lwork         (local input) integer
*               Amount of local work space, may be different on
*               each processor.
*
*
* mb            (global input) integer
*               Row block size.
*
* nb            (global input) integer
*               Column block size.
*
*
*
*     .. Scalar Arguments ..
      CHARACTER          ROWORCOL
      INTEGER            ICONTXT, INFO, LWORK, M, MB, N, NB, P0, Q0
*     ..
*     .. Local Scalars ..
      LOGICAL            ISVALID
      INTEGER            CA, CDEST, CSRC, IVAL, MYPCOL, MYPROW, NCOL,
     $                   NPCOL, NPROW, NROW, RA, RCFLAG, RDEST, RSRC
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            NUMROC, NUMROCINV
      EXTERNAL           LSAME, NUMROC, NUMROCINV
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, IGAMN2D
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
      ISVALID = ( LSAME( ROWORCOL, 'R' ) .OR. LSAME( ROWORCOL, 'C' ) )
      IF( .NOT.ISVALID ) THEN
         INFO = -1
         RETURN
      ENDIF
      ISVALID = ( LWORK.GE.1 )
      IF( .NOT.ISVALID ) THEN
         INFO = -2
         RETURN
      ENDIF
      CALL BLACS_GRIDINFO( ICONTXT, NPROW, NPCOL, MYPROW, MYPCOL )
*
*    Overestimate storage requirements, just to be safe.
*
      RSRC = MYPROW
      CSRC = MYPCOL
      IF( LSAME( ROWORCOL, 'R' ) ) THEN
* n is given, generate m
         NCOL = NUMROC( N, NB, MYPCOL, CSRC, NPCOL )
         IF( NCOL.GE.1 ) THEN
            NROW = LWORK / NCOL
            IVAL = NUMROCINV( NROW, MB, MYPROW, RSRC, NPROW )
         ELSE
* this processor is not participating,
* use large value for global min operation
            IVAL = MAX( M, N )*MAX( M, N )
         ENDIF
      ELSE
         IF( LSAME( ROWORCOL, 'C' ) ) THEN
* m is given, generate n
            NROW = NUMROC( M, MB, MYPROW, RSRC, NPROW )
            IF( NROW.GE.1 ) THEN
               NCOL = LWORK / NROW
               IVAL = NUMROCINV( NCOL, NB, MYPCOL, CSRC, NPCOL )
            ELSE
* this processor is not participating,
* use large value for global min operation
               IVAL = MAX( M, N )*MAX( M, N )
            ENDIF
         ENDIF
      ENDIF
      RA = -1
      CA = -1
      RCFLAG = -1
      RDEST = -1
      CDEST = -1
      CALL IGAMN2D( ICONTXT, 'All', ' ', 1, 1, IVAL, 1, RA, CA, RCFLAG,
     $              RDEST, CDEST )
      IF( LSAME( ROWORCOL, 'R' ) ) THEN
         M = IVAL
      ELSE
         IF( LSAME( ROWORCOL, 'C' ) ) THEN
            N = IVAL
         ENDIF
      ENDIF
* cleanup and exit
      INFO = 0
      RETURN
      END
