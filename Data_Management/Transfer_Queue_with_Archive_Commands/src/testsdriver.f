      PROGRAM TESTSDRIVER
*
*
*  -- ScaLAPACK auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     Oct 10, 1996
*
*
*   Purpose:
*   ========
*   Main test program for out-of-core factorizations.
*
*   The program must be driven by a short data file.
*   An annotated example of a data file follows:
*
*
* ScaLAPACK out-of-core LU,QR,LL^t factorization input file
* testdriver.out
* 6                                     device out
* 3                                     number of factorizations
* QR
* LU
* LL
* 1                                     number of problem sizes
* 551                                   values of M
* 551                                   values of N
* 101                                   values of nrhs
* 28000                                         values of Asize
* 1                                     number of MB's and NB's
* 25                                    values of MB
* 25                                    values of NB
* 3                                     number of process grids
* 1 2 1                                 values of P
* 4 2 4                                 values of Q
*
*==========================
*
* The driver generates a (m by nrhs) random matrix  X and
* a (m by n) random matrix A
*
* The driver computes
*
*       (i) B = A * X
*       (ii) factorzation of A
*               QR              QR factorization
*               LU              LU factorization with pivoting
*               LL              LL' Cholesky factorization
*       (iii) Xbar = Solve(A, B)
*
*       For each column i=1:nrhs,
*
*       (iv) error(i)=norm(Xbar(:,i)-X(:,i),2)/
*                               (norm(X(:,i),2)*eps)
*       (v) residual(i)=norm(A*Xbar(:,i)-B(:,i),2)/
*                               (m*eps*norm(B(:,i),2))
*
*       maxerr = max( error(:) )
*       maxresidual = max( residual(:) )
*
*
*
*
* parameters
* ==========
*
*
*  TOTMEM   INTEGER, default = 2000000
*           TOTMEM is a machine-specific parameter indicating the
*           maximum amount of available memory in bytes.
*           The user should customize TOTMEM to his platform.  Remember
*           to leave room in memory for the operating system, the BLACS
*           buffer, etc.
*           Some experimenting with the maximum value of
*           TOTMEM may be required.
*
*  TOTIMEM  INTEGER, maximum amount of integer work space
*
*  MEM      dimension TOTMEM
*
*  IMEM     dimension TOTIMEM
*
*           All arrays used by SCALAPACK routines are allocated from
*           this array and referenced by pointers.  The integer IPA,
*           for example, is a pointer to the starting element of MEM for
*           the matrix A.
*
*
*  ASIZE    controls the amount of temporary storage allocated
*          for out-of-core operations.
*
*
*          Mininum storage is about space for two copies ScaLAPACK
*          m x nnb matrix panels for factorization.
*
*          Storage for one ScaLAPACK m x nnb matrix panel is required
*          for solvers such as PFxGEQRS.
*
*
*
*     .. Parameters ..
      INTEGER            DLEN_
      PARAMETER          ( DLEN_ = 9 )
      INTEGER            CTXT_, M_, N_, MB_, NB_
      PARAMETER          ( CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6 )
      INTEGER            RSRC_, CSRC_, LLD_
      PARAMETER          ( RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      INTEGER            IODEV_, SIZE_
      PARAMETER          ( IODEV_ = 10, SIZE_ = 11 )
      INTEGER            FLEN_
      PARAMETER          ( FLEN_ = 12 )
      INTEGER            IFACTOR, JFACTOR
      PARAMETER          ( IFACTOR = 3, JFACTOR = 1 )
      INTEGER            MAXNRHS
      PARAMETER          ( MAXNRHS = 1000 )
      LOGICAL            USE_MISALIGNED
      PARAMETER          ( USE_MISALIGNED = .false. )
      INTEGER            TOTMEM
      PARAMETER          ( TOTMEM = 6*1000*1000 )
      INTEGER            TOTIMEM
      PARAMETER          ( TOTIMEM = 600*1000 )
      INTEGER            NTESTS
      PARAMETER          ( NTESTS = 100 )
      INTEGER            LDFACT, LDMVAL, LDNVAL, LDMBVAL, LDNBVAL,
     $                   LDPVAL, LDQVAL, LDNRHSVAL, LDASIZEVAL
      PARAMETER          ( LDFACT = NTESTS, LDMVAL = NTESTS,
     $                   LDNVAL = NTESTS, LDMBVAL = NTESTS,
     $                   LDNBVAL = NTESTS, LDPVAL = NTESTS,
     $                   LDQVAL = NTESTS, LDNRHSVAL = NTESTS,
     $                   LDASIZEVAL = NTESTS )
*     ..
*     .. Local Scalars ..
      LOGICAL            ASIZEQUERY, ISVALID
      CHARACTER          AFORM, DIAG, FILETYPE, UPLO
      CHARACTER*2        FACT
      CHARACTER*80       FILENAME, OUTFILE, ROUT
      INTEGER            ANEED, ASIZE, CSRC, EL, I, IA, IAM, IB, IBEGIN,
     $                   ICNUM, ICOFF, ICTXT, IDUMMY, IFREE, INC, INFO,
     $                   INFO1, INFO2, IODEV, IODEVAOLD, IP_A, IP_AOLD,
     $                   IP_B, IP_IPIV, IP_R, IP_TAU, IP_WORK, IP_X, IR,
     $                   IRNUM, IROFF, ISEED, IX, J, JA, JB, JFREE, JR,
     $                   JX, K, KFAIL, KPASS, KSKIP, KTESTS, LDB, LDR,
     $                   LDX, LOCP, LOCQ, LWORK, M, MAXMN, MB,
     $                   MBNBCOUNT, MINMN, MMB, MYID, MYPCOL, MYPROW, N,
     $                   NB, NBRHS, NCOUNT, NFACT, NFREE, NFREE2,
     $                   NGRIDS, NMAT, NNB, NOUT, NPCOL, NPROCS, NPROW,
     $                   NRHS, P0, Q0, RSRC
      REAL               ALPHA, BETA, DZERO, MAXERR, MAXRESID, NORM2,
     $                   NORMB, NORMX, ONE, ZERO
      DOUBLE PRECISION   DEPS
*     ..
*     .. Local Arrays ..
      CHARACTER*2        FACTOR( LDFACT )
      INTEGER            ASIZEVAL( LDASIZEVAL ), DESCA( FLEN_ ),
     $                   DESCAOLD( FLEN_ ), DESCB( DLEN_ ),
     $                   DESCR( DLEN_ ), DESCX( DLEN_ ), IERR( 1 ),
     $                   IMEM( TOTIMEM ), MBVAL( LDMBVAL ),
     $                   MVAL( LDMVAL ), NBVAL( LDNBVAL ),
     $                   NRHSVAL( LDNRHSVAL ), NVAL( LDNVAL ),
     $                   PVAL( LDPVAL ), QVAL( LDQVAL )
      REAL               BNORM( MAXNRHS ), MEM( TOTMEM ),
     $                   XNORM( MAXNRHS )
      DOUBLE PRECISION   CTIME( 10 ), WTIME( 10 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      INTEGER            NUMROC
      REAL               PSLAMCH
      EXTERNAL           LSAMEN, NUMROC, PSLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_BARRIER, BLACS_EXIT, BLACS_GET,
     $                   BLACS_GRIDEXIT, BLACS_GRIDINFO, BLACS_GRIDINIT,
     $                   BLACS_PINFO, DESCINIT, DRIVERINFO, IGAMN2D,
     $                   IGSUM2D, LACLOSE, PFDESCINIT, PFSGEMM,
     $                   PFSGEQRF, PFSGEQRS, PFSGETRF, PFSGETRS,
     $                   PFSMATGEN, PFSPOTRF, PFSPOTRS, PSMATADD,
     $                   PSMATGEN, PSNRM2, SGAMX2D, SLBOOT, SLCOMBINE,
     $                   SLTIMER
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, CHAR, DBLE, INT, MAX, MIN, REAL
*     ..
*     .. Executable Statements ..
*
*   Get starting information
*
      ONE = REAL( 1 )
      ZERO = REAL( 0 )
      DZERO = DBLE( 0 )
      KTESTS = 0
      KPASS = 0
      KFAIL = 0
      KSKIP = 0
      CALL BLACS_PINFO( IAM, NPROCS )
      MYID = IAM
      CALL DRIVERINFO( OUTFILE, NOUT, NFACT, FACTOR, LDFACT, NMAT, MVAL,
     $                 LDMVAL, NVAL, LDNVAL, NRHSVAL, LDNRHSVAL,
     $                 ASIZEVAL, LDASIZEVAL, MBNBCOUNT, MBVAL, LDMBVAL,
     $                 NBVAL, LDNBVAL, NGRIDS, PVAL, LDPVAL, QVAL,
     $                 LDQVAL, IAM, NPROCS )
*
*   Loop over the different factorization types
*
      DO 160 I = 1, NFACT
         FACT = FACTOR( I )
*
*   Print headings.
*
         IF( IAM.EQ.0 ) THEN
            WRITE( NOUT, FMT = 9999 )
 9999       FORMAT( / / 'Note: ', /
     $            ' maximum error is norm(X - Xbar,2)/(eps*norm(X,2)) ',
     $            /
     $      ' maximum residual is norm(A*Xbar - B,2)/(eps*m*norm(B,2)) '
     $            , / ' where right hand side B = A*X ',
     $            / ' and solution Xbar = A \\ B ',
     $            / ' A is an m by m  random matrix ',
     $            / ' and X is an m by 1 random vector ', / / )
            WRITE( NOUT, FMT = * )
            IF( LSAMEN( 2, FACT, 'QR' ) ) THEN
               ROUT = 'PFxGEQRF'
               WRITE( NOUT, FMT = '(A)' )'QR factorization tests.'
            ELSE
               IF( LSAMEN( 2, FACT, 'LU' ) ) THEN
                  ROUT = 'PFxGETRF'
                  WRITE( NOUT, FMT = '(A)' )'LU factorzation tests.'
               ELSE
                  IF( LSAMEN( 2, FACT, 'LL' ) ) THEN
                     ROUT = 'PFxPOTRF'
                     WRITE( NOUT, FMT = '(A)' )
     $                  'LL^t factorization tests.'
                  ENDIF
               ENDIF
            ENDIF
            WRITE( NOUT, FMT = * )
            WRITE( NOUT, FMT = 9998 )
            WRITE( NOUT, FMT = 9997 )
            WRITE( NOUT, FMT = * )
         ENDIF
 9998    FORMAT(
     $      'TIME     M      N  MB  NB NRHS    P     Q Fact/Solve Time '
     $         , '     Error  Residual' )
 9997    FORMAT(
     $         '---- ------ ------ --- --- ----- ----- --------------- '
     $         , '----------- --------' )
*
*   Loop over different process grids
*
         DO 140 J = 1, NGRIDS
            NPROW = PVAL( J )
            NPCOL = QVAL( J )
*
*   Make sure grid information is correct.
*
            IERR( 1 ) = 0
            IF( NPROW.LT.1 ) THEN
               IERR( 1 ) = 1
               IF( IAM.EQ.0 ) THEN
                  WRITE( NOUT, FMT = 9994 )'grid', 'nprow', NPROW
               ENDIF
            ELSE
               IF( NPCOL.LT.1 ) THEN
                  IERR( 1 ) = 1
                  IF( IAM.EQ.0 ) THEN
                     WRITE( NOUT, FMT = 9994 )'grid', NPCOL, NPCOL
                  ENDIF
               ELSE
                  IF( NPROW*NPCOL.GT.NPROCS ) THEN
                     IERR( 1 ) = 1
                     IF( IAM.EQ.0 ) THEN
                        WRITE( NOUT, FMT = 9993 )NPROW*NPCOL, NPROCS
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
            IF( IERR( 1 ).NE.0 ) THEN
               IF( IAM.EQ.0 ) THEN
                  WRITE( NOUT, FMT = 9992 )'grid'
               ENDIF
               KSKIP = KSKIP + 1
               GOTO 140
            ENDIF
*
*   Define process grid.
*
            CALL BLACS_GET( -1, 0, ICTXT )
            CALL BLACS_GRIDINIT( ICTXT, 'Row-major', NPROW, NPCOL )
            CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYPROW, MYPCOL )
            DEPS = DBLE( PSLAMCH( ICTXT, 'eps' ) )
            IF( DEPS.EQ.REAL( 0 ) ) THEN
               WRITE( *, FMT = * )' error: PSLAMCH(eps) returns zero '
               STOP '** error ** '
            ENDIF
*
*   Go to bottom of loop if this case does not use my process.
*
            IF( ( MYPROW.GE.NPROW ) .OR. ( MYPCOL.GE.NPCOL ) ) THEN
               GOTO 140
            ENDIF
            DO 120 K = 1, NMAT
               M = MVAL( K )
               N = NVAL( K )
               ASIZE = ASIZEVAL( K )
               ASIZEQUERY = ( ASIZE.EQ.-1 )
               NRHS = NRHSVAL( K )
               NRHS = MAX( 1, MIN( MAXNRHS, NRHS ) )
*
*   Make sure matrix information is correct.
*
               IERR( 1 ) = 0
               IF( M.LT.1 ) THEN
                  IERR( 1 ) = 1
                  IF( IAM.EQ.0 ) THEN
                     WRITE( NOUT, FMT = 9994 )'MATRIX', 'M', M
                  ENDIF
               ELSE
                  IF( N.LT.1 ) THEN
                     IERR( 1 ) = 1
                     IF( IAM.EQ.0 ) THEN
                        WRITE( NOUT, FMT = 9994 )'MATRIX', 'N', N
                     ENDIF
                  ELSE
                     IF( M.NE.N ) THEN
                        IERR( 1 ) = 1
                        IF( IAM.EQ.0 ) THEN
                           WRITE( NOUT, FMT = * )'Require M .eq. N ',
     $                        ' M = ', M, ' N = ', N
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
*
*   Make sure no one had error.
*
               CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1, -1 )
               IF( IERR( 1 ).NE.0 ) THEN
                  IF( IAM.EQ.0 ) THEN
                     WRITE( NOUT, FMT = 9992 )'MATRIX'
                  ENDIF
                  KSKIP = KSKIP + 1
                  GOTO 120
               ENDIF
*
*   Loop over different blocking sizes.
*
               DO 100 EL = 1, MBNBCOUNT
                  MB = MBVAL( EL )
                  NB = NBVAL( EL )
*
*   Make sure mb is legal.
*
                  IERR( 1 ) = 0
                  IF( MB.LT.1 ) THEN
                     IERR( 1 ) = 1
                     IF( IAM.EQ.0 ) THEN
                        WRITE( NOUT, FMT = 9994 )'MB', 'MB', MB
                     ENDIF
                  ENDIF
                  CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1,
     $                          -1 )
                  IF( IERR( 1 ).NE.0 ) THEN
                     IF( IAM.EQ.0 ) THEN
                        WRITE( NOUT, FMT = 9992 )'MB'
                     ENDIF
                     KSKIP = KSKIP + 1
                     GOTO 100
                  ENDIF
*
*   Make sure nb is legal.
*
                  IERR( 1 ) = 0
                  IF( NB.LT.1 ) THEN
                     IERR( 1 ) = 1
                     IF( IAM.EQ.0 ) THEN
                        WRITE( NOUT, FMT = 9994 )'NB', 'NB', NB
                     ENDIF
                  ELSE
                     IF( NB.NE.MB ) THEN
                        IERR( 1 ) = 1
                        IF( IAM.EQ.0 ) THEN
                           WRITE( NOUT, FMT = * )'Require MB .eq. NB ',
     $                        ' MB = ', MB, ' NB = ', NB
                        ENDIF
                     ENDIF
                  ENDIF
                  CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1,
     $                          -1 )
                  IF( IERR( 1 ).NE.0 ) THEN
                     IF( IAM.EQ.0 ) THEN
                        WRITE( NOUT, FMT = 9992 )'NB'
                     ENDIF
                     KSKIP = KSKIP + 1
                     GOTO 100
                  ENDIF
*
*   Allocate storage.
*
                  CSRC = 0
                  RSRC = 0
                  IFREE = 1
                  JFREE = 1
                  IF( USE_MISALIGNED ) THEN
                     P0 = NPROW - 1
                     Q0 = NPCOL - 1
                  ELSE
                     P0 = 0
                     Q0 = 0
                  ENDIF
                  MAXMN = MAX( M, N )
                  LOCP = NUMROC( MAXMN, MB, MYPROW, P0, NPROW )
                  LOCQ = NUMROC( NRHS, NB, MYPCOL, Q0, NPCOL )
                  IP_B = IFREE
                  IFREE = IFREE + LOCP*LOCQ
                  IP_X = IFREE
                  IFREE = IFREE + LOCP*LOCQ
                  IP_R = IFREE
                  IFREE = IFREE + LOCP*LOCQ
                  IF( .NOT.ASIZEQUERY ) THEN
                     IP_A = IFREE
                     IFREE = IFREE + ASIZE
*
*                   Reuse buffer storage for Aold.
*
                     IP_AOLD = IP_A
                  ENDIF
*
*          Overestimate storage.
*
                  P0 = MYPROW
                  Q0 = MYPCOL
                  LOCQ = NUMROC( N, NB, MYPCOL, P0, NPCOL )
                  LOCP = NUMROC( M, MB, MYPROW, Q0, NPROW )
                  IF( LSAMEN( 2, FACT, 'QR' ) ) THEN
                     IP_TAU = IFREE
                     IFREE = IFREE + LOCQ
                  ELSE
                     IF( LSAMEN( 2, FACT, 'LU' ) ) THEN
                        IP_IPIV = JFREE
                        JFREE = JFREE + MAX( LOCP, LOCQ ) +
     $                          MAX( MB, NB )
                     ENDIF
                  ENDIF
                  NFREE = TOTMEM - IFREE + 1
                  IF( .NOT.ASIZEQUERY ) THEN
                     IP_WORK = IFREE
                     LWORK = NFREE
                  ENDIF
                  IERR( 1 ) = 0
                  IF( NFREE.LE.0 ) THEN
                     IERR( 1 ) = 1
                  ENDIF
                  CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1,
     $                          -1 )
                  IF( IERR( 1 ).NE.0 ) THEN
                     IF( IAM.EQ.0 ) THEN
                        WRITE( NOUT, FMT = 9992 )'totmem'
                     ENDIF
                     KSKIP = KSKIP + 1
                     GOTO 100
                  ENDIF
                  IERR( 1 ) = 0
                  NFREE2 = TOTIMEM - JFREE + 1
                  IF( NFREE2.LT.0 ) THEN
                     IERR( 1 ) = 1
                  ENDIF
                  CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1,
     $                          -1 )
                  IF( IERR( 1 ).NE.0 ) THEN
                     IF( IAM.EQ.0 ) THEN
                        WRITE( NOUT, FMT = 9992 )'totimem'
                     ENDIF
                     KSKIP = KSKIP + 1
                     GOTO 100
                  ENDIF
*
*   setup  descriptor.
*
                  INFO = 0
                  IODEV = 11
                  IODEVAOLD = 12
                  FILETYPE = 'D'
                  FILENAME = FACT( 1: 2 ) // '.dat' // CHAR( 0 )
                  MMB = IFACTOR*MB*NPROW
                  NNB = JFACTOR*NB*NPCOL
*
*                Note the use of p0,q0 to overestimate storage usage.
*
                  P0 = MYPROW
                  Q0 = MYPCOL
                  LOCP = NUMROC( MMB, MB, MYPROW, P0, NPROW )
                  LOCQ = NUMROC( NNB, NB, MYPCOL, Q0, NPCOL )
                  IF( .NOT.ASIZEQUERY ) THEN
                     IERR( 1 ) = 0
                     ANEED = 2*LOCP*LOCQ
                     ISVALID = ( ASIZE.GE.ANEED )
                     IF( .NOT.ISVALID ) THEN
                        IERR( 1 ) = 1
                     ENDIF
                     CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1,
     $                             -1 )
                     IF( IERR( 1 ).NE.0 ) THEN
                        IF( IAM.EQ.0 ) THEN
                           WRITE( NOUT, FMT = 9992 )'Asize'
                        ENDIF
                        KSKIP = KSKIP + 1
                        GOTO 100
                     ENDIF
                     IF( ( IAM.EQ.0 ) .AND. ( IERR( 1 ).NE.0 ) ) THEN
                        WRITE( NOUT, FMT = * )'Asize needs to be ',
     $                     ANEED
                     ENDIF
                  ENDIF
                  CALL PFDESCINIT( DESCA, M, N, MB, NB, RSRC, CSRC,
     $                             ICTXT, IODEV, FILETYPE, MMB, NNB,
     $                             ASIZE, FILENAME, INFO )
                  ISVALID = ( INFO.EQ.0 ) .OR.
     $                      ( ( ASIZEQUERY ) .AND. ( INFO.EQ.-13 ) )
                  IERR( 1 ) = 0
                  IF( .NOT.ISVALID ) THEN
                     IERR( 1 ) = 1
                  ENDIF
                  CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1,
     $                          -1 )
                  IF( IERR( 1 ).NE.0 ) THEN
                     IF( IAM.EQ.0 ) THEN
                        WRITE( NOUT, FMT = 9992 )'descriptor'
                     ENDIF
                     KSKIP = KSKIP + 1
                     GOTO 100
                  ENDIF
                  FILENAME = 'Aold.dat' // CHAR( 0 )
                  CALL PFDESCINIT( DESCAOLD, M, N, MB, NB, CSRC, RSRC,
     $                             ICTXT, IODEVAOLD, FILETYPE, MMB, NNB,
     $                             ASIZE, FILENAME, INFO )
                  ISVALID = ( INFO.EQ.0 ) .OR.
     $                      ( ( ASIZEQUERY ) .AND. ( INFO.EQ.-13 ) )
                  IERR( 1 ) = 0
                  IF( .NOT.ISVALID ) THEN
                     IERR( 1 ) = 1
                  ENDIF
                  CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1,
     $                          -1 )
                  IF( IERR( 1 ).NE.0 ) THEN
                     IF( IAM.EQ.0 ) THEN
                        WRITE( NOUT, FMT = 9992 )'descriptor'
                     ENDIF
                     KSKIP = KSKIP + 1
                     GOTO 100
                  ENDIF
                  IF( USE_MISALIGNED ) THEN
                     P0 = NPROW - 1
                     Q0 = NPCOL - 1
                  ELSE
                     P0 = 0
                     Q0 = 0
                  ENDIF
                  MAXMN = MAX( M, N )
                  MINMN = MIN( M, N )
                  LDX = NUMROC( MAXMN, MB, MYPROW, P0, NPROW )
                  LDB = LDX
                  LDR = LDX
                  CALL DESCINIT( DESCX, MAXMN, NRHS, MB, NB, P0, Q0,
     $                           ICTXT, LDX, INFO )
                  CALL DESCINIT( DESCB, MAXMN, NRHS, MB, NB, P0, Q0,
     $                           ICTXT, LDB, INFO )
                  CALL DESCINIT( DESCR, MAXMN, NRHS, MB, NB, P0, Q0,
     $                           ICTXT, LDR, INFO )
*
*   Generate out-of-core matrix A
*
                  AFORM = 'N'
                  DIAG = 'N'
                  IF( LSAMEN( 2, FACT, 'LL' ) ) THEN
                     DIAG = 'D'
                     AFORM = 'S'
                  ENDIF
                  ISEED = 13
                  INFO = 0
                  IF( ASIZEQUERY ) THEN
                     IP_A = IFREE
                     DESCA( SIZE_ ) = -1
                     MEM( IP_A ) = 0
                     CALL PFSMATGEN( DESCA( IODEV_ ), DESCA( CTXT_ ),
     $                               AFORM, DIAG, DESCA( M_ ),
     $                               DESCA( N_ ), DESCA( MB_ ),
     $                               DESCA( NB_ ), MEM( IP_A ),
     $                               DESCA( SIZE_ ), DESCA( RSRC_ ),
     $                               DESCA( CSRC_ ), ISEED, MYPROW,
     $                               MYPCOL, NPROW, NPCOL, INFO )
                     ANEED = INT( DBLE( 100 )*DEPS+ABS( MEM( IP_A ) ) )
                     DESCA( SIZE_ ) = ANEED
                     IFREE = IFREE + ANEED
                     ISVALID = ( IFREE.LE.TOTMEM )
                     IERR( 1 ) = 0
                     IF( .NOT.ISVALID ) THEN
                        IERR( 1 ) = 1
                     ENDIF
                     CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1,
     $                             -1 )
                     IF( IERR( 1 ).NE.0 ) THEN
                        IF( IAM.EQ.0 ) THEN
                           WRITE( NOUT, FMT = 9992 )'Asize'
                        ENDIF
                        KSKIP = KSKIP + 1
                        GOTO 100
                     ENDIF
                     CALL PFSMATGEN( DESCA( IODEV_ ), DESCA( CTXT_ ),
     $                               AFORM, DIAG, DESCA( M_ ),
     $                               DESCA( N_ ), DESCA( MB_ ),
     $                               DESCA( NB_ ), MEM( IP_A ),
     $                               DESCA( SIZE_ ), DESCA( RSRC_ ),
     $                               DESCA( CSRC_ ), ISEED, MYPROW,
     $                               MYPCOL, NPROW, NPCOL, INFO )
                     IFREE = IFREE - ANEED
                  ELSE
                     CALL PFSMATGEN( DESCA( IODEV_ ), DESCA( CTXT_ ),
     $                               AFORM, DIAG, DESCA( M_ ),
     $                               DESCA( N_ ), DESCA( MB_ ),
     $                               DESCA( NB_ ), MEM( IP_A ),
     $                               DESCA( SIZE_ ), DESCA( RSRC_ ),
     $                               DESCA( CSRC_ ), ISEED, MYPROW,
     $                               MYPCOL, NPROW, NPCOL, INFO )
                  ENDIF
                  IERR( 1 ) = INFO
                  CALL IGAMN2D( ICTXT, 'All', ' ', 1, 1, IERR, 1,
     $                          IDUMMY, IDUMMY, -1, IDUMMY, IDUMMY )
                  IF( IERR( 1 ).NE.0 ) THEN
                     IF( IAM.EQ.0 ) THEN
                        WRITE( NOUT, FMT = * )'PFxMATGEN returns info ',
     $                     IERR( 1 )
                     ENDIF
                     KSKIP = KSKIP + 1
                     GOTO 90
                  ENDIF
*
*       Generate another copy for Aold.
*
                  ISEED = 13
                  INFO = 0
                  IF( ASIZEQUERY ) THEN
                     IP_AOLD = IFREE
                     DESCAOLD( SIZE_ ) = -1
                     MEM( IP_AOLD ) = 0
                     CALL PFSMATGEN( DESCAOLD( IODEV_ ),
     $                               DESCAOLD( CTXT_ ), AFORM, DIAG,
     $                               DESCAOLD( M_ ), DESCAOLD( N_ ),
     $                               DESCAOLD( MB_ ), DESCAOLD( NB_ ),
     $                               MEM( IP_AOLD ), DESCAOLD( SIZE_ ),
     $                               DESCAOLD( RSRC_ ),
     $                               DESCAOLD( CSRC_ ), ISEED, MYPROW,
     $                               MYPCOL, NPROW, NPCOL, INFO )
                     ANEED = INT( DBLE( 100 )*DEPS+
     $                       ABS( MEM( IP_AOLD ) ) )
                     DESCAOLD( SIZE_ ) = ANEED
                     IFREE = IFREE + ANEED
                     ISVALID = ( IFREE.LE.TOTMEM )
                     IERR( 1 ) = 0
                     IF( .NOT.ISVALID ) THEN
                        IERR( 1 ) = 1
                     ENDIF
                     CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1,
     $                             -1 )
                     IF( IERR( 1 ).NE.0 ) THEN
                        IF( IAM.EQ.0 ) THEN
                           WRITE( NOUT, FMT = 9992 )'Asize'
                        ENDIF
                        KSKIP = KSKIP + 1
                        GOTO 100
                     ENDIF
                     CALL PFSMATGEN( DESCAOLD( IODEV_ ),
     $                               DESCAOLD( CTXT_ ), AFORM, DIAG,
     $                               DESCAOLD( M_ ), DESCAOLD( N_ ),
     $                               DESCAOLD( MB_ ), DESCAOLD( NB_ ),
     $                               MEM( IP_A ), DESCAOLD( SIZE_ ),
     $                               DESCAOLD( RSRC_ ),
     $                               DESCAOLD( CSRC_ ), ISEED, MYPROW,
     $                               MYPCOL, NPROW, NPCOL, INFO )
                     IFREE = IFREE - ANEED
                  ELSE
                     CALL PFSMATGEN( DESCAOLD( IODEV_ ),
     $                               DESCAOLD( CTXT_ ), AFORM, DIAG,
     $                               DESCAOLD( M_ ), DESCAOLD( N_ ),
     $                               DESCAOLD( MB_ ), DESCAOLD( NB_ ),
     $                               MEM( IP_AOLD ), DESCAOLD( SIZE_ ),
     $                               DESCAOLD( RSRC_ ),
     $                               DESCAOLD( CSRC_ ), ISEED, MYPROW,
     $                               MYPCOL, NPROW, NPCOL, INFO )
                  ENDIF
                  IERR( 1 ) = INFO
                  CALL IGAMN2D( ICTXT, 'All', ' ', 1, 1, IERR, 1,
     $                          IDUMMY, IDUMMY, -1, IDUMMY, IDUMMY )
                  IF( IERR( 1 ).NE.0 ) THEN
                     IF( IAM.EQ.0 ) THEN
                        WRITE( NOUT, FMT = * )'PFxMATGEN returns info ',
     $                     IERR( 1 )
                     ENDIF
                     KSKIP = KSKIP + 1
                     GOTO 90
                  ENDIF
*
*   Generate solution X.
*
                  AFORM = 'N'
                  DIAG = 'N'
                  IROFF = 0
                  ICOFF = 0
                  IRNUM = NUMROC( DESCX( M_ ), DESCX( MB_ ), MYPROW,
     $                    DESCX( RSRC_ ), NPROW )
                  ICNUM = NUMROC( DESCX( N_ ), DESCX( NB_ ), MYPCOL,
     $                    DESCX( CSRC_ ), NPCOL )
                  ISEED = 999999
                  CALL PSMATGEN( DESCX( CTXT_ ), AFORM, DIAG,
     $                           DESCX( M_ ), DESCX( N_ ), DESCX( MB_ ),
     $                           DESCX( NB_ ), MEM( IP_X ),
     $                           DESCX( LLD_ ), DESCX( RSRC_ ),
     $                           DESCX( CSRC_ ), ISEED, IROFF, IRNUM,
     $                           ICOFF, ICNUM, MYPROW, MYPCOL, NPROW,
     $                           NPCOL )
*
*  Compute norm of X.
*
                  IX = 1
                  INC = 1
                  DO 10 JX = 1, NRHS
                     NORM2 = -9999
                     CALL PSNRM2( MIN( M, N ), NORM2, MEM( IP_X ), IX,
     $                            JX, DESCX, INC )
                     CALL SGAMX2D( ICTXT, 'A', ' ', 1, 1, NORM2, 1, 0,
     $                             0, -1, -1, -1 )
                     XNORM( JX ) = NORM2
   10             CONTINUE
   20             CONTINUE
                  IX = 1
                  JX = 1
*
*   Generate rhs B <- A*X
*
*   A is m x n, X is n x nrhs, B is m x nrhs
*
                  ALPHA = ONE
                  BETA = ZERO
                  IA = 1
                  JA = 1
                  IX = 1
                  JX = 1
                  IB = 1
                  JB = 1
                  IF( ASIZEQUERY ) THEN
                     IP_A = IFREE
                     DESCAOLD( SIZE_ ) = -1
                     MEM( IP_A ) = 0
                     CALL PFSGEMM( 'NotransA', 'NotransX', M, NRHS, N,
     $                             ALPHA, MEM( IP_A ), IA, JA, DESCAOLD,
     $                             MEM( IP_X ), IX, JX, DESCX, BETA,
     $                             MEM( IP_B ), IB, JB, DESCB, INFO )
                     ANEED = INT( DBLE( 100 )*DEPS+ABS( MEM( IP_A ) ) )
                     DESCAOLD( SIZE_ ) = ANEED
                     IFREE = IFREE + ANEED
                     ISVALID = ( IFREE.LE.TOTMEM )
                     IERR( 1 ) = 0
                     IF( .NOT.ISVALID ) THEN
                        IERR( 1 ) = 1
                     ENDIF
                     CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1,
     $                             -1 )
                     IF( IERR( 1 ).NE.0 ) THEN
                        IF( IAM.EQ.0 ) THEN
                           WRITE( NOUT, FMT = 9992 )'Asize'
                        ENDIF
                        KSKIP = KSKIP + 1
                        GOTO 100
                     ENDIF
                     CALL PFSGEMM( 'NotransA', 'NotransX', M, NRHS, N,
     $                             ALPHA, MEM( IP_A ), IA, JA, DESCAOLD,
     $                             MEM( IP_X ), IX, JX, DESCX, BETA,
     $                             MEM( IP_B ), IB, JB, DESCB, INFO )
                     IFREE = IFREE - ANEED
                  ELSE
                     CALL PFSGEMM( 'NotransA', 'NotransX', M, NRHS, N,
     $                             ALPHA, MEM( IP_A ), IA, JA, DESCAOLD,
     $                             MEM( IP_X ), IX, JX, DESCX, BETA,
     $                             MEM( IP_B ), IB, JB, DESCB, INFO )
                  ENDIF
                  IERR( 1 ) = INFO
                  CALL IGAMN2D( ICTXT, 'All', ' ', 1, 1, IERR, 1,
     $                          IDUMMY, IDUMMY, -1, IDUMMY, IDUMMY )
                  IF( IERR( 1 ).NE.0 ) THEN
                     IF( IAM.EQ.0 ) THEN
                        WRITE( NOUT, FMT = * )'PFxGEMM returns info ',
     $                     IERR( 1 )
                     ENDIF
                     KSKIP = KSKIP + 1
                     GOTO 90
                  ENDIF
*
*  Compute norm of B.
*
                  DO 30 JB = 1, NRHS
                     NORM2 = -9999
                     CALL PSNRM2( MIN( M, N ), NORM2, MEM( IP_B ), IB,
     $                            JB, DESCB, INC )
                     CALL SGAMX2D( ICTXT, 'A', ' ', 1, 1, NORM2, 1, 0,
     $                             0, -1, -1, -1 )
                     BNORM( JB ) = NORM2
   30             CONTINUE
   40             CONTINUE
                  IB = 1
                  JB = 1
*
*   Overwrite R with B.
*
                  ALPHA = ONE
                  BETA = ZERO
                  IR = 1
                  JR = 1
                  CALL PSMATADD( M, NRHS, ALPHA, MEM( IP_B ), IB, JB,
     $                           DESCB, BETA, MEM( IP_R ), IR, JR,
     $                           DESCR )
*
*   Generate factorization and solve
*
                  CALL SLBOOT
                  CALL BLACS_BARRIER( ICTXT, 'All' )
                  IF( LSAMEN( 2, FACT, 'QR' ) ) THEN
*
*                    QR factorization.
*
                     INFO = 0
                     CALL SLTIMER( 1 )
                     IF( ASIZEQUERY ) THEN
                        IP_A = IFREE
                        DESCA( SIZE_ ) = -1
                        IP_WORK = IFREE + 1
                        LWORK = -1
*
*                        Query storage requirement.
*
                        CALL PFSGEQRF( M, N, MEM( IP_A ), IA, JA, DESCA,
     $                                 MEM( IP_TAU ), MEM( IP_WORK ),
     $                                 LWORK, INFO )
                        LWORK = INT( DBLE( 100 )*DEPS+
     $                          ABS( MEM( IP_WORK ) ) )
                        ANEED = INT( DBLE( 100 )*DEPS+
     $                          ABS( MEM( IP_A ) ) )
                        DESCA( SIZE_ ) = ANEED
                        IP_A = IFREE
                        IFREE = IFREE + ANEED
                        IP_WORK = IFREE
                        IFREE = IFREE + LWORK
                        ISVALID = ( IFREE.LE.TOTMEM )
                        IERR( 1 ) = 0
                        IF( .NOT.ISVALID ) THEN
                           IERR( 1 ) = 1
                        ENDIF
                        CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1,
     $                                -1, -1 )
                        IF( IERR( 1 ).NE.0 ) THEN
                           IF( IAM.EQ.0 ) THEN
                              WRITE( NOUT, FMT = 9992 )'Asize'
                           ENDIF
                           KSKIP = KSKIP + 1
                           GOTO 100
                        ENDIF
                        CALL PFSGEQRF( M, N, MEM( IP_A ), IA, JA, DESCA,
     $                                 MEM( IP_TAU ), MEM( IP_WORK ),
     $                                 LWORK, INFO )
*
*                        Deallocate storage.
*
                        IFREE = IFREE - LWORK
                        IFREE = IFREE - ANEED
                     ELSE
                        CALL PFSGEQRF( M, N, MEM( IP_A ), IA, JA, DESCA,
     $                                 MEM( IP_TAU ), MEM( IP_WORK ),
     $                                 LWORK, INFO )
                     ENDIF
                     CALL SLTIMER( 1 )
                     IERR( 1 ) = INFO
                     CALL IGAMN2D( ICTXT, 'All', ' ', 1, 1, IERR, 1,
     $                             IDUMMY, IDUMMY, -1, IDUMMY, IDUMMY )
                     IF( IERR( 1 ).NE.0 ) THEN
                        IF( IAM.EQ.0 ) THEN
                           WRITE( NOUT, FMT = * )
     $                        'PFxGEQRF returns info ', IERR( 1 )
                        ENDIF
                        KSKIP = KSKIP + 1
                        GOTO 90
                     ENDIF
                     CALL BLACS_BARRIER( ICTXT, 'All' )
                     CALL SLTIMER( 2 )
                     IF( ASIZEQUERY ) THEN
*
*                        Query storage requirements.
*
                        IP_A = IFREE
                        DESCA( SIZE_ ) = -1
                        IP_WORK = IFREE + 1
                        LWORK = -1
                        MEM( IP_A ) = 0
                        MEM( IP_WORK ) = 0
                        CALL PFSGEQRS( M, N, NRHS, MEM( IP_A ), IA, JA,
     $                                 DESCA, MEM( IP_TAU ),
     $                                 MEM( IP_R ), IR, JR, DESCR,
     $                                 MEM( IP_WORK ), LWORK, INFO )
                        ANEED = INT( DBLE( 100 )*DEPS+
     $                          ABS( MEM( IP_A ) ) )
                        DESCA( SIZE_ ) = ANEED
                        LWORK = INT( DBLE( 100 )*DEPS+
     $                          ABS( MEM( IP_WORK ) ) )
                        IP_A = IFREE
                        IFREE = IFREE + ANEED
                        IP_WORK = IFREE
                        IFREE = IFREE + LWORK
                        ISVALID = ( IFREE.LE.TOTMEM )
                        IERR( 1 ) = 0
                        IF( .NOT.ISVALID ) THEN
                           IERR( 1 ) = 1
                        ENDIF
                        CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1,
     $                                -1, -1 )
                        IF( IERR( 1 ).NE.0 ) THEN
                           IF( IAM.EQ.0 ) THEN
                              WRITE( NOUT, FMT = 9992 )'Asize'
                           ENDIF
                           KSKIP = KSKIP + 1
                           GOTO 100
                        ENDIF
                        CALL PFSGEQRS( M, N, NRHS, MEM( IP_A ), IA, JA,
     $                                 DESCA, MEM( IP_TAU ),
     $                                 MEM( IP_R ), IR, JR, DESCR,
     $                                 MEM( IP_WORK ), LWORK, INFO )
*
*                        Deallocate storage.
*
                        IFREE = IFREE - LWORK
                        IFREE = IFREE - ANEED
                     ELSE
                        CALL PFSGEQRS( M, N, NRHS, MEM( IP_A ), IA, JA,
     $                                 DESCA, MEM( IP_TAU ),
     $                                 MEM( IP_R ), IR, JR, DESCR,
     $                                 MEM( IP_WORK ), LWORK, INFO )
                     ENDIF
                     CALL SLTIMER( 2 )
                     IERR( 1 ) = INFO
                     CALL IGAMN2D( ICTXT, 'All', ' ', 1, 1, IERR, 1,
     $                             IDUMMY, IDUMMY, -1, IDUMMY, IDUMMY )
                     IF( IERR( 1 ).NE.0 ) THEN
                        IF( IAM.EQ.0 ) THEN
                           WRITE( NOUT, FMT = * )
     $                        'PFxGEQRS returns info ', IERR( 1 )
                        ENDIF
                        KSKIP = KSKIP + 1
                        GOTO 90
                     ENDIF
                  ELSE
                     IF( LSAMEN( 2, FACT, 'LL' ) ) THEN
*
*                    Cholesky factorization.
*
                        INFO = 0
                        UPLO = 'Lower'
                        CALL SLTIMER( 1 )
                        IF( ASIZEQUERY ) THEN
*
*                         Query storage requirements.
*
                           IP_A = IFREE
                           DESCA( SIZE_ ) = -1
                           MEM( IP_A ) = 0
                           CALL PFSPOTRF( UPLO, MIN( M, N ),
     $                                    MEM( IP_A ), IA, JA, DESCA,
     $                                    INFO )
                           ANEED = INT( DBLE( 100 )*DEPS+
     $                             ABS( MEM( IP_A ) ) )
                           DESCA( SIZE_ ) = ANEED
                           IFREE = IFREE + ANEED
                           ISVALID = ( IFREE.LE.TOTMEM )
                           IERR( 1 ) = 0
                           IF( .NOT.ISVALID ) THEN
                              IERR( 1 ) = 1
                           ENDIF
                           CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR,
     $                                   1, -1, -1 )
                           IF( IERR( 1 ).NE.0 ) THEN
                              IF( IAM.EQ.0 ) THEN
                                 WRITE( NOUT, FMT = 9992 )'Asize'
                              ENDIF
                              KSKIP = KSKIP + 1
                              GOTO 100
                           ENDIF
                           CALL PFSPOTRF( UPLO, MIN( M, N ),
     $                                    MEM( IP_A ), IA, JA, DESCA,
     $                                    INFO )
                           IFREE = IFREE - ANEED
                        ELSE
                           CALL PFSPOTRF( UPLO, MIN( M, N ),
     $                                    MEM( IP_A ), IA, JA, DESCA,
     $                                    INFO )
                        ENDIF
                        CALL SLTIMER( 1 )
                        IERR( 1 ) = INFO
                        CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1,
     $                                -1, 0 )
                        IF( IERR( 1 ).NE.0 ) THEN
                           IF( IAM.EQ.0 ) THEN
                              WRITE( NOUT, FMT = * )
     $                           'PFxPOTRF returns info ', IERR( 1 )
                           ENDIF
                           KSKIP = KSKIP + 1
                           GOTO 90
                        ENDIF
                        CALL BLACS_BARRIER( ICTXT, 'All' )
                        CALL SLTIMER( 2 )
                        IF( ASIZEQUERY ) THEN
                           IP_A = IFREE
                           DESCA( SIZE_ ) = -1
                           MEM( IP_A ) = 0
                           CALL PFSPOTRS( UPLO, MIN( M, N ), NRHS,
     $                                    MEM( IP_A ), IA, JA, DESCA,
     $                                    MEM( IP_R ), IR, JR, DESCR,
     $                                    INFO )
                           ANEED = INT( DBLE( 100 )*DEPS+
     $                             ABS( MEM( IP_A ) ) )
                           DESCA( SIZE_ ) = ANEED
                           IFREE = IFREE + ANEED
                           ISVALID = ( IFREE.LE.TOTMEM )
                           IERR( 1 ) = 0
                           IF( .NOT.ISVALID ) THEN
                              IERR( 1 ) = 1
                           ENDIF
                           CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR,
     $                                   1, -1, -1 )
                           IF( IERR( 1 ).NE.0 ) THEN
                              IF( IAM.EQ.0 ) THEN
                                 WRITE( NOUT, FMT = 9992 )'Asize'
                              ENDIF
                              KSKIP = KSKIP + 1
                              GOTO 100
                           ENDIF
                           CALL PFSPOTRS( UPLO, MIN( M, N ), NRHS,
     $                                    MEM( IP_A ), IA, JA, DESCA,
     $                                    MEM( IP_R ), IR, JR, DESCR,
     $                                    INFO )
                           IFREE = IFREE - ANEED
                        ELSE
                           CALL PFSPOTRS( UPLO, MIN( M, N ), NRHS,
     $                                    MEM( IP_A ), IA, JA, DESCA,
     $                                    MEM( IP_R ), IR, JR, DESCR,
     $                                    INFO )
                        ENDIF
                        CALL SLTIMER( 2 )
                        IERR( 1 ) = INFO
                        CALL IGAMN2D( ICTXT, 'All', ' ', 1, 1, IERR, 1,
     $                                IDUMMY, IDUMMY, -1, IDUMMY,
     $                                IDUMMY )
                        IF( IERR( 1 ).NE.0 ) THEN
                           IF( IAM.EQ.0 ) THEN
                              WRITE( NOUT, FMT = * )
     $                           'PFxPOTRS returns info ', IERR( 1 )
                           ENDIF
                           KSKIP = KSKIP + 1
                           GOTO 90
                        ENDIF
                     ELSE
                        IF( LSAMEN( 2, FACT, 'LU' ) ) THEN
*
*                    LU factorization.
*
                           INFO = 0
                           CALL SLTIMER( 1 )
                           IF( ASIZEQUERY ) THEN
                              IP_A = IFREE
                              DESCA( SIZE_ ) = -1
                              MEM( IP_A ) = 0
                              CALL PFSGETRF( M, N, MEM( IP_A ), IA, JA,
     $                                       DESCA, IMEM( IP_IPIV ),
     $                                       INFO )
                              ANEED = INT( DBLE( 100 )*DEPS+
     $                                ABS( MEM( IP_A ) ) )
                              DESCA( SIZE_ ) = ANEED
                              IP_A = IFREE
                              IFREE = IFREE + ANEED
                              ISVALID = ( IFREE.LE.TOTMEM )
                              IERR( 1 ) = 0
                              IF( .NOT.ISVALID ) THEN
                                 IERR( 1 ) = 1
                              ENDIF
                              CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1,
     $                                      IERR, 1, -1, -1 )
                              IF( IERR( 1 ).NE.0 ) THEN
                                 IF( IAM.EQ.0 ) THEN
                                    WRITE( NOUT, FMT = 9992 )'Asize'
                                 ENDIF
                                 KSKIP = KSKIP + 1
                                 GOTO 100
                              ENDIF
                              CALL PFSGETRF( M, N, MEM( IP_A ), IA, JA,
     $                                       DESCA, IMEM( IP_IPIV ),
     $                                       INFO )
                              IFREE = IFREE - ANEED
                           ELSE
                              CALL PFSGETRF( M, N, MEM( IP_A ), IA, JA,
     $                                       DESCA, IMEM( IP_IPIV ),
     $                                       INFO )
                           ENDIF
                           CALL SLTIMER( 1 )
                           IERR( 1 ) = INFO
                           CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR,
     $                                   1, -1, 0 )
                           IF( IERR( 1 ).NE.0 ) THEN
                              IF( IAM.EQ.0 ) THEN
                                 WRITE( NOUT, FMT = * )
     $                              'PFxGETRF returns info ', IERR( 1 )
                              ENDIF
                              KSKIP = KSKIP + 1
                              GOTO 90
                           ENDIF
                           CALL BLACS_BARRIER( ICTXT, 'All' )
                           CALL SLTIMER( 2 )
                           IF( ASIZEQUERY ) THEN
                              IP_A = IFREE
                              DESCA( SIZE_ ) = -1
                              MEM( IP_A ) = 0
                              CALL PFSGETRS( 'NotransA', M, NRHS,
     $                                       MEM( IP_A ), IA, JA, DESCA,
     $                                       IMEM( IP_IPIV ),
     $                                       MEM( IP_R ), IR, JR, DESCR,
     $                                       INFO )
                              ANEED = INT( DBLE( 100 )*DEPS+
     $                                ABS( MEM( IP_A ) ) )
                              DESCA( SIZE_ ) = ANEED
                              IP_A = IFREE
                              IFREE = IFREE + ANEED
                              ISVALID = ( IFREE.LE.TOTMEM )
                              IERR( 1 ) = 0
                              IF( .NOT.ISVALID ) THEN
                                 IERR( 1 ) = 1
                              ENDIF
                              CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1,
     $                                      IERR, 1, -1, -1 )
                              IF( IERR( 1 ).NE.0 ) THEN
                                 IF( IAM.EQ.0 ) THEN
                                    WRITE( NOUT, FMT = 9992 )'Asize'
                                 ENDIF
                                 KSKIP = KSKIP + 1
                                 GOTO 100
                              ENDIF
                              CALL PFSGETRS( 'NotransA', M, NRHS,
     $                                       MEM( IP_A ), IA, JA, DESCA,
     $                                       IMEM( IP_IPIV ),
     $                                       MEM( IP_R ), IR, JR, DESCR,
     $                                       INFO )
                              IFREE = IFREE - ANEED
                           ELSE
                              CALL PFSGETRS( 'NotransA', M, NRHS,
     $                                       MEM( IP_A ), IA, JA, DESCA,
     $                                       IMEM( IP_IPIV ),
     $                                       MEM( IP_R ), IR, JR, DESCR,
     $                                       INFO )
                           ENDIF
                           CALL SLTIMER( 2 )
                           IERR( 1 ) = INFO
                           CALL IGAMN2D( ICTXT, 'All', ' ', 1, 1, IERR,
     $                                   1, IDUMMY, IDUMMY, -1, IDUMMY,
     $                                   IDUMMY )
                           IF( IERR( 1 ).NE.0 ) THEN
                              IF( IAM.EQ.0 ) THEN
                                 WRITE( NOUT, FMT = * )
     $                              'PFxGETRS returns info ', IERR( 1 )
                              ENDIF
                              KSKIP = KSKIP + 1
                              GOTO 90
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
*
*   Compute error vector.
*   X <- X - R
*
                  ALPHA = -ONE
                  BETA = ONE
                  CALL PSMATADD( MIN( M, N ), NRHS, ALPHA, MEM( IP_R ),
     $                           IR, JR, DESCR, BETA, MEM( IP_X ), IX,
     $                           JX, DESCX )
                  MAXERR = 0
                  IX = 1
                  INC = 1
                  DO 50 JX = 1, NRHS
                     NORM2 = -9999
                     CALL PSNRM2( MIN( M, N ), NORM2, MEM( IP_X ), IX,
     $                            JX, DESCX, INC )
                     CALL SGAMX2D( ICTXT, 'A', ' ', 1, 1, NORM2, 1, 0,
     $                             0, -1, -1, -1 )
                     NORMX = XNORM( JX )
                     IF( NORMX.GT.DEPS ) THEN
                        MAXERR = MAX( MAXERR, NORM2 / NORMX )
                     ENDIF
                     IF( IAM.EQ.0 ) THEN
                        WRITE( NOUT, FMT = * )'error for ', JX,
     $                     '-th col is ', REAL( NORM2 )
                     ENDIF
   50             CONTINUE
   60             CONTINUE
*
*       Normalize error.
*
                  MAXERR = DBLE( MAXERR ) / ( DEPS )
*
*   Compute residual error
*   B <- B - A*R
*
                  IR = 1
                  JR = 1
                  IB = 1
                  JB = 1
                  IA = 1
                  JA = 1
                  ALPHA = -ONE
                  BETA = ONE
                  IF( ASIZEQUERY ) THEN
                     IP_A = IFREE
                     DESCAOLD( SIZE_ ) = -1
                     MEM( IP_A ) = 0
                     CALL PFSGEMM( 'NotransA', 'NotransR', M, NRHS, N,
     $                             ALPHA, MEM( IP_A ), IA, JA, DESCAOLD,
     $                             MEM( IP_R ), IR, JR, DESCR, BETA,
     $                             MEM( IP_B ), IB, JB, DESCB, INFO )
                     ANEED = INT( DBLE( 100 )*DEPS+ABS( MEM( IP_A ) ) )
                     DESCAOLD( SIZE_ ) = ANEED
                     IFREE = IFREE + ANEED
                     ISVALID = ( IFREE.LE.TOTMEM )
                     IERR( 1 ) = 0
                     IF( .NOT.ISVALID ) THEN
                        IERR( 1 ) = 1
                     ENDIF
                     CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1,
     $                             -1 )
                     IF( IERR( 1 ).NE.0 ) THEN
                        IF( IAM.EQ.0 ) THEN
                           WRITE( NOUT, FMT = 9992 )'Asize'
                        ENDIF
                        KSKIP = KSKIP + 1
                        GOTO 100
                     ENDIF
                     CALL PFSGEMM( 'NotransA', 'NotransR', M, NRHS, N,
     $                             ALPHA, MEM( IP_A ), IA, JA, DESCAOLD,
     $                             MEM( IP_R ), IR, JR, DESCR, BETA,
     $                             MEM( IP_B ), IB, JB, DESCB, INFO )
                     IFREE = IFREE - ANEED
                  ELSE
                     CALL PFSGEMM( 'NotransA', 'NotransR', M, NRHS, N,
     $                             ALPHA, MEM( IP_A ), IA, JA, DESCAOLD,
     $                             MEM( IP_R ), IR, JR, DESCR, BETA,
     $                             MEM( IP_B ), IB, JB, DESCB, INFO )
                  ENDIF
                  IERR( 1 ) = INFO
                  CALL IGAMN2D( ICTXT, 'All', ' ', 1, 1, IERR, 1,
     $                          IDUMMY, IDUMMY, -1, IDUMMY, IDUMMY )
                  IF( IERR( 1 ).NE.0 ) THEN
                     IF( IAM.EQ.0 ) THEN
                        WRITE( NOUT, FMT = * )'PFxGEMM returns info ',
     $                     IERR( 1 )
                     ENDIF
                     KSKIP = KSKIP + 1
                     GOTO 90
                  ENDIF
                  MAXRESID = 0
                  IB = 1
                  INC = 1
                  DO 70 JB = 1, NRHS
                     NORM2 = -9999
                     CALL PSNRM2( MIN( M, N ), NORM2, MEM( IP_B ), IB,
     $                            JB, DESCB, INC )
                     CALL SGAMX2D( ICTXT, 'A', ' ', 1, 1, NORM2, 1, 0,
     $                             0, -1, -1, -1 )
                     NORMB = BNORM( JB )
                     IF( NORMB.GT.DEPS ) THEN
                        MAXRESID = MAX( MAXRESID, NORM2 / NORMB )
                     ENDIF
   70             CONTINUE
   80             CONTINUE
*
*                Normalize by size of matrix.
*
                  MAXRESID = DBLE( MAXRESID ) / ( DEPS*DBLE( M ) )
*
*   Gather maximum of all cpu and wall clock timings
*
                  NCOUNT = 2
                  IBEGIN = 1
                  CALL SLCOMBINE( ICTXT, 'All', '>', 'W', NCOUNT,
     $                            IBEGIN, WTIME )
                  CALL SLCOMBINE( ICTXT, 'All', '>', 'C', NCOUNT,
     $                            IBEGIN, CTIME )
                  IF( WTIME( 1 ).GE.DBLE( 0 ) ) THEN
                     NBRHS = NB
                     IF( IAM.EQ.0 ) THEN
                        WRITE( NOUT, FMT = 9995 )'WALL', M, N, MB, NB,
     $                     NRHS, NPROW, NPCOL, WTIME( 1 ), WTIME( 2 ),
     $                     REAL( MAXERR ), REAL( MAXRESID )
                     ENDIF
                  ENDIF
                  IF( CTIME( 1 ).GE.DBLE( 0 ) ) THEN
                     IF( IAM.EQ.0 ) THEN
                        WRITE( NOUT, FMT = 9995 )'CPU ', M, N, MB, NB,
     $                     NRHS, NPROW, NPCOL, CTIME( 1 ), CTIME( 2 ),
     $                     REAL( MAXERR ), REAL( MAXRESID )
                     ENDIF
                  ENDIF
*
*   Clean up.
*
   90             CONTINUE
                  CALL LACLOSE( IODEVAOLD, 'Nokeep', MYID, NPROCS,
     $                          INFO1 )
                  CALL LACLOSE( IODEV, 'Nokeep', MYID, NPROCS, INFO2 )
                  IERR( 1 ) = 0
                  IF( ( INFO1.NE.0 ) .OR. ( INFO2.NE.0 ) ) THEN
                     IERR( 1 ) = 1
                  ENDIF
                  CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1,
     $                          0 )
                  IF( IERR( 1 ).NE.0 ) THEN
                     IF( IAM.EQ.0 ) THEN
                        WRITE( NOUT, FMT = * )
     $                     'laclose returns info1,info2', INFO1, INFO2
                     ENDIF
                     KSKIP = KSKIP + 1
                     GOTO 100
                  ENDIF
                  IF( IAM.EQ.0 ) THEN
                  ENDIF
  100          CONTINUE
  110          CONTINUE
* end do el
  120       CONTINUE
  130       CONTINUE
* end do k
  140    CONTINUE
  150    CONTINUE
* end do j
  160 CONTINUE
  170 CONTINUE
* end do i
      CALL BLACS_GRIDEXIT( ICTXT )
      IF( IAM.EQ.0 ) THEN
         WRITE( NOUT, FMT = 9996 )KSKIP
 9996    FORMAT( I5, ' tests skipped because of illegal input values.' )
      ENDIF
      IF( ( NOUT.NE.0 ) .AND. ( NOUT.NE.6 ) ) THEN
         CLOSE ( NOUT )
      ENDIF
      CALL BLACS_EXIT( 0 )
      STOP '** all done **'
 9995 FORMAT( A4, 1X, I5, 1X, I5, 1X, I3, 1X, I3, 1X, I4, 1X, I4, 1X,
     $      I4, 1X, F8.2, 1X, F8.2, 1X, 1P, E8.2, 1X, 1P, E8.2 )
 9994 FORMAT( 'ILLEGAL ', A6, ': ', A5, ' = ', I3,
     $      '; It should be at least 1' )
 9993 FORMAT( 'ILLEGAL GRID: nprow*npcol = ', I4, '. It can be at most',
     $      I4 )
 9992 FORMAT( 'Bad ', A6, ' parameters: going on to next test case.' )
      END
