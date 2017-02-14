      SUBROUTINE DRIVERINFO( SUMMARY, NOUT, NFACT, FACTOR, LDFACT, NMAT,
     $                       MVAL, LDMVAL, NVAL, LDNVAL, NRHSVAL,
     $                       LDNRHSVAL, ASIZEVAL, LDASIZEVAL, NNB,
     $                       MBVAL, LDMBVAL, NBVAL, LDNBVAL, NGRIDS,
     $                       PVAL, LDPVAL, QVAL, LDQVAL, IAM, NPROCS )
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
* Get startup information for test driver.
*
*
*     .. Parameters ..
      INTEGER            NIN
      PARAMETER          ( NIN = 11 )
*     ..
*     .. Scalar Arguments ..
      CHARACTER*( * )    SUMMARY
      INTEGER            IAM, LDASIZEVAL, LDFACT, LDMBVAL, LDMVAL,
     $                   LDNBVAL, LDNRHSVAL, LDNVAL, LDPVAL, LDQVAL,
     $                   NFACT, NGRIDS, NMAT, NNB, NOUT, NPROCS
*     ..
*     .. Array Arguments ..
      CHARACTER*2        FACTOR( LDFACT )
      INTEGER            ASIZEVAL( LDASIZEVAL ), MBVAL( LDMBVAL ),
     $                   MVAL( LDMVAL ), NBVAL( LDNBVAL ),
     $                   NRHSVAL( LDNRHSVAL ), NVAL( LDNVAL ),
     $                   PVAL( LDPVAL ), QVAL( LDQVAL )
*     ..
*     .. Local Scalars ..
      LOGICAL            ISVALID
      CHARACTER          CH
      CHARACTER*2        FACT
      CHARACTER*10       ERRMSG
      CHARACTER*80       FILENAME
      INTEGER            I, ICTXT, II, IVAL, LINENO, MYID
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_ABORT, BLACS_GET, BLACS_GRIDEXIT,
     $                   BLACS_GRIDINIT, BLACS_PINFO, BLACS_SETUP,
     $                   IGEBR2D, IGEBS2D
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, LEN, MAX
*     ..
*     .. Executable Statements ..
      IF( IAM.EQ.0 ) THEN
         LINENO = 0
*
*       Open file and skip data file header.
*
         FILENAME = 'testdriver.in'
         OPEN( NIN, FILE = FILENAME, STATUS = 'old', FORM = 'formatted',
     $       ACCESS = 'sequential', ERR = 150 )
         REWIND ( NIN )
         LINENO = LINENO + 1
         READ( NIN, FMT = '(A)', ERR = 160 )SUMMARY
         SUMMARY = ' '
*
*       Read name and unit number for summary output file.
*
         LINENO = LINENO + 1
         READ( NIN, FMT = '(A)', ERR = 160 )SUMMARY
         LINENO = LINENO + 1
         READ( NIN, FMT = *, ERR = 160 )NOUT
         IF( ( NOUT.NE.0 ) .AND. ( NOUT.NE.6 ) ) THEN
            FILENAME = SUMMARY
            OPEN( NOUT, FILE = FILENAME, ACCESS = 'sequential',
     $          STATUS = 'unknown', FORM = 'formatted', ERR = 150 )
            REWIND ( NOUT )
         ENDIF
*
*       Read parameters for test.
*
         FILENAME = 'testdriver.in'
         LINENO = LINENO + 1
         READ( NIN, FMT = *, ERR = 160 )NFACT
         ISVALID = ( 1.LE.NFACT ) .AND. ( NFACT.LE.LDFACT )
         IF( .NOT.ISVALID ) THEN
            ERRMSG = 'NFACT'
            GOTO 170
         ENDIF
         DO 10 I = 1, NFACT
            LINENO = LINENO + 1
            READ( NIN, FMT = '(A)', ERR = 160 )FACT
            FACTOR( I ) = FACT
   10    CONTINUE
   20    CONTINUE
         LINENO = LINENO + 1
         READ( NIN, FMT = *, ERR = 160 )NMAT
         ISVALID = ( 1.LE.NMAT ) .AND. ( NMAT.LE.LDNVAL ) .AND.
     $             ( NMAT.LE.LDMVAL )
         IF( .NOT.ISVALID ) THEN
            ERRMSG = 'NMAT'
            GOTO 170
         ENDIF
*
* Get matrix and dimensions.
*
         LINENO = LINENO + 1
         READ( NIN, FMT = *, ERR = 160 )( MVAL( I ), I = 1, NMAT )
         LINENO = LINENO + 1
         READ( NIN, FMT = *, ERR = 160 )( NVAL( I ), I = 1, NMAT )
         LINENO = LINENO + 1
         READ( NIN, FMT = *, ERR = 160 )( NRHSVAL( I ), I = 1, NMAT )
         LINENO = LINENO + 1
         READ( NIN, FMT = *, ERR = 160 )( ASIZEVAL( I ), I = 1, NMAT )
*
* Get values of nb.
*
         LINENO = LINENO + 1
         READ( NIN, FMT = *, ERR = 160 )NNB
         ISVALID = ( 1.LE.NNB ) .AND. ( NNB.LE.LDMBVAL ) .AND.
     $             ( NNB.LE.LDNBVAL )
         IF( .NOT.ISVALID ) THEN
            ERRMSG = 'NNB'
            GOTO 170
         ENDIF
         LINENO = LINENO + 1
         READ( NIN, FMT = *, ERR = 160 )( MBVAL( I ), I = 1, NNB )
         LINENO = LINENO + 1
         READ( NIN, FMT = *, ERR = 160 )( NBVAL( I ), I = 1, NNB )
*
* Get number of grids.
*
         LINENO = LINENO + 1
         READ( NIN, FMT = *, ERR = 160 )NGRIDS
         ISVALID = ( 1.LE.NGRIDS ) .AND. ( NGRIDS.LE.LDPVAL ) .AND.
     $             ( NGRIDS.LE.LDQVAL )
         IF( .NOT.ISVALID ) THEN
            ERRMSG = 'NGRIDS'
            GOTO 170
         ENDIF
         LINENO = LINENO + 1
         READ( NIN, FMT = *, ERR = 160 )( PVAL( I ), I = 1, NGRIDS )
         LINENO = LINENO + 1
         READ( NIN, FMT = *, ERR = 160 )( QVAL( I ), I = 1, NGRIDS )
         CLOSE ( NIN )
*
**
**        Temporarily define blacs grid to include all processes so
**        all processes have needed startup information
**
*
         CALL BLACS_PINFO( MYID, NPROCS )
         IF( NPROCS.LT.1 ) THEN
            NPROCS = 0
            DO 30 I = 1, NGRIDS
               NPROCS = MAX( NPROCS, PVAL( I )*QVAL( I ) )
   30       CONTINUE
   40       CONTINUE
            CALL BLACS_SETUP( IAM, NPROCS )
         ENDIF
         CALL BLACS_GET( -1, 0, ICTXT )
         CALL BLACS_GRIDINIT( ICTXT, 'Row-major', 1, NPROCS )
      ELSE
*
**
**        If in pvm, must participate setting up virtual machine
**
*
         IF( NPROCS.LT.1 ) THEN
            CALL BLACS_SETUP( IAM, NPROCS )
         ENDIF
*
**
**        Temporarily define blacs grid to include all processes so
**        all processes have needed startup information
**
*
         CALL BLACS_GET( -1, 0, ICTXT )
         CALL BLACS_GRIDINIT( ICTXT, 'Row-major', 1, NPROCS )
      ENDIF
*
* Perform Broadcasts.
*
      IF( IAM.EQ.0 ) THEN
         CALL IGEBS2D( ICTXT, 'All', ' ', 1, 1, NOUT, 1 )
      ELSE
         CALL IGEBR2D( ICTXT, 'All', ' ', 1, 1, NOUT, 1, 0, 0 )
      ENDIF
      IF( IAM.EQ.0 ) THEN
         CALL IGEBS2D( ICTXT, 'All', ' ', 1, 1, NFACT, 1 )
      ELSE
         CALL IGEBR2D( ICTXT, 'All', ' ', 1, 1, NFACT, 1, 0, 0 )
      ENDIF
      DO 70 I = 1, NFACT
         FACT = FACTOR( I )
         DO 50 II = 1, LEN( FACT )
            IF( IAM.EQ.0 ) THEN
               CH = FACT( II: II )
               IVAL = ICHAR( CH )
            ENDIF
            IF( IAM.EQ.0 ) THEN
               CALL IGEBS2D( ICTXT, 'All', ' ', 1, 1, IVAL, 1 )
            ELSE
               CALL IGEBR2D( ICTXT, 'All', ' ', 1, 1, IVAL, 1, 0, 0 )
            ENDIF
            CH = CHAR( IVAL )
            FACT( II: II ) = CH
   50    CONTINUE
   60    CONTINUE
         FACTOR( I ) = FACT
   70 CONTINUE
   80 CONTINUE
      IF( IAM.EQ.0 ) THEN
         CALL IGEBS2D( ICTXT, 'All', ' ', 1, 1, NMAT, 1 )
      ELSE
         CALL IGEBR2D( ICTXT, 'All', ' ', 1, 1, NMAT, 1, 0, 0 )
      ENDIF
      DO 90 I = 1, NMAT
         IF( IAM.EQ.0 ) THEN
            CALL IGEBS2D( ICTXT, 'All', ' ', 1, 1, MVAL( I ), 1 )
         ELSE
            CALL IGEBR2D( ICTXT, 'All', ' ', 1, 1, MVAL( I ), 1, 0, 0 )
         ENDIF
         IF( IAM.EQ.0 ) THEN
            CALL IGEBS2D( ICTXT, 'All', ' ', 1, 1, NVAL( I ), 1 )
         ELSE
            CALL IGEBR2D( ICTXT, 'All', ' ', 1, 1, NVAL( I ), 1, 0, 0 )
         ENDIF
         IF( IAM.EQ.0 ) THEN
            CALL IGEBS2D( ICTXT, 'All', ' ', 1, 1, NRHSVAL( I ), 1 )
         ELSE
            CALL IGEBR2D( ICTXT, 'All', ' ', 1, 1, NRHSVAL( I ), 1, 0,
     $                    0 )
         ENDIF
         IF( IAM.EQ.0 ) THEN
            CALL IGEBS2D( ICTXT, 'All', ' ', 1, 1, ASIZEVAL( I ), 1 )
         ELSE
            CALL IGEBR2D( ICTXT, 'All', ' ', 1, 1, ASIZEVAL( I ), 1, 0,
     $                    0 )
         ENDIF
   90 CONTINUE
  100 CONTINUE
      IF( IAM.EQ.0 ) THEN
         CALL IGEBS2D( ICTXT, 'All', ' ', 1, 1, NNB, 1 )
      ELSE
         CALL IGEBR2D( ICTXT, 'All', ' ', 1, 1, NNB, 1, 0, 0 )
      ENDIF
      DO 110 I = 1, NNB
         IF( IAM.EQ.0 ) THEN
            CALL IGEBS2D( ICTXT, 'All', ' ', 1, 1, MBVAL( I ), 1 )
         ELSE
            CALL IGEBR2D( ICTXT, 'All', ' ', 1, 1, MBVAL( I ), 1, 0, 0 )
         ENDIF
         IF( IAM.EQ.0 ) THEN
            CALL IGEBS2D( ICTXT, 'All', ' ', 1, 1, NBVAL( I ), 1 )
         ELSE
            CALL IGEBR2D( ICTXT, 'All', ' ', 1, 1, NBVAL( I ), 1, 0, 0 )
         ENDIF
  110 CONTINUE
  120 CONTINUE
      IF( IAM.EQ.0 ) THEN
         CALL IGEBS2D( ICTXT, 'All', ' ', 1, 1, NGRIDS, 1 )
      ELSE
         CALL IGEBR2D( ICTXT, 'All', ' ', 1, 1, NGRIDS, 1, 0, 0 )
      ENDIF
      DO 130 I = 1, NGRIDS
         IF( IAM.EQ.0 ) THEN
            CALL IGEBS2D( ICTXT, 'All', ' ', 1, 1, PVAL( I ), 1 )
         ELSE
            CALL IGEBR2D( ICTXT, 'All', ' ', 1, 1, PVAL( I ), 1, 0, 0 )
         ENDIF
         IF( IAM.EQ.0 ) THEN
            CALL IGEBS2D( ICTXT, 'All', ' ', 1, 1, QVAL( I ), 1 )
         ELSE
            CALL IGEBR2D( ICTXT, 'All', ' ', 1, 1, QVAL( I ), 1, 0, 0 )
         ENDIF
  130 CONTINUE
  140 CONTINUE
      CALL BLACS_GRIDEXIT( ICTXT )
      RETURN
  150 CONTINUE
      WRITE( *, FMT = * )' error in opening file ', FILENAME
      CALL BLACS_ABORT( ICTXT, 1 )
      STOP '** error ** '
  160 CONTINUE
      WRITE( *, FMT = * )' error in reading file ', FILENAME
      WRITE( *, FMT = * )' around line ', LINENO
      STOP '** error ** '
  170 CONTINUE
      WRITE( *, FMT = * )' error with parameter ', ERRMSG
      STOP '** error ** '
      END
