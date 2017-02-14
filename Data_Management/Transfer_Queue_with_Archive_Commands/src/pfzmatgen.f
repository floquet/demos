      SUBROUTINE PFZMATGEN( IODEV, ICTXT, AFORM, DIAG, M, N, MB, NB, A,
     $                      ASIZE, IAROW, IACOL, ISEED, MYPROW, MYPCOL,
     $                      NPROW, NPCOL, INFO )
*
*
*  -- ScaLAPACK auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     Oct 10, 1996
*
*
*
*  Purpose
*  =======
*
*  PFxMATGEN : Parallel Complex Double precision MATrix GENerator.
*  Generate (or regenerate) a distributed matrix A (or sub-matrix of A).
*
*  Arguments
*  =========
*
*  ICTXT   (global input) INTEGER
*          The BLACS context handle, indicating the global context of
*          the operation. The context itself is global.
*
*  AFORM   (global input) CHARACTER*1
*          if AFORM = 'S' : A is returned is a symmetric matrix.
*          if AFORM = 'H' : A is returned is a Hermitian matrix.
*          if AFORM = 'T' : A is overwritten with the transpose of
*                           what would normally be generated.
*          if AFORM = 'C' : A is overwritten with the conjugate trans-
*                           pose of what would normally be generated.
*          otherwise a random matrix is generated.
*
*  DIAG    (global input) CHARACTER*1
*          if DIAG = 'D' : A is diagonally dominant.
*
*  M       (global input) INTEGER
*          The number of rows in the generated distributed matrix.
*
*  N       (global input) INTEGER
*          The number of columns in the generated distributed
*          matrix.
*
*  MB      (global input) INTEGER
*          The row blocking factor of the distributed matrix A.
*
*  NB      (global input) INTEGER
*
*  A       (local output) COMPLEX*16, pointer into the local memory
*          to an array of dimension ( LDA, * ) containing the local
*          pieces of the distributed matrix.
*
*  ASIZE   (local input) INTEGER
*          The amount of local work space available.
*
*  IAROW   (global input) INTEGER
*          The row processor coordinate which holds the first block
*          of the distributed matrix A.
*
*  IACOL   (global input) INTEGER
*          The column processor coordinate which holds the first
*          block of the distributed matrix A.
*
*  ISEED   (global input) INTEGER
*          The seed number to generate the distributed matrix A.
*
*  IROFF   (local input) INTEGER
*          The number of local rows of A that have already been
*          generated.  It should be a multiple of MB.
*
*  IRNUM   (local input) INTEGER
*          The number of local rows to be generated.
*
*  ICOFF   (local input) INTEGER
*          The number of local columns of A that have already been
*          generated.  It should be a multiple of NB.
*
*  ICNUM   (local input) INTEGER
*          The number of local columns to be generated.
*
*  MYROW   (local input) INTEGER
*          The row process coordinate of the calling process.
*
*  MYCOL   (local input) INTEGER
*          The column process coordinate of the calling process.
*
*  NPROW   (global input) INTEGER
*          The number of process rows in the grid.
*
*  NPCOL   (global input) INTEGER
*          The number of process columns in the grid.
*
*
*  Notes
*  =====
*  The code internally uses  PxMATGEN.
*
*  See documentation on PxMATGEN for more details.
*
*     .. Parameters ..
      INTEGER            DLEN_
      PARAMETER          ( DLEN_ = 9 )
      INTEGER            LLD_
      PARAMETER          ( LLD_ = 9 )
*     ..
*     .. Scalar Arguments ..
      CHARACTER          AFORM, DIAG
      INTEGER            ASIZE, IACOL, IAROW, ICTXT, INFO, IODEV, ISEED,
     $                   M, MB, MYPCOL, MYPROW, N, NB, NPCOL, NPROW
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( * )
*     ..
*     .. Local Scalars ..
      LOGICAL            ASIZEQUERY, HASWORK, ISDIAG, ISMINE, ISUPPER,
     $                   ISVALID, NEEDSQUARE
      CHARACTER          LAFORM, LDIAG
      INTEGER            ANEED, CSRC, CSRCIO, I, IA, ICNUM, ICOFF, ICOL,
     $                   ICONTXT, IDX, IEND, IERR, II, IIHI, IILO,
     $                   IISIZE, INC, IRNUM, IROFF, IROW, ISIZE, ISTART,
     $                   JA, JEND, JJ, JJHI, JJLO, JJSIZE, JSIZE,
     $                   JSTART, LCINDX, LDA, LOCP, LOCQ, LRINDX, LSEED,
     $                   MM, MMB, NN, NNB, P0, Q0, RSRC, RSRCIO
      COMPLEX*16         ALPHA
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA( DLEN_ )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL, ILCM, INDXG2P, NUMROC
      EXTERNAL           LSAME, ICEIL, ILCM, INDXG2P, NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           ASSERT, BLACS_GRIDINFO, DESCINIT, IGSUM2D,
     $                   INFOG2L, LAIO_INFO, PFMAX2SIZE, PXERBLA,
     $                   PZMATGEN, ZLAWRITE
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, INT, MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
      INFO = 0
      RSRC = IAROW
      CSRC = IACOL
      ASIZEQUERY = ( ASIZE.EQ.-1 )
      CALL LAIO_INFO( IODEV, MM, NN, MMB, NNB, MB, NB, CSRCIO, RSRCIO,
     $                ICONTXT )
      CALL BLACS_GRIDINFO( ICONTXT, NPROW, NPCOL, MYPROW, MYPCOL )
      IF( ASIZEQUERY ) THEN
*
*         Make recommendation for amount of storage required.
*         May not be minimal requirement since this may lead
*          to unacceptable performance.
*
         NEEDSQUARE = ( LSAME( AFORM, 'S' ) .OR. LSAME( AFORM, 'H' ) )
         IF( NEEDSQUARE ) THEN
            MM = ILCM( MMB, NNB )
            NN = MM
         ELSE
            MM = MMB
            NN = NNB
         ENDIF
         LOCP = NUMROC( MM, MB, MYPROW, RSRC, NPROW )
         LOCQ = NUMROC( NN, NB, MYPCOL, CSRC, NPCOL )
         ANEED = MAX( LOCP, 1 )*MAX( LOCQ, 1 )
         A( 1 ) = ANEED
         INFO = -10
         RETURN
      ENDIF
*
* Use square patches to correctly handle
* symmetric diagonal blocks.
*
      IF( LSAME( AFORM, 'S' ) .OR. LSAME( AFORM, 'H' ) ) THEN
*
*       need mb == nb, square patches
*
         IF( MB.NE.NB ) THEN
            INFO = -7
            CALL PXERBLA( ICTXT, 'PFZMATGEN', 7 )
            RETURN
         ENDIF
         MM = ILCM( MMB, NNB )
         NN = MM
*
*       Try to reduce storage requirements.
*
         DO 10 I = 1, MAX( 1, INT( MM / MB ) )
            NN = MM
            IF( MM.LE.MB ) THEN
               GOTO 20
            ENDIF
            LOCP = NUMROC( MM, MB, MYPROW, RSRC, NPROW )
            LOCQ = NUMROC( NN, NB, MYPCOL, CSRC, NPCOL )
            ANEED = LOCP*LOCQ
            IERR = 0
            IF( ANEED.GT.ASIZE ) THEN
               IERR = 1
            ENDIF
            CALL IGSUM2D( ICTXT, 'A', ' ', 1, 1, IERR, 1, -1, -1 )
            IF( IERR.EQ.0 ) THEN
*
*               All processor agree they have sufficient storage.
*
               GOTO 20
            ENDIF
            MM = MM - MB
   10    CONTINUE
   20    CONTINUE
         MM = MAX( MM, MB )
         NN = MM
      ELSE
         CALL PFMAX2SIZE( ASIZE, M, N, MM, NN, IODEV )
      ENDIF
      IF( MM.NE.M ) THEN
         IF( MM.GE.MMB ) THEN
            MM = MM - MOD( MM, MMB )
         ENDIF
      ENDIF
      IF( NN.NE.N ) THEN
         IF( NN.GE.NNB ) THEN
            NN = NN - MOD( NN, NNB )
         ENDIF
      ENDIF
      IA = 1
      JA = 1
*
*   Align for local I/O.
*
      P0 = INDXG2P( IA, MB, MYPROW, RSRCIO, NPROW )
      Q0 = INDXG2P( JA, NB, MYPCOL, CSRCIO, NPCOL )
*
* Check if sufficient storage
*
      LOCP = NUMROC( MM, MB, MYPROW, P0, NPROW )
      LOCQ = NUMROC( NN, NB, MYPCOL, Q0, NPCOL )
      ANEED = LOCP*LOCQ
      IF( ASIZE.LT.ANEED ) THEN
         A( 1 ) = ANEED
         INFO = -10
         IF( .NOT.ASIZEQUERY ) THEN
            CALL PXERBLA( ICTXT, 'PFMATGEN', ( -INFO ) )
         ENDIF
         RETURN
      ENDIF
      LDA = MAX( 1, LOCP )
      INFO = 0
      CALL DESCINIT( DESCA, MM, NN, MB, NB, P0, Q0, ICTXT, LDA, INFO )
      CALL ASSERT( INFO.EQ.0, '** PFMATGEN : descinit returns info ',
     $             INFO )
      LSEED = ISEED
      IROFF = 0
      ICOFF = 0
      JJLO = MAX( 0, INT( JA / NN ) )
      JJHI = ICEIL( JA+N-1, NN )
      IILO = MAX( 0, INT( IA / MM ) )
      IIHI = ICEIL( IA+M-1, MM )
      IISIZE = IIHI - IILO + 1
      JJSIZE = JJHI - JJLO + 1
      JJ = JJLO
   30 CONTINUE
      IF( JJ.LE.JJHI ) THEN
         JSTART = MAX( JA, JJ*NN+1 )
         JEND = MIN( JA+( N-1 ), ( JJ+1 )*NN )
         JSIZE = JEND - JSTART + 1
         IF( JSIZE.LE.0 ) THEN
            GOTO 80
         ENDIF
         II = IILO
   40    CONTINUE
         IF( II.LE.IIHI ) THEN
            ISTART = MAX( IA, II*MM+1 )
            IEND = MIN( IA+( M-1 ), ( II+1 )*MM )
            ISIZE = IEND - ISTART + 1
            IF( ISIZE.LE.0 ) THEN
               GOTO 70
            ENDIF
            ISDIAG = ( ISTART.EQ.JSTART )
            IF( ISDIAG ) THEN
               LAFORM = AFORM
               LDIAG = DIAG
            ELSE
               LAFORM = 'N'
               LDIAG = 'N'
            ENDIF
            ISUPPER = ( JSTART.GT.ISTART )
            LAFORM = AFORM
            IF( LSAME( AFORM, 'S' ) .OR. LSAME( AFORM, 'H' ) ) THEN
               LSEED = ISEED + ( MIN( II, JJ )-1 ) +
     $                 ( MAX( II, JJ )-1 )*MAX( JJSIZE, IISIZE )
            ELSE
               LSEED = ISEED + ( II-1 ) + ( JJ-1 )*MAX( JJSIZE, IISIZE )
            ENDIF
            IF( LSAME( AFORM, 'S' ) ) THEN
               IF( ISDIAG ) THEN
                  LAFORM = AFORM
               ELSE
                  IF( ISUPPER ) THEN
                     LAFORM = 'T'
                  ELSE
                     LAFORM = 'N'
                  ENDIF
               ENDIF
            ENDIF
            IF( LSAME( AFORM, 'H' ) ) THEN
               IF( ISDIAG ) THEN
                  LAFORM = AFORM
               ELSE
                  IF( ISUPPER ) THEN
                     LAFORM = 'C'
                  ELSE
                     LAFORM = 'N'
                  ENDIF
               ENDIF
            ENDIF
            IRNUM = NUMROC( ISIZE, MB, MYPROW, RSRC, NPROW )
            ICNUM = NUMROC( JSIZE, NB, MYPCOL, CSRC, NPCOL )
            HASWORK = ( ISIZE.GE.1 ) .AND. ( JSIZE.GE.1 )
            IF( HASWORK ) THEN
               CALL PZMATGEN( ICTXT, LAFORM, LDIAG, ISIZE, JSIZE, MB,
     $                        NB, A, LDA, RSRC, CSRC, LSEED, IROFF,
     $                        IRNUM, ICOFF, ICNUM, MYPROW, MYPCOL,
     $                        NPROW, NPCOL )
               IF( LSAME( DIAG, 'D' ) .AND. ISDIAG ) THEN
*
*                       Scale diagonal entry  to make sure
*                       matrix is positive definite.
*
                  ALPHA = DCMPLX( DBLE( MAX( M, N ) ) )
                  INC = 1
                  DO 50 I = 1, MIN( ISIZE, JSIZE )
                     CALL INFOG2L( I, I, DESCA, NPROW, NPCOL, MYPROW,
     $                             MYPCOL, LRINDX, LCINDX, IROW, ICOL )
                     ISMINE = ( IROW.EQ.MYPROW ) .AND.
     $                        ( ICOL.EQ.MYPCOL )
                     IF( ISMINE ) THEN
                        IDX = LRINDX + ( LCINDX-1 )*DESCA( LLD_ )
                        ISVALID = ( 1.LE.IDX ) .AND. ( IDX.LE.ASIZE )
                        IF( .NOT.ISVALID ) THEN
                           WRITE( *, FMT = 9999 )MYPROW, MYPCOL, IDX,
     $                        I, LRINDX, LCINDX, IROW, ICOL
 9999                      FORMAT( 'PFMATGEN: ', 'myprow,mypcol ',
     $                           2( 1X, I4 ), ' invalid idx ', I7,
     $                           'i,lrindx,lcindx ', 3( 1X, I6 ),
     $                           'irow,icol ', 2( 1X, I4 ) )
                           STOP '** error ** '
                        ENDIF
                        A( IDX ) = ABS( A( IDX )*ALPHA )
                     ENDIF
   50             CONTINUE
   60             CONTINUE
               ENDIF
               CALL ZLAWRITE( IODEV, ISIZE, JSIZE, ISTART, JSTART, A, 1,
     $                        1, DESCA, INFO )
               CALL ASSERT( INFO.EQ.0,
     $                      '** PFMATGEN: LAWRITE returns info ', INFO )
            ENDIF
            II = II + 1
            GOTO 40
         ENDIF
   70    CONTINUE
         JJ = JJ + 1
         GOTO 30
      ENDIF
   80 CONTINUE
      A( 1 ) = ANEED
      INFO = 0
      RETURN
      END
