      SUBROUTINE PCBMATGEN( ICTXT, AFORM, DIAG, BWL, BWU, N, MB, NB, A,
     $                      LDA, IAROW, IACOL, ISEED, MYROW, MYCOL,
     $                      NPROW, NPCOL )
*
*  -- ScaLAPACK routine (version 1.2 ALPHA) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 10, 1996
*
*
*  Purpose
*  =======
*
*  PCBMATGEN : Parallel Complex Single precision Band MATrix GENerator.
*  (Re)Generate a distributed Band matrix A (or sub-matrix of A).
*
*  Arguments
*  =========
*
*  ICTXT   (global input) INTEGER
*          The BLACS context handle, indicating the global context of
*          the operation. The context itself is global.
*
*  AFORM   (global input) CHARACTER*1
*          if AFORM = 'L' : A is returned as a hermitian lower
*            triangular matrix.
*          if AFORM = 'U' : A is returned as a hermitian upper
*            triangular matrix.
*          if AFORM = 'T' : A is returned as the transpose of a general
*             matrix.
*          if AFORM = 'G' : A is returned as a general matrix.
*          otherwise a random matrix is generated.
*
*  DIAG    (global input) CHARACTER*1
*          if DIAG = 'D' : A is diagonally dominant.
*
*  M       (global input) INTEGER
*          The number of nonzero rows in the generated distributed
*           band matrix.
*
*  N       (global input) INTEGER
*          The number of columns in the generated distributed
*          matrix.
*
*  MB      (global input) INTEGER
*          The row blocking factor of the distributed matrix A.
*
*  NB      (global input) INTEGER
*          The column blocking factor of the distributed matrix A.
*
*  A       (local output) COMPLEX, pointer into the local memory to
*          an array of dimension ( LDA, * ) containing the local
*          pieces of the distributed matrix.
*
*  LDA     (local input) INTEGER
*          The leading dimension of the array containing the local
*          pieces of the distributed matrix A.
*
*  IAROW   (global input) INTEGER
*          The row processor coordinate which holds the first block
*          of the distributed matrix A.
*            A( DIAG_INDEX, I ) = A( DIAG_INDEX, I ) + BWL+BWU
*
*  IACOL   (global input) INTEGER
*          The column processor coordinate which holds the first
*          block of the distributed matrix A.
*
*  ISEED   (global input) INTEGER
*          The seed number to generate the distributed matrix A.
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
*  Notes
*  =====
*
*  This code is a simple wrapper around PCMATGEN, for band matrices.
*
*  =====================================================================
*
*  Implemented by: Andrew J. Cleary, May 10, 1996
*
*  =====================================================================
*
*     .. Scalar Arguments ..
      CHARACTER          AFORM, DIAG
      INTEGER            BWL, BWU, IACOL, IAROW, ICTXT, ISEED, LDA, MB,
     $                   MYCOL, MYROW, N, NB, NPCOL, NPROW
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * )
*     ..
*     .. Local Scalars ..
      INTEGER            DIAG_INDEX, I, M_MATGEN, NQ, N_MATGEN,
     $                   START_INDEX
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            NUMROC
      EXTERNAL           LSAME, NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           PCMATGEN
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, REAL
*     ..
*     .. Executable Statements ..
*
*
      IF( LSAME( AFORM, 'L' ) .OR. LSAME( AFORM, 'U' ) ) THEN
         M_MATGEN = BWL + 1
         N_MATGEN = N
         START_INDEX = 1
         IF( LSAME( AFORM, 'L' ) ) THEN
            DIAG_INDEX = 1
         ELSE
            DIAG_INDEX = BWL + 1
         ENDIF
      ELSE
         IF( LSAME( DIAG, 'D' ) ) THEN
            M_MATGEN = BWL + BWU + 1
            N_MATGEN = N
            DIAG_INDEX = BWU + 1
            START_INDEX = 1
         ELSE
            M_MATGEN = BWL + BWU + 1
            N_MATGEN = N
            DIAG_INDEX = BWL + BWU + 1
            START_INDEX = BWL + 1
         ENDIF
      ENDIF
*
      NQ = NUMROC( N, NB, MYCOL, IACOL, NPCOL )
*
*
*     Generate a random matrix initially
*
      IF( LSAME( AFORM, 'T' ) ) THEN
*
         CALL PCMATGEN( ICTXT, 'T', 'N', N_MATGEN, M_MATGEN, NB,
     $                  M_MATGEN, A( START_INDEX, 1 ), NB+10, IAROW,
     $                  IACOL, ISEED, 0, NQ, 0, M_MATGEN, MYCOL, MYROW,
     $                  NPCOL, NPROW )
*
      ELSE
*
         CALL PCMATGEN( ICTXT, 'N', 'N', M_MATGEN, N_MATGEN, M_MATGEN,
     $                  NB, A( START_INDEX, 1 ), LDA, IAROW, IACOL,
     $                  ISEED, 0, M_MATGEN, 0, NQ, MYROW, MYCOL, NPROW,
     $                  NPCOL )
*
      ENDIF
*
*
      IF( LSAME( DIAG, 'D' ) ) THEN
*
*       Loop over diagonal elements stored on this processor.
*
*
         DO 10 I = 1, NQ
            IF( LSAME( AFORM, 'T' ) ) THEN
               IF( NPROW.EQ.1 ) THEN
                  A( I, DIAG_INDEX ) = CMPLX( REAL( A( I,
     $                                 DIAG_INDEX ) )+2*( BWL+BWU+1 ) )
               ENDIF
            ELSE
               IF( NPROW.EQ.1 ) THEN
                  A( DIAG_INDEX, I ) = CMPLX( REAL( A( DIAG_INDEX,
     $                                 I ) )+2*( BWL+BWU+1 ) )
               ENDIF
            ENDIF
   10    CONTINUE
*
*
      ENDIF
*
      RETURN
*
*     End of PCBMATGEN
*
      END
