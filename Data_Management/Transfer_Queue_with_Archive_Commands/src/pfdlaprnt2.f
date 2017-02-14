      SUBROUTINE PFDLAPRNT2( M, N, IAOFF, JAOFF, A, IA, JA, DESCA,
     $                       IRPRNT, ICPRNT, CMATNM, NOUT, WORK )
*
*  -- ScaLAPACK routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     Oct 10, 1996
*
*
*
*
* Purpose
* =======
*
* PFxLAPRNT2 prints to the standard output a distributed matrix sub( a )
* denoting a(ia:ia+m-1,ja:ja+n-1). the local pieces are sent and
* printed by the process of coordinates (irprnt, icprnt).
*
* notes
* =====
*
* each global data object is described by an associated description
* vector.  this vector stores the information required to establish
* the mapping between an object element and its corresponding process
* and memory location.
*
* let a be a generic term for any 2d block cyclicly distributed array.
* such a global array has an associated description vector desca.
* in the following comments, the character _ should be read as
* "of the global array".
*
* notation        stored in      explanation
* --------------- -------------- --------------------------------------
* dtype_a(global) desca( dtype_ )the descriptor type.  in this case,
*                                dtype_a = 1.
* ctxt_a (global) desca( ctxt_ ) the blacs context handle, indicating
*                                the blacs process grid a is distribu-
*                                ted over. the context itself is glo-
*                                bal, but the handle (the integer
*                                value) may vary.
* m_a    (global) desca( m_ )    the number of rows in the global
*                                array a.
* n_a    (global) desca( n_ )    the number of columns in the global
*                                array a.
* mb_a   (global) desca( mb_ )   the blocking factor used to distribute
*                                the rows of the array.
* nb_a   (global) desca( nb_ )   the blocking factor used to distribute
*                                the columns of the array.
* rsrc_a (global) desca( rsrc_ ) the process row over which the first
*                                row of the array a is distributed.
* csrc_a (global) desca( csrc_ ) the process column over which the
*                                first column of the array a is
*                                distributed.
* lld_a  (local)  desca( lld_ )  the leading dimension of the local
*                                array.  lld_a >= max(1,locp(m_a)).
*
* let k be the number of rows or columns of a distributed matrix,
* and assume that its process grid has dimension p x q.
* locp( k ) denotes the number of elements of k that a process
* would receive if k were distributed over the p processes of its
* process column.
* similarly, locq( k ) denotes the number of elements of k that a
* process would receive if k were distributed over the q processes of
* its process row.
* the values of locp() and locq() may be determined via a call to the
* scalapack tool function, numroc:
*         locp( m ) = numroc( m, mb_a, myrow, rsrc_a, nprow ),
*         locq( n ) = numroc( n, nb_a, mycol, csrc_a, npcol ).
* an upper bound for these quantities may be computed by:
*         locp( m ) <= ceil( ceil(m/mb_a)/nprow )*mb_a
*         locq( n ) <= ceil( ceil(n/nb_a)/npcol )*nb_a
*
* arguments
* =========
*
* m       (global input) integer
*         the number of rows to be operated on i.e the number of rows
*         of the distributed submatrix sub( a ). m >= 0.
*
* n       (global input) integer
*         the number of columns to be operated on i.e the number of
*         columns of the distributed submatrix sub( a ). n >= 0.
*
* a       (local input) DTYPE pointer into the local memory to a
*         local array of dimension (lld_a, locq(ja+n-1) ) containing
*         the local pieces of the distributed matrix sub( a ).
*
* ia      (global input) integer
*         the row index in the global array a indicating the first
*         row of sub( a ).
*
* ja      (global input) integer
*         the column index in the global array a indicating the
*         first column of sub( a ).
*
* desca   (global and local input) integer array of dimension dlen_.
*         the array descriptor for the distributed matrix a.
*
* irprnt  (global input) integer
*         the row index of the printing process.
*
* icprnt  (global input) integer
*         the column index of the printing process.
*
* cmatnm  (global input) character*(*)
*         identifier of the distributed matrix to be printed.
*
* nout    (global input) integer
*         the unit number for output file. nout = 6, ouput to screen,
*         nout = 0, output to stderr.
*
* work    (local workspace)
*         working array of minimum size equal to mb_a.
*
* =====================================================================
*
*
*     .. Parameters ..
      INTEGER            CTXT_, MB_, NB_
      PARAMETER          ( CTXT_ = 2, MB_ = 5, NB_ = 6 )
      INTEGER            LLD_
      PARAMETER          ( LLD_ = 9 )
*     ..
*     .. Scalar Arguments ..
      CHARACTER*( * )    CMATNM
      INTEGER            IA, IAOFF, ICPRNT, IRPRNT, JA, JAOFF, M, N,
     $                   NOUT
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      DOUBLE PRECISION   A( * ), WORK( * )
*     ..
*     .. Local Scalars ..
      INTEGER            H, I, IACOL, IAROW, IB, ICTXT, ICURCOL,
     $                   ICURROW, II, IIA, IN, J, JB, JJ, JJA, JN, K,
     $                   LDA, MYCOL, MYROW, NPCOL, NPROW
*     ..
*     .. External Functions ..
      INTEGER            ICEIL
      EXTERNAL           ICEIL
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_BARRIER, BLACS_GRIDINFO, DGERV2D,
     $                   DGESD2D, INFOG2L
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MIN, MOD
*     ..
*     .. Executable Statements ..
*
*      get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, IIA, JJA,
     $              IAROW, IACOL )
      ICURROW = IAROW
      ICURCOL = IACOL
      II = IIA
      JJ = JJA
      LDA = DESCA( LLD_ )
*
*      handle the first block of column separately
*
      JN = MIN( ICEIL( JA, DESCA( NB_ ) )*DESCA( NB_ ), JA+N-1 )
      JB = JN - JA + 1
      DO 110 H = 0, JB - 1
         IN = MIN( ICEIL( IA, DESCA( MB_ ) )*DESCA( MB_ ), IA+M-1 )
         IB = IN - IA + 1
         IF( ICURROW.EQ.IRPRNT .AND. ICURCOL.EQ.ICPRNT ) THEN
            IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
               DO 10 K = 0, IB - 1
                  WRITE( NOUT, FMT = 9999 )CMATNM, IAOFF + IA + K,
     $               JAOFF + JA + H, DBLE( A( II+K+( JJ+H-1 )*LDA ) )
   10          CONTINUE
   20          CONTINUE
            ENDIF
         ELSE
            IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
               CALL DGESD2D( ICTXT, IB, 1, A( II+( JJ+H-1 )*LDA ), LDA,
     $                       IRPRNT, ICPRNT )
            ELSE
               IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
                  CALL DGERV2D( ICTXT, IB, 1, WORK, DESCA( MB_ ),
     $                          ICURROW, ICURCOL )
                  DO 30 K = 1, IB
                     WRITE( NOUT, FMT = 9999 )CMATNM,
     $                  IAOFF + IA + K - 1, JAOFF + JA + H,
     $                  DBLE( WORK( K ) )
   30             CONTINUE
   40             CONTINUE
               ENDIF
            ENDIF
         ENDIF
         IF( MYROW.EQ.ICURROW ) THEN
            II = II + IB
         ENDIF
         ICURROW = MOD( ICURROW+1, NPROW )
         CALL BLACS_BARRIER( ICTXT, 'all' )
*
*         loop over remaining block of rows
*
         DO 90 I = IN + 1, IA + M - 1, DESCA( MB_ )
            IB = MIN( DESCA( MB_ ), IA+M-I )
            IF( ICURROW.EQ.IRPRNT .AND. ICURCOL.EQ.ICPRNT ) THEN
               IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
                  DO 50 K = 0, IB - 1
                     WRITE( NOUT, FMT = 9999 )CMATNM, IAOFF + I + K,
     $                  JAOFF + JA + H, DBLE( A( II+K+( JJ+H-1 )*LDA ) )
   50             CONTINUE
   60             CONTINUE
               ENDIF
            ELSE
               IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
                  CALL DGESD2D( ICTXT, IB, 1, A( II+( JJ+H-1 )*LDA ),
     $                          LDA, IRPRNT, ICPRNT )
               ELSE
                  IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
                     CALL DGERV2D( ICTXT, IB, 1, WORK, DESCA( MB_ ),
     $                             ICURROW, ICURCOL )
                     DO 70 K = 1, IB
                        WRITE( NOUT, FMT = 9999 )CMATNM,
     $                     IAOFF + I + K - 1, JAOFF + JA + H,
     $                     DBLE( WORK( K ) )
   70                CONTINUE
   80                CONTINUE
                  ENDIF
               ENDIF
            ENDIF
            IF( MYROW.EQ.ICURROW ) THEN
               II = II + IB
            ENDIF
            ICURROW = MOD( ICURROW+1, NPROW )
            CALL BLACS_BARRIER( ICTXT, 'all' )
   90    CONTINUE
  100    CONTINUE
         II = IIA
         ICURROW = IAROW
  110 CONTINUE
  120 CONTINUE
      IF( MYCOL.EQ.ICURCOL ) THEN
         JJ = JJ + JB
      ENDIF
      ICURCOL = MOD( ICURCOL+1, NPCOL )
      CALL BLACS_BARRIER( ICTXT, 'all' )
*
*      loop over remaining column blocks
*
      DO 250 J = JN + 1, JA + N - 1, DESCA( NB_ )
         JB = MIN( DESCA( NB_ ), JA+N-J )
         DO 230 H = 0, JB - 1
            IN = MIN( ICEIL( IA, DESCA( MB_ ) )*DESCA( MB_ ), IA+M-1 )
            IB = IN - IA + 1
            IF( ICURROW.EQ.IRPRNT .AND. ICURCOL.EQ.ICPRNT ) THEN
               IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
                  DO 130 K = 0, IB - 1
                     WRITE( NOUT, FMT = 9999 )CMATNM, IAOFF + IA + K,
     $                  JAOFF + J + H, DBLE( A( II+K+( JJ+H-1 )*LDA ) )
  130             CONTINUE
  140             CONTINUE
               ENDIF
            ELSE
               IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
                  CALL DGESD2D( ICTXT, IB, 1, A( II+( JJ+H-1 )*LDA ),
     $                          LDA, IRPRNT, ICPRNT )
               ELSE
                  IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
                     CALL DGERV2D( ICTXT, IB, 1, WORK, DESCA( MB_ ),
     $                             ICURROW, ICURCOL )
                     DO 150 K = 1, IB
                        WRITE( NOUT, FMT = 9999 )CMATNM,
     $                     IAOFF + IA + K - 1, JAOFF + J + H,
     $                     DBLE( WORK( K ) )
  150                CONTINUE
  160                CONTINUE
                  ENDIF
               ENDIF
            ENDIF
            IF( MYROW.EQ.ICURROW ) THEN
               II = II + IB
            ENDIF
            ICURROW = MOD( ICURROW+1, NPROW )
            CALL BLACS_BARRIER( ICTXT, 'all' )
*
*            loop over remaining block of rows
*
            DO 210 I = IN + 1, IA + M - 1, DESCA( MB_ )
               IB = MIN( DESCA( MB_ ), IA+M-I )
               IF( ICURROW.EQ.IRPRNT .AND. ICURCOL.EQ.ICPRNT ) THEN
                  IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
                     DO 170 K = 0, IB - 1
                        WRITE( NOUT, FMT = 9999 )CMATNM, IAOFF + I + K,
     $                     JAOFF + J + H, DBLE( A( II+K+( JJ+H-1 )*
     $                     LDA ) )
  170                CONTINUE
  180                CONTINUE
                  ENDIF
               ELSE
                  IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
                     CALL DGESD2D( ICTXT, IB, 1, A( II+( JJ+H-1 )*LDA ),
     $                             LDA, IRPRNT, ICPRNT )
                  ELSE
                     IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
                        CALL DGERV2D( ICTXT, IB, 1, WORK, DESCA( MB_ ),
     $                                ICURROW, ICURCOL )
                        DO 190 K = 1, IB
                           WRITE( NOUT, FMT = 9999 )CMATNM,
     $                        IAOFF + I + K - 1, JAOFF + J + H,
     $                        DBLE( WORK( K ) )
  190                   CONTINUE
  200                   CONTINUE
                     ENDIF
                  ENDIF
               ENDIF
               IF( MYROW.EQ.ICURROW ) THEN
                  II = II + IB
               ENDIF
               ICURROW = MOD( ICURROW+1, NPROW )
               CALL BLACS_BARRIER( ICTXT, 'all' )
  210       CONTINUE
  220       CONTINUE
            II = IIA
            ICURROW = IAROW
  230    CONTINUE
  240    CONTINUE
         IF( MYCOL.EQ.ICURCOL ) THEN
            JJ = JJ + JB
         ENDIF
         ICURCOL = MOD( ICURCOL+1, NPCOL )
         CALL BLACS_BARRIER( ICTXT, 'all' )
  250 CONTINUE
  260 CONTINUE
 9999 FORMAT( A, '(', I6, ',', I6, ')=', D30.18 )
      RETURN
*
**     end of routine
*
      END
