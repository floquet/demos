      INTEGER          FUNCTION NUMROC2( N, INDXGLOB, NB, IPROC,
     $                 ISRCPROC, NPROCS )
*
*
*  -- ScaLAPACK auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     Oct 10, 1996
*
*
*  Purpose
*  =======
*
*  NUMROC2 computes the NUMber of Rows Or Columns of a distributed
*  matrix owned by the process indicated by IPROC starting at
*  global index indxglob.
*
*  Index range in   indxglb:(indxglb+n-1).
*
*     .. Scalar Arguments ..
      INTEGER            INDXGLOB, IPROC, ISRCPROC, N, NB, NPROCS
*     ..
*     .. Local Scalars ..
      INTEGER            ANS, ANS1, ANS2, NN
*     ..
*     .. External Functions ..
      INTEGER            NUMROC
      EXTERNAL           NUMROC
*     ..
*     .. Executable Statements ..
      IF( INDXGLOB.EQ.1 ) THEN
         ANS = NUMROC( N, NB, IPROC, ISRCPROC, NPROCS )
      ELSE
         NN = INDXGLOB - 1
         ANS1 = 0
         IF( NN.GE.1 ) THEN
            ANS1 = NUMROC( NN, NB, IPROC, ISRCPROC, NPROCS )
         ENDIF
         NN = INDXGLOB + N - 1
         ANS2 = 0
         IF( NN.GE.1 ) THEN
            ANS2 = NUMROC( NN, NB, IPROC, ISRCPROC, NPROCS )
         ENDIF
         ANS = ANS2 - ANS1
      ENDIF
      NUMROC2 = ( ANS )
      RETURN
      END
