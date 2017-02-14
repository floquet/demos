      SUBROUTINE PFDESCINIT( DESCA, M, N, MB, NB, RSRC, CSRC, ICTXT,
     $                       IODEV, FILETYPE, MMB, NNB, ASIZE, FILENAME,
     $                       INFO )
*
*
*  -- ScaLAPACK auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     Oct 10, 1996
*
*
*  Purpose:
*  =======
*
*  Initialize out-of-core array descriptor.
*
* Notes
* =====
*
* I/O is performed in fixed size records. Each record can
* be considered an MMB x NNB distributed ScaLAPACK array.
*
* MMB should be a multiple of MB*nprow and
* NNB should be a multiple of NB*npcol.
*
* where  ICTXT is associated with an nprow x npcol processor grid.
*
* Arguments
* =========
*
* DESC    (output) INTEGER array of dimension DLEN_.
*         The array descriptor of a distributed matrix to be set.
*
* M       (global input) INTEGER
*         The number of rows in the distributed matrix. M >= 0.
*
* N       (global input) INTEGER
*         The number of columns in the distributed matrix. N >= 0.
*
* MB      (global input) INTEGER
*         The blocking factor used to distribute the rows of the
*         matrix. MB >= 1.
*
* NB      (global input) INTEGER
*         The blocking factor used to distribute the columns of the
*         matrix. NB >= 1.
*
* IRSRC   (global input) INTEGER
*         The process row over which the first row of the matrix is
*         distributed. 0 <= IRSRC < NPROW.
*
* ICSRC   (global input) INTEGER
*         The process column over which the first column of the
*         matrix is distributed. 0 <= ICSRC < NPCOL.
*
* ICTXT   (global input) INTEGER
*         The BLACS context handle, indicating the global context of
*         the operation on the matrix. The context itself is global.
*
* LLD     (local input)  INTEGER
*         The leading dimension of the local array storing the local
*         blocks of the distributed matrix. LLD >= MAX(1,LOCp(M)).
*
*
*
* IODEV   (global input) INTEGER
*        The device number (between 1..99) to be associated with
*         the out-of-core matrix.
*
* FILETYPE (global input) CHARACTER
*
*        = 'D': distributed files, each processor uses  a uniqe file
*
*        = 'S': A 'shared' file is used by all processors.
*               The file system must support concurrent read/write
*               operations, similar to M_ASYNC mode on the Paragon.
*               Data associated with each processor is stored in
*               contiguous blocks.
*
*               *NOTE* Most NFS implmentations will NOT support
*               concurrent read/write operations on shared files.
*
*        = 'I': A 'shared' file is used (similar to 'S') but
*               the data from each processor is 'interlaced'.
*
*
* MMB   (global input) INTEGER
*
*       The number of rows in I/O record.
*
* NNB   (global input) INTEGER
*
*       The number of columns in I/O record.
*
* ASIZE (local input) INTEGER
*
*       Size of local buffer associated with out-of-core matrix.
*
* FILENAME (global input) CHARACTER*(*)
*
*       Name of file associated with the out-of-core matrix.
*
*       An absolute path is used  if FILENAME starts with a '/'
*       (FILENAME(1:1) .eq. '/')
*
*       Otherwise, the file will reside on '/tmp' or '/pfs' directory.
*
*
* INFO    (output) INTEGER
*         = 0: successful exit
*         < 0: if INFO = -i, the i-th argument had an illegal value
*
*
*
*     .. Parameters ..
      INTEGER            DT_
      PARAMETER          ( DT_ = 1 )
      INTEGER            MB_, NB_
      PARAMETER          ( MB_ = 5, NB_ = 6 )
      INTEGER            RSRC_, CSRC_
      PARAMETER          ( RSRC_ = 7, CSRC_ = 8 )
      INTEGER            IODEV_, SIZE_, DISK_BLOCK_CYCLIC_2D
      PARAMETER          ( IODEV_ = 10, SIZE_ = 11,
     $                   DISK_BLOCK_CYCLIC_2D = 601 )
      INTEGER            FLEN_
      PARAMETER          ( FLEN_ = 12 )
*     ..
*     .. Scalar Arguments ..
      CHARACTER          FILETYPE
      CHARACTER*( * )    FILENAME
      INTEGER            ASIZE, CSRC, ICTXT, INFO, IODEV, M, MB, MMB, N,
     $                   NB, NNB, RSRC
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( FLEN_ )
*     ..
*     .. Local Scalars ..
      LOGICAL            ASIZEQUERY, ISVALID
      INTEGER            LDA, LINFO, LMMB, LNNB, MYID, MYPCOL, MYPROW,
     $                   NPCOL, NPROC, NPROW
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            NUMROC
      EXTERNAL           LSAME, NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           ASSERT, BLACS_GRIDINFO, BLACS_PINFO, DESCINIT,
     $                   LAOPEN, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MOD
*     ..
*     .. Executable Statements ..
      INFO = 0
      CALL BLACS_PINFO( MYID, NPROC )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYPROW, MYPCOL )
      LDA = NUMROC( M, MB, MYPROW, RSRC, NPROW )
      CALL DESCINIT( DESCA, M, N, MB, NB, CSRC, RSRC, ICTXT,
     $               MAX( 1, LDA ), INFO )
      IF( INFO.NE.0 ) THEN
         RETURN
      ENDIF
      IF( MOD( MMB, NPROW*MB ).NE.0 ) THEN
         INFO = -11
      ENDIF
      IF( MOD( NNB, NPCOL*NB ).NE.0 ) THEN
         INFO = -12
      ENDIF
      ISVALID = LSAME( FILETYPE, 'I' ) .OR. LSAME( FILETYPE, 'D' ) .OR.
     $          LSAME( FILETYPE, 'S' )
      IF( .NOT.ISVALID ) THEN
         INFO = -10
         CALL PXERBLA( ICTXT, 'PFDESCINIT', ( -INFO ) )
         RETURN
      ENDIF
      LINFO = 0
      LMMB = NUMROC( MMB, DESCA( MB_ ), MYPROW, DESCA( RSRC_ ), NPROW )
      LNNB = NUMROC( NNB, DESCA( NB_ ), MYPCOL, DESCA( CSRC_ ), NPCOL )
      CALL LAOPEN( IODEV, FILETYPE, FILENAME, DESCA, LMMB, LNNB, MMB,
     $             NNB, MYID, NPROC, LINFO )
      IF( LINFO.EQ.-3 ) THEN
         INFO = -14
      ELSE
         IF( LINFO.EQ.-1 ) THEN
            INFO = -9
         ENDIF
      ENDIF
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PFDESCINIT', ( -INFO ) )
         RETURN
      ENDIF
      CALL ASSERT( INFO.EQ.0, '** PFDESCINIT: laopen returns info ',
     $             INFO )
      DESCA( DT_ ) = DISK_BLOCK_CYCLIC_2D
      DESCA( IODEV_ ) = IODEV
      DESCA( SIZE_ ) = ASIZE
      ASIZEQUERY = ( ASIZE.EQ.-1 )
      IF( ( ASIZE.LE.0 ) .AND. ( .NOT.ASIZEQUERY ) ) THEN
         INFO = -13
         CALL PXERBLA( ICTXT, 'PFDESCINIT', ( -INFO ) )
      ENDIF
      RETURN
      END
