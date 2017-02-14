      SUBROUTINE ASSERT( LCOND, MESG, IVAL )
*     .. Scalar Arguments ..
      LOGICAL            LCOND
      CHARACTER*( * )    MESG
      INTEGER            IVAL
*     ..
*     .. Executable Statements ..
      IF( .NOT.LCOND ) THEN
         WRITE( *, FMT = * )MESG, IVAL
         STOP '** assertion error ** '
      ENDIF
      RETURN
      END
