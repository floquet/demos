       program picalc
C       use ISO_FORTRAN_ENV
C2345678
C       use IFPORT
CC Intel-posix includes
C       use IFPOSIX
       include "mpif.h"
       double precision PI25DT
       parameter (PI25DT = 3.1415926535898d0, n=1100)
       double precision mypi,pi,h,sum,x,f,a
       character*(MPI_MAX_PROCESSOR_NAME) PNAME
       integer NAMESIZE
C            functions to integrate
       f(a) = 4.d0 / (1.d0 + a*a)
       call MPI_INIT(ierr)
       Call MPI_COMM_RANK (MPI_COMM_WORLD, myid, ierr)
       call MPI_COMM_SIZE (MPI_COMM_WORLD, numprocs, ierr)
       call MPI_Get_processor_name(PNAME,NAMESIZE, IERR)
C              calculate the interval size
       do i=0,numprocs-1
         if (i.eq.myid) then
           write(6,*) 'numprocs=',numprocs,' myid=',myid,' ',
     &'processor=',PNAME(1:NAMESIZE)
C            ierr=commitqq(6)
c            ierr=commitqq(6)
C           call PXFTCDRAIN(6,ierr)
            FLUSH(6)
         END IF
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       end do
C
       h = 1.d0 / n
       sum = 0.0d0
       do i=myid+1,n,numprocs
         x = h * (dble(i) - 0.5d0)
	 sum = sum + f(x)
       end do
       mypi = h * sum
C              collect all the partial sums
       call MPI_REDUCE (mypi,pi,1,MPI_DOUBLE_PRECISION,
     $         MPI_SUM,0,MPI_COMM_WORLD,ierr)
       CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
       if (myid .eq. 0) then
         print *, 'pi is ', pi, 'Error is', dabs(pi - PI25DT)
c            ierr=commitqq(6)
C           call PXFTCDRAIN(6,ierr)
         call MPI_FINALIZE (ierr)
       else 
         call MPI_FINALIZE (ierr)
       endif
       stop
       end
