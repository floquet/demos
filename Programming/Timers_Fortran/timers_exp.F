      program test

      implicit none
#ifdef USEMPI
      include "mpif.h"
      integer ierr
#endif
      integer length, i, ntimes
      parameter (length = 90000000)
      real*8 a(length), sumof

!!! declarations for each of the Four Timers..
! for f77 intrinsic SECNDS
      real t1, dsecnds

! for f90/95 intrinsic SYSTEM_CLOCK
      integer*4 icount, irate
      real*8 t2, dsysclk

! for f90/95 intrinsic CPU_TIME
      real*4 t3, t4, dcputim

! for f95/2000+ intrinsic DATE_AND_TIME
      integer t7(8), t8(8)
      real*4 dandttim
      character*8 cdate
      character*10 ctime
      character*5 czone

#ifdef USEMPI
! for MPI timer MPI_WTIME()
      real*8 t5, t6, mpitim
#endif

!! start  SECNDS  Timing
      t1 = secnds(0.0)

!! start  SYSTEM_CLOCK  Timing
      call system_clock (icount, irate)
      t2 = dble(icount) / dble(irate)

!! start  CPU_TIME  Timing
      call cpu_time(t3)

#ifdef USEMPI
!! start  MPI_WTIME  Timing
      call mpi_init(ierr)
      t5 = mpi_wtime()
#endif

!! start  DATE_AND_TIME  Timing
      call date_and_time (cdate, ctime, czone, t7)

!!!!!!!!!!!!!!!!!!!!!
!  WORK to TIME 
      do i = 1, length
         a(i) = 1.0d0/dble(i)
      enddo
      sumof = 0.0d0
      do ntimes = 1, 100
         a(ntimes) = a(ntimes) + a(length - ntimes)
         do i = 1, length
            sumof = sumof + a(i)
         enddo
      enddo
!!!!!!!!!!!!!!!!!!!!!

!! end  SECNDS  Timing
      dsecnds = secnds(t1)

!! end  SYSTEM_CLOCK  Timing
      call system_clock (icount, irate)
      dsysclk = ( dble(icount) / dble(irate) ) - t2

!! end  CPU_TIME  Timing
      call cpu_time(t4)
      dcputim = t4 - t3

#ifdef USEMPI
!! end  MPI_WTIME  Timing
      t6 = mpi_wtime()
      mpitim = t6 - t5
#endif

!! end  DATE_AND_TIME  Timing
      call date_and_time (cdate, ctime, czone, t8)
      dandttim = real(60*(t8(6) - t7(6))) + real(t8(7) - t7(7))
     2           + real(t8(8) - t7(8))/1000.0

!!! Write sum of work and each Timers results:
      write(*,*) ' Date: ',cdate
      write(*,*) ' Time: ',ctime
      write(*,*) ' Zone: ',czone
      write(*,*) ' Sum=', sumof
      write(*,*) '  Took ',dsecnds,' Seconds  using SECNDS'
      write(*,*) '  Took ',dsysclk,' Seconds  using SYSTEM_CLOCK'
      write(*,*) '  Took ',dcputim,' Seconds  using CPU_TIME'
#ifdef USEMPI
      write(*,*) '  Took ',mpitim, ' Seconds  using MPI_WTIME'
#endif
      write(*,*) '  Took ',dandttim, ' Seconds  using DATE_AND_TIME'
      write(*,*) ' Began: ',t7
      write(*,*) ' Ended: ',t8
#ifdef USEMPI
      call mpi_finalize(ierr)
#endif

      end

