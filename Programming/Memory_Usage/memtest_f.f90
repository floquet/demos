      program memtest
      implicit none
      include 'mpif.h'
      integer myproc, numprocs, which_proc, i, ierr
      integer buffer_size, buffer_size_base
      integer my_pid
      integer getpid, pid_rsm, pid_rsm_retval
      integer(8) memory_used, mem_before_malloc, mem_after_malloc
      integer, dimension(:), allocatable :: buffer
      integer(8) sum

      call MPI_Init(ierr)
      call MPI_Comm_rank (MPI_COMM_WORLD, myproc, ierr)
      call MPI_Comm_size (MPI_COMM_WORLD, numprocs, ierr)

      buffer_size_base = 1024*512
      my_pid = getpid()

      do which_proc = 0, (numprocs - 1)
         buffer_size = myproc*buffer_size_base
         if(which_proc .eq. myproc) then

            mem_before_malloc = 0;
            pid_rsm_retval = pid_rsm(my_pid, memory_used)
            if(pid_rsm_retval .ne. 0) then
               write(*,*) ' Error returned from pid_rsm'
            else
               mem_before_malloc = memory_used
            endif
            if(mod(myproc, 4) .eq. 0) then
               write(*,*) ' Proc ',myproc,', &
                    mem before_malloc ', mem_before_malloc
            endif

            allocate(buffer(buffer_size))
            ! Allocate'd memory not allocated until used.
            do i = 1, buffer_size
               buffer(i) = i
            enddo
            ! Avoid initialization being optimized away by using result.
            sum = 0
            do i = 1, buffer_size
               sum = sum + buffer(i)
            enddo
            if(sum .eq. 3) then
               write(*,*) 'sum is three'
            endif

            mem_after_malloc = 0;
            pid_rsm_retval = pid_rsm(my_pid, memory_used)
            if(pid_rsm_retval .ne. 0) then
               write(*,*) ' Error returned from pid_rsm'
            else
               mem_after_malloc = memory_used
            endif
            if(mod(myproc, 4) .eq. 0) then
               write(*,*) ' Proc ',myproc,', &
                    mem after_malloc ', mem_after_malloc
            endif

         endif
         
         call MPI_Barrier(MPI_COMM_WORLD, ierr)
      enddo

      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call MPI_Finalize(ierr)
      stop
      end
