      program  mpitest
      include 'mpif.h'
      integer rank, ssize, ierr

      call MPI_Init (ierr)
      call MPI_Comm_rank (MPI_COMM_WORLD, rank, ierr)
      call MPI_Comm_size (MPI_COMM_WORLD, ssize, ierr)

      if (rank.eq.1) then
        write(*,*)" Hello I am compute node 1 from job 1"
        call system (" hostname ")
      endif

      call sleep(45)  ! waste some time  
 
90    call MPI_FINALIZE (ierr)

      end

