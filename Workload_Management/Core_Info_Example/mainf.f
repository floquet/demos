      program  mpitest
      include 'mpif.h'
      integer ierr

      call MPI_Init(ierr)
      call print_core_info(MPI_COMM_WORLD)
      call MPI_FINALIZE(ierr)
      end
