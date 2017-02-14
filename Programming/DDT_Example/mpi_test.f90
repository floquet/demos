      program  mpitest
      implicit none
      include 'mpif.h'
      REAL sum, asum, xi
 
      integer ierr, rank, i, j, ist, iend, size
 
      call MPI_Init(ierr)
      call MPI_Comm_rank (MPI_COMM_WORLD, rank, ierr)
      call MPI_Comm_size (MPI_COMM_WORLD, size, ierr)
 
        iend = 5000000 / size
        ist =  rank * iend
        iend = ist + iend
        ist  = ist + 1
 
        do j=1,800
          xi = 1.0
 
          sum = 0.0

          if( j .eq. 191 ) xi = 0.0  ! will generate divide by zero FP exception
          do i=ist,iend
 
            sum = sum + (sqrt(float(i) * float(i))) /(5000000.0*xi) + float(j)
 
          enddo
 
          call MPI_AllReduce(sum,asum,1,MPI_REAL,MPI_SUM, MPI_COMM_WORLD,ierr)
 
          write(*,*)' j and sum ',j,asum
 
        enddo
 
90      call MPI_FINALIZE(ierr)
      STOP
      end
