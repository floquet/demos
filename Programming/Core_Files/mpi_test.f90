      program  mpitest
      implicit none
      include 'mpif.h'
      REAL tsum, asum, y
 
      integer ierr, rank, i, j, ist, iend, isize
 
      ierr = 0
      rank = 0
      isize = 0
      call MPI_Init(ierr)
      call MPI_Comm_rank (MPI_COMM_WORLD, rank, ierr)
      call MPI_Comm_size (MPI_COMM_WORLD, isize, ierr)
 
      iend = 8000000 / isize
 
      ist =  rank * iend
      iend = ist + iend
      ist  = ist + 1

!! when j=47  then y=0 and divide by zero error occurs 

      do j=1,50
        y = real(j)
        if (j.eq.47) y = 0.0
 
        tsum = 0.0
 
        do i=ist,iend
  
          tsum = tsum + (sqrt(real(i) * real(i))) /(y*50000.0)
 
        enddo
 
        call MPI_AllReduce(tsum,asum,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
 
        write(*,*)' j and sum ',j,asum
 
      enddo
 
90    call MPI_FINALIZE(ierr)
      STOP
      end
