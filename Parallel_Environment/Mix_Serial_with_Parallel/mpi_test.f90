      program  mpitest
      include 'mpif.h'
      integer ierr, rank, comm_size
      integer i, j
      real*8 z(524288)
      real*8 summ
 
      call MPI_Init(ierr)
      call MPI_Comm_rank (MPI_COMM_WORLD, rank, ierr)
      call MPI_Comm_size (MPI_COMM_WORLD, comm_size, ierr)
 
      do i = 1, 524288
        z(i) = dble(i + rank + comm_size)
      enddo

      do j=1,25000

        summ = 0.0d0

        do i=1,524288

          summ = summ + dble(j + rank)  + z(i) / 500000.0d0

        enddo

        if (mod(j,500).eq.0) write(*,200) j,summ,rank

      enddo

200   format(' j and sum ',i5,f16.2, ' by process ',i3)

      call system(' hostname ') 

90    call MPI_FINALIZE(ierr)
      end
