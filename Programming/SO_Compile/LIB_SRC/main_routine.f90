      program  testnode
      implicit none
       include 'mpif.h'
        integer rank, size, ierr
         real*8 sum


     
      call MPI_Init(ierr)
      call MPI_Comm_rank (MPI_COMM_WORLD, rank, ierr)
      call MPI_Comm_size (MPI_COMM_WORLD, size, ierr)

!  call SO one
          call so1(rank,size,sum)
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
              call sleep(2)
                write(6,*) " (rank+size)x1xPie=  ",rank, sum

!  call SO two
          call so2(rank,size,sum)
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
              call sleep(2)
                write(6,*) " (rank+size)x2xPie=  ",rank, sum

!  call SO three
          call so3(rank,size,sum)
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
              call sleep(2)
                write(6,*) " (rank+size)x3xPie=  ",rank, sum



      call MPI_FINALIZE(ierr)
      end
