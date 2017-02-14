      subroutine so3(rank,size,sum) 
      implicit none
      include 'mpif.h'
      real*8 sum
       integer*4 ierr, rank, size


         sum = dfloat(rank+size) * 3.14_8 * 3.0_8

             call MPI_Barrier(MPI_COMM_WORLD, ierr)


      return
      end
