      module dummy 
      
      contains

      subroutine consumemem(A, mem)
        
      implicit none

      INTEGER mem, rank, pagesz, i
      PARAMETER ( pagesz = 4096/16 )
      complex*16 A(mem)
      do i=1, mem, pagesz 
         if (abs(sin(A(i))) .EQ. 5.0 ) then
            print*,"A(i)= ",A(i) 
         endif 
         A(i) = A(i) + cos(sqrt(5.0d0+0))
      enddo 

      end subroutine consumemem

      end module dummy  
