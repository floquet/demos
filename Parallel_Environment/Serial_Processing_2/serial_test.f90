      program  serial_test
      integer i, j
      real*8 z(524288)
      real*8 summ
 
      do i = 1, 524288
        z(i) = dble(i)
      enddo

      do j=1,3000

        summ = 0.0d0

        do i=1, 524288

          summ = summ + dble(j)  + z(i) / 5000000.0d0

        enddo

        write(*,200) j,summ

      enddo

200   format(' j and sum ',i4,f16.2)

      STOP
      end
