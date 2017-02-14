      program  serial_test
      integer i, j
      real*8 z(524288)
      real*8 summ

      do i = 1, 524288
        z(i) = dble(i)
      enddo

      do j=1,25000

        summ = 0.0d0

        do i=1,524288

          summ = summ + dble(j)  + z(i) / 5000000.0d0

        enddo

        if (mod(j,500).eq.0) write(*,200) j,summ

      enddo
       
      call system(' hostname ')

200   format(' j and sum ',i5,f16.2)

      call sleep(120)
      end
