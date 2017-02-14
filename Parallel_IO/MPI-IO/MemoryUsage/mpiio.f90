
subroutine decompose(num_elements, num_slices, my_slice, &
     start_elem, end_elem)
  ! my_slice has values 0 to (num_slices - 1)
  ! element index begins at 1
  ! indices for a slice include end_elem, unlike C++ standard template
  implicit none
  integer, intent(in) :: num_elements, num_slices, my_slice
  integer, intent(out) :: start_elem, end_elem

  integer :: deficit, nlocal
    
  deficit = mod(num_elements, num_slices)
  nlocal = num_elements / num_slices
  start_elem = my_slice*nlocal + min(my_slice, deficit) + 1
  if(my_slice .LT. deficit) then
     nlocal = nlocal + 1
  endif
  end_elem = start_elem + nlocal - 1
  if((end_elem .GT. num_elements) .OR. (my_slice .EQ. (num_slices - 1))) then
     end_elem = num_elements
  endif
end subroutine decompose

subroutine calc_coords(overlap, block_size_in, slice_dims, &
     i_total_max, j_total_max, k_total_max, &
     cart_slice, &
     i_global_start, j_global_start, k_global_start, &
     i_global_end, j_global_end, k_global_end, &
     i_global_size, j_global_size, k_global_size, &
     i_globe_start, j_globe_start, k_globe_start, &
     i_globe_end, j_globe_end, k_globe_end, &
     i_globe_size, j_globe_size, k_globe_size)
  implicit none
  integer, intent(in) :: overlap, block_size_in
  integer, intent(in) :: slice_dims(3)
  integer, intent(in) :: cart_slice(3)
  integer, intent(out) :: i_total_max, j_total_max, k_total_max
  integer, intent(out) :: i_global_start, j_global_start, k_global_start
  integer, intent(out) :: i_global_end, j_global_end, k_global_end
  integer, intent(out) :: i_global_size, j_global_size, k_global_size
  integer, intent(out) :: i_globe_start, j_globe_start, k_globe_start
  integer, intent(out) :: i_globe_end, j_globe_end, k_globe_end
  integer, intent(out) :: i_globe_size, j_globe_size, k_globe_size

  ! Array cart_slice values begins at zero
  ! For a parallel program, each process would have a different cart_slice
  ! global_start, global_end, global_size,
  ! globe_start, globe_end and globe_size are
  ! global indices for a given cart_slice

  ! Grid size.
  i_total_max = slice_dims(1)*block_size_in
  j_total_max = slice_dims(2)*block_size_in
  k_total_max = slice_dims(3)*block_size_in

  ! Calculate global coordinates for a particular cart_slice
  call decompose(i_total_max, slice_dims(1), cart_slice(1), &
       i_global_start, i_global_end)
  call decompose(j_total_max, slice_dims(2), cart_slice(2), &
       j_global_start, j_global_end)
  call decompose(k_total_max, slice_dims(3), cart_slice(3), &
       k_global_start, k_global_end)
  i_global_size = i_global_end - i_global_start + 1
  j_global_size = j_global_end - j_global_start + 1
  k_global_size = k_global_end - k_global_start + 1
  ! Use globe (instead of global) for overlap
  i_globe_start = i_global_start - overlap
  j_globe_start = j_global_start - overlap
  k_globe_start = k_global_start - overlap
  if(i_globe_start.LT.1)then
     i_globe_start = 1
  endif
  if(j_globe_start.LT.1)then
     j_globe_start = 1
  endif
  if(k_globe_start.LT.1)then
     k_globe_start = 1
  endif
  i_globe_end = i_global_end + overlap
  j_globe_end = j_global_end + overlap
  k_globe_end = k_global_end + overlap
  if(i_globe_end.GT.i_total_max)then
     i_globe_end = i_total_max
  endif
  if(j_globe_end.GT.j_total_max)then
     j_globe_end = j_total_max
  endif
  if(k_globe_end.GT.k_total_max)then
     k_globe_end = k_total_max
  endif
  i_globe_size = i_globe_end - i_globe_start + 1
  j_globe_size = j_globe_end - j_globe_start + 1
  k_globe_size = k_globe_end - k_globe_start + 1
end subroutine calc_coords

subroutine init_data(lattice, header, &
     i_globe_size, j_globe_size, k_globe_size, &
     myid, header_size)
  implicit none
  integer, intent(in) :: i_globe_size, j_globe_size, k_globe_size
  integer, intent(in) :: myid, header_size
  double precision, dimension(i_globe_size, j_globe_size, k_globe_size), intent(out) :: lattice
  integer, dimension(header_size), intent(out) :: header
  
  integer :: i, j, k

  forall(i=1:i_globe_size,j=1:j_globe_size,k=1:k_globe_size)
     lattice(i,j,k) = (((myid + 1.0)*i)*j)*k
  end forall
  call set_header(header, header_size)
end subroutine init_data

subroutine set_header(header, header_size)
  implicit none
  integer, intent(in) :: header_size
  integer, dimension(header_size), intent(out) :: header
  integer :: m
  forall(m=1:header_size)
     header(m) = mod(m,128)
  end forall
end subroutine set_header

subroutine prog_die(message, myid)
  use iso_fortran_env
  implicit none
  include "mpif.h"
  character(len = *), intent(in) :: message
  integer, intent(in) :: myid
  integer :: err_code, ierr
  logical :: flag
  
  err_code = 1
  
  write(error_unit,'(i0,1x,a)') myid, trim(message)
  
  call mpi_initialized(flag, ierr)
  
  if (flag) then
     call mpi_abort(MPI_COMM_WORLD, err_code, ierr)
  else
     stop 'program stop in prog_die'
  end if
end subroutine prog_die

subroutine globe2local(il, jl, kl, ig, jg, kg, &
     i_globe_start, j_globe_start, k_globe_start)
  ! The local indices are for the lattice with overlap.
  implicit none
  integer, intent(out) :: il, jl, kl
  integer, intent(in) :: ig, jg, kg
  integer, intent(in) :: i_globe_start, j_globe_start, k_globe_start
  il = ig - i_globe_start + 1
  jl = jg - j_globe_start + 1
  kl = kg - k_globe_start + 1
end subroutine globe2local

subroutine define_data_types(subarray_type1, subarray_type2, &
     i_total_max, j_total_max, k_total_max, &
     i_global_size, j_global_size, k_global_size, &
     i_global_start, j_global_start, k_global_start, &
     i_globe_start, j_globe_start, k_globe_start, &
     i_globe_size, j_globe_size, k_globe_size)
  ! Subarray data types for input and output
  implicit none
  include "mpif.h"
  integer, intent(out) :: subarray_type1, subarray_type2
  integer, intent(in) :: i_total_max, j_total_max, k_total_max
  integer, intent(in) :: i_global_size, j_global_size, k_global_size
  integer, intent(in) :: i_global_start, j_global_start, k_global_start
  integer, intent(in) :: i_globe_start, j_globe_start, k_globe_start
  integer, intent(in) :: i_globe_size, j_globe_size, k_globe_size

  integer, dimension(3) :: big_sizes, sub_sizes, starts
  integer :: ierr, ilocal, jlocal, klocal

  !
  ! Points that this block owns, wrt global array (entire grid).
  !
  
  ! Entire grid.
  big_sizes(:) = (/ i_total_max, j_total_max, k_total_max /)
  ! Local lattice, excluding overlap.
  sub_sizes(:) = (/ i_global_size, j_global_size, k_global_size /)
  starts = (/ i_global_start, j_global_start, k_global_start /) - 1
  call mpi_type_create_subarray(3, big_sizes, sub_sizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, subarray_type1, ierr)
  call mpi_type_commit(subarray_type1, ierr)
  
  !
  ! Points that this block owns, wrt local array with overlap
  !
  
  ! Given the start of the local owned (not overlap) vertices
  ! in terms of the global indices, *_global_start,
  ! calculate the local indices in terms of the lattice with overlap.
  ! Variables with global and globe are both global indices,
  ! globe refers to lattice with overlap
  call globe2local(ilocal, jlocal, klocal, &
       i_global_start, j_global_start, k_global_start, &
       i_globe_start, j_globe_start, k_globe_start)
  
  ! local sizes with overlap
  big_sizes(:) = (/ i_globe_size, j_globe_size, k_globe_size /)
  ! Local lattice, excluding overlap.
  sub_sizes(:) = (/ i_global_size, j_global_size, k_global_size /)
  starts = (/ ilocal, jlocal, klocal /) - 1
  call mpi_type_create_subarray(3, big_sizes, sub_sizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, subarray_type2, ierr)
  call mpi_type_commit(subarray_type2, ierr)
  
end subroutine define_data_types

subroutine show_hints(mpiHints)
  use iso_fortran_env
  implicit none
  include "mpif.h"
  integer :: mpiHints
  integer :: i, nkeys, ierr
  logical flag
  character(len = MPI_MAX_INFO_VAL) :: key, value
  call MPI_Info_get_nkeys(mpiHints, nkeys, ierr)
  if(ierr.NE.MPI_SUCCESS) then
     write(error_unit, '(a)') 'Error from MPI_Info_get_nkeys'
     return
  endif
  do i = 0, nkeys - 1
     call MPI_Info_get_nthkey(mpiHints, i, key, ierr)
     if(ierr.NE.MPI_SUCCESS) then
        write(error_unit, '(a)') 'Error from MPI_Info_get_nthkey'
        return
     endif
     call MPI_Info_get(mpiHints, key, MPI_MAX_INFO_VAL-1, &
                       value, flag, ierr)
     if(ierr.NE.MPI_SUCCESS) then
        write(error_unit, '(a)') 'Error from MPI_Info_get'
        return
     endif
     write(error_unit, '(3a)') TRIM(key), ' = ', TRIM(value)
  enddo
end subroutine show_hints

program test_mpiio
  use iso_fortran_env
  implicit none
  include "mpif.h"

  integer, parameter :: block_size = 200
  integer, parameter :: header_size = 160

  integer :: slice_dims(3), cart_slice(3)
  integer :: i_total_max, j_total_max, k_total_max
  integer :: i_global_start, j_global_start, k_global_start
  integer :: i_global_end, j_global_end, k_global_end
  integer :: i_global_size, j_global_size, k_global_size
  integer :: i_globe_start, j_globe_start, k_globe_start
  integer :: i_globe_end, j_globe_end, k_globe_end
  integer :: i_globe_size, j_globe_size, k_globe_size

  integer :: num_elements, num_slices, total_slices
  integer :: my_slice
  integer :: start_elem, end_elem
  logical :: debug_decompose

  integer :: overlap
  integer :: block_size_in
  logical :: debug_calc_coords
  logical :: oops

  double precision, allocatable, dimension(:,:,:) :: lattice
  integer, allocatable, dimension(:) :: header
  integer :: header_type_size
  integer :: ntimes, n

  integer :: myid, nprocs, ierr

  logical :: cart_reorder, cart_isperiodic
  integer :: dimensionality
  integer :: my_communicator

  character(len=255) :: filename, string
  integer :: nargs
  integer(kind=mpi_offset_kind) :: disp
  integer :: mpiHints
  logical :: start_ok

  integer :: subarray_type1, subarray_type2
  integer :: file_handle

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

  call MPI_Type_size(MPI_INTEGER, header_type_size, ierr)
  disp = header_size*header_type_size

  ! Read command-line parameters.
  if(myid.EQ.0) then
     start_ok = .TRUE.
     nargs = command_argument_count()
     if(nargs.NE.2) then
        start_ok = .FALSE.
     else
        call get_command_argument(1,string)
        read(string,*) slice_dims(1), slice_dims(2), slice_dims(3)
        call get_command_argument(2,filename)
     endif
     if(.NOT.start_ok)then
        write(*,*) '  Error, start not OK'
     endif
  endif ! if(myid.EQ.0)
  call MPI_Bcast(start_ok, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  if(.NOT.start_ok) then
     call MPI_Finalize(ierr)
  endif
  call MPI_Bcast(slice_dims, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(filename, 255, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

  ! below, test decompose
  if(myid.EQ.0) then
     debug_decompose = .FALSE.
     if(debug_decompose) then
        do num_elements = 1, 10
           do num_slices = 1, num_elements
              write(*,*) ' num_elements = ', num_elements, ', &
                   num_slices = ', num_slices
              do my_slice = 0, (num_slices - 1)
                 call decompose(num_elements, num_slices, my_slice, &
                      start_elem, end_elem)
                 write(*,*) my_slice, ') ', start_elem, end_elem
              enddo
           enddo
        enddo
     endif ! if(debug_decompose)
  endif
  ! above, test decompose

  overlap = 3

  total_slices = slice_dims(1)*slice_dims(2)*slice_dims(3)
  if (myid.EQ.0) then
     if(total_slices.NE.nprocs) then
        call prog_die('>>> program requires number_tasks = isplit*jsplit*ksplit', myid)
     endif
     write(error_unit,'(4(a,i0))') '>>> MPI tasks: ', nprocs, ' &
          split1: ', slice_dims(1), &
          ' split2: ', slice_dims(2), ' split3: ', slice_dims(3)
  end if

  cart_reorder = .TRUE.
  cart_isperiodic = .FALSE.
  dimensionality = 3
  call MPI_Cart_create(MPI_COMM_WORLD, dimensionality, slice_dims, &
       cart_isperiodic, cart_reorder, my_communicator, ierr)
  ! Note that myid is redefined.
  call MPI_Comm_rank(my_communicator, myid, ierr)
  call MPI_Cart_get(my_communicator, dimensionality, slice_dims, &
       cart_isperiodic, cart_slice, ierr)

  block_size_in = block_size
  call calc_coords(overlap, block_size_in, slice_dims, &
       i_total_max, j_total_max, k_total_max, &
       cart_slice, &
       i_global_start, j_global_start, k_global_start, &
       i_global_end, j_global_end, k_global_end, &
       i_global_size, j_global_size, k_global_size, &
       i_globe_start, j_globe_start, k_globe_start, &
       i_globe_end, j_globe_end, k_globe_end, &
       i_globe_size, j_globe_size, k_globe_size)
  
  allocate(lattice(i_globe_size, j_globe_size, k_globe_size))
  allocate(header(header_size))

  ! Note that data at overlap are not consistent.
  call init_data(lattice, header, &
     i_globe_size, j_globe_size, k_globe_size, &
     myid, header_size)
  call define_data_types(subarray_type1, subarray_type2, &
       i_total_max, j_total_max, k_total_max, &
       i_global_size, j_global_size, k_global_size, &
       i_global_start, j_global_start, k_global_start, &
       i_globe_start, j_globe_start, k_globe_start, &
       i_globe_size, j_globe_size, k_globe_size)

  ! below, test sublattices
  debug_calc_coords = .FALSE.
  if(debug_calc_coords) then
     oops = .FALSE.
     do n = 0, (nprocs - 1)
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        if(myid.EQ.n) then
!           write(*,'(1x,i6,2x,3(i4,2x))') myid, &
!                i_global_size, j_global_size, k_global_size
           if(  (i_global_size.NE.block_size).OR.&
                (j_global_size.NE.block_size).OR.&
                (k_global_size.NE.block_size) ) then
              write(*,*) ' oops'
              write(*,'(1x,4(i4,2x))') myid, &
                   i_global_size, j_global_size, k_global_size
              oops = .TRUE.
           endif
           if(oops)then
              call prog_die('>>> error in testing sublattices', myid)
           endif
        endif ! if(myid.EQ.n) then
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
     enddo
  endif ! if(debug_calc_coords)
  ! above, test sublattices
  
  ntimes = 10
  if (myid.EQ.0) then
     write(error_unit,'(a)') '>>> writing'
  endif
  if (myid.EQ.0) then
     call MPI_File_delete(filename, MPI_INFO_NULL, ierr)
  endif
  call MPI_File_open(my_communicator, filename, &
       MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, file_handle, ierr)
  if (myid.EQ.0) then
     call MPI_File_get_info(file_handle, mpiHints, ierr)
     if (ierr.NE.mpi_success) call prog_die('>>> MPI_File_get_info failed', myid)
     write(error_unit,'(a)') 'hints passed to MPI_File_open() {'
     call show_hints(mpiHints)
     write(error_unit,'(a)') '}'
  endif

  if (myid.EQ.0) then
     call mpi_file_write(file_handle, header, header_size, &
          MPI_INTEGER, MPI_STATUSES_IGNORE, ierr)
     if (ierr.NE.MPI_SUCCESS) call prog_die('>>> header write failed', myid)
  endif

  ! Data type subarray_type1 is a subdomain of the overall grid,
  ! not including overlap
  call mpi_file_set_view(file_handle, disp, MPI_DOUBLE_PRECISION, &
       subarray_type1, 'NATIVE', MPI_INFO_NULL, ierr)

  do n = 1, ntimes

     ! Data type subarray_type2 is the vertices owned by the subdomain
     ! within the lattice that includes overlap.
     call mpi_file_write_all(file_handle, lattice, 1, &
          subarray_type2, MPI_STATUSES_IGNORE, ierr)
     if (ierr.NE.mpi_success) call prog_die('>>> write failed', myid)
     if (myid.EQ.0) then
        write(error_unit,'(2(a,i0))') &
             '>>> finished writing step ', n, ' of ', ntimes
     endif
  enddo

  call MPI_File_close(file_handle, ierr)

  if (myid.EQ.0) then
     write(error_unit,'(a)') '>>> done'
  endif

  deallocate(lattice)
  deallocate(header)

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  call MPI_Finalize(ierr)
  stop
end program test_mpiio

