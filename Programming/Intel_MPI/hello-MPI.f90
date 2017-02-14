program mympi
    use mpi

    integer ierror
    integer rank
    integer total

    call MPI_INIT(ierror)
 
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,total,ierror)

    PRINT *, rank, " OF ", total

    call MPI_FINALIZE(ierror)

end program

