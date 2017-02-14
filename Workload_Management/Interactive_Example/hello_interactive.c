#include <mpi.h>
#include <omp.h>

extern void print_hello_info(MPI_Comm);

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    print_hello_info(MPI_COMM_WORLD);
    MPI_Finalize();
    return(0);
}
