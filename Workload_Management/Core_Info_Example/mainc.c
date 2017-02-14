#include <mpi.h>
#include <omp.h>

extern void print_core_info(MPI_Comm);

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    print_core_info(MPI_COMM_WORLD);
    MPI_Finalize();
    return(0);
}
