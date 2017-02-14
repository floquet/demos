
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <mpi.h>

int main(int argc, char * argv[]) {

  int namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME + 1];
  int myproc, numprocs;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myproc);

  MPI_Get_processor_name(processor_name, &namelen);
  processor_name[namelen] = '\0';

  fprintf(stdout,"Hello from proc %i on node %s\n",
	  myproc, processor_name);

  MPI_Finalize();
  return 0;
}
