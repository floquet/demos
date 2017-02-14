
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>

#include <mpi.h>
#include "measure_memory.h"


int main(int argc, char * argv[]) {
  int myproc, numprocs, which_proc, i;
  size_t sum;
  int buffer_size, buffer_size_base;
  int * buffer;
  int pid_rsm_retval;
  pid_t my_pid;
  int status;
  int64_t memory_used, mem_before_malloc, mem_after_malloc;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myproc);

  status = 0;
  buffer_size_base = 1024*512;
  my_pid = getpid();

  for(which_proc = 0; which_proc < numprocs; which_proc++) {
    buffer_size = myproc*buffer_size_base;
    if(which_proc == myproc) {

      mem_before_malloc = 0;
      pid_rsm_retval = pid_rsm(&my_pid, &memory_used);
      if(pid_rsm_retval != 0) {
        fprintf(stderr,"Error returned from pid_rsm\n");
        fprintf(stderr,"line %d file %s\n", __LINE__, __FILE__);
        status = 1;
      }
      else { mem_before_malloc = memory_used; }
      if((myproc % 4) == 0) {
	fprintf(stdout," Proc %d, mem before malloc %ld\n",
		myproc, mem_before_malloc);
      }

      buffer = (int *) malloc(buffer_size*sizeof(int));

      // Malloc'ed memory not allocated until used.
      for(i = 0; i < buffer_size; i++) {
	buffer[i] = i;
      }
      // Avoid initialization being optimized away by using result.
      sum = 0;
      for(i = 0; i < buffer_size; i++) {
	sum += buffer[i];
      }
      if(sum == 3) { fprintf(stderr,"sum is three\n"); }

      mem_after_malloc = 0;
      pid_rsm_retval = pid_rsm(&my_pid, &memory_used);
      if(pid_rsm_retval != 0) {
        fprintf(stderr,"Error returned from pid_rsm\n");
        fprintf(stderr,"line %d file %s\n", __LINE__, __FILE__);
        status = 1;
      }
      else { mem_after_malloc = memory_used; }
      if((myproc % 4) == 0) {
	fprintf(stdout," Proc %d, mem after malloc  %ld\n",
		myproc, mem_after_malloc);
      }

    } // if(which_proc == myproc)
    MPI_Barrier(MPI_COMM_WORLD);
  } // for(which_proc = 0; which_proc < numprocs; which_proc++)

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return status;
}
